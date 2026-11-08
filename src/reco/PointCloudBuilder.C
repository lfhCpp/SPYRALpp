#include "PointCloudBuilder.h"
#include "NumpyArgsort.h"

#include <cmath>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace rk4 {

PointCloudBuilder::~PointCloudBuilder() = default;

bool PointCloudBuilder::Init(const DetectorConfig& cfg, const PadPlane* pads,
                             const PeakFinderParameters& peak_params,
                             double baseline_window_scale) {
    pads_ = nullptr;
    wave_.reset();
    if (pads == nullptr) return false;

    cfg_ = cfg;
    pads_ = pads;
    peak_params_ = peak_params;

    // 示例数据的波形长度固定为 512；后续可从配置读取
    trace_length_ = 512;
    wave_.reset(WaveProcessor::create(trace_length_, baseline_window_scale));
    if (!wave_) {
        std::fprintf(stderr, "[PointCloudBuilder] 波形处理器初始化失败\n");
        return false;
    }
    wave_->setIntegerInput(true);

    input_.assign(static_cast<size_t>(trace_length_), 0.0);
    filtered_.assign(static_cast<size_t>(trace_length_), 0.0);
    return true;
}

void PointCloudBuilder::ConfigureTemplateFit(bool enable, const PulseShapeConfig& rc,
                                            const TemplateFitConfig& tc) {
    template_fit_ = enable;
    response_cfg_ = rc;
    template_cfg_ = tc;
    if (enable) response_.Configure(rc);
}

bool PointCloudBuilder::Build(const TraceEvent& ev, PointCloud& cloud) {
    cloud=PointCloud{};
    accepted_traces_=0;total_peaks_=0;
    cloud.event = ev.number;
    cloud.metadata = ev.metadata;
    if(cloud.metadata.orig_event<0)cloud.metadata.orig_event=ev.number;
    if(!ev.ic_trace.empty()){
        cloud.metadata.ic_amplitude=-1;cloud.metadata.ic_centroid=-1;
        cloud.metadata.ic_integral=-1;cloud.metadata.ic_multiplicity=-1;
        if(ev.ic_trace.size()!=2048)return false;
        if(!frib_wave_){
            frib_wave_=WaveProcessor::createUnique(2048,frib_config_.baseline_scale);
            if(!frib_wave_)return false;
            frib_wave_->setIntegerInput(true);
        }
        std::vector<double> raw(ev.ic_trace.begin(),ev.ic_trace.end()),filtered(2048);
        if(!frib_wave_->process(raw.data(),filtered.data()))return false;
        for(auto& value:filtered)value=std::trunc(value);
        PeakFinder finder(frib_config_.peaks);finder.setTrace(filtered);
        const auto peaks=finder.findPeaksWithDither([]{return 0.;});
        const Peak* trigger=nullptr;int multiplicity=0;
        for(const auto& peak:peaks)if(peak.centroid>frib_config_.delay){
            ++multiplicity;if(!trigger)trigger=&peak;
        }
        if(trigger&&multiplicity<=frib_config_.max_multiplicity){
            cloud.metadata.ic_amplitude=trigger->amplitude;
            cloud.metadata.ic_centroid=trigger->centroid;
            cloud.metadata.ic_integral=trigger->integral;
            cloud.metadata.ic_multiplicity=multiplicity;
        }
    }
    cloud.points.clear();
    cloud.dither_draws = 0;
    // Keep the same run-wide NumPy dithering stream as the reference path.
    // Template fitting must not change the stochastic z coordinate used by
    // clustering; otherwise tiny charge-model differences can reshuffle tracks.
    cloud.pending_dither = defer_dither_;
    accepted_traces_ = 0;
    total_peaks_ = 0;

    if (!wave_ || pads_ == nullptr) return false;

    const int num_pads = pads_->NumPads();

    // Workers defer the run-wide NumPy random stream until results are ordered.

    for (const auto& trace : ev.traces) {
        if (static_cast<int>(trace.samples.size()) != trace_length_) continue;
        if (trace.pad == -1) continue;

        for (int i = 0; i < trace_length_; ++i) {
            input_[static_cast<size_t>(i)] = static_cast<double>(trace.samples[i]);
        }

        if (!wave_->process(input_.data(), filtered_.data())) continue;
        ++accepted_traces_;
        for (double& value : filtered_) value = std::trunc(value);

        PeakFinder finder(peak_params_);
        finder.setTrace(filtered_);
        const size_t first_draw = cloud.dither_draws;
        const std::vector<Peak> peaks = finder.findPeaksWithDither([&] {
            return cloud.pending_dither ? 0.0 : rng_.Random();
        });
        cloud.dither_draws += finder.DitherDraws();
        // Python GetTrace consumes draws before PointCloud filters these pads.
        if (trace.pad < 0 || trace.pad >= num_pads || pads_->IsBeamPad(trace.pad)) continue;
        if (peaks.empty()) continue;

        // ---- 可选：GET 响应模板拟合（含饱和修复）----
        if (template_fit_) {
            // Match the standard point-cloud acceptance exactly. The response
            // model is not channel calibrated, so it must not invent peaks or
            // remove otherwise valid spatial measurements.
            if (peaks.size() > 5) continue;
            const SaturationInfo sat =
                DetectSaturation(input_, response_cfg_.sat_adc, 2);
            if (sat.saturated) ++saturated_traces_;

            fit_candidates_.clear();
            for (const Peak& p : peaks) fit_candidates_.push_back(p.centroid);
            // The analytic response is not channel calibrated. Apply its
            // charge correction only where clipping makes the measured charge
            // unusable; preserve baseline charge for ordinary pulses.
            const std::vector<FittedPeak> fitted = sat.saturated
                ? FitPeaksWithTemplate(filtered_, response_, sat, fit_candidates_, template_cfg_)
                : std::vector<FittedPeak>{};

            // Associate template fits with their original peak candidates.
            // Every accepted standard peak remains one point; failed, merged,
            // or sub-threshold fits fall back to the measured peak values.
            std::vector<const FittedPeak*> match(peaks.size(), nullptr);
            std::vector<bool> used(fitted.size(), false);
            for (size_t i = 0; i < peaks.size(); ++i) {
                double best = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < fitted.size(); ++j) {
                    if (used[j]) continue;
                    const double distance = std::fabs(fitted[j].time - peaks[i].centroid);
                    if (distance < best) { best = distance; match[i] = &fitted[j]; }
                }
                if (match[i] && best <= 3.0) {
                    used[static_cast<size_t>(match[i] - fitted.data())] = true;
                } else {
                    match[i] = nullptr;
                }
            }

            const double offset = pads_->TimeOffset(trace.pad);
            const double x = pads_->X(trace.pad);
            const double y = pads_->Y(trace.pad);
            if (peak_sink_ != nullptr) {
                for (size_t i = 0; i < peaks.size(); ++i) {
                    TracePeakInfo info;
                    info.pad = trace.pad;
                    info.index = static_cast<int>(i);
                    info.dither_index = first_draw + peaks[i].dither_index;
                    info.npeaks = static_cast<int>(peaks.size());
                    info.time = peaks[i].centroid + offset;
                    info.amplitude = match[i] ? match[i]->amplitude : peaks[i].amplitude;
                    info.charge = match[i] ? match[i]->charge : peaks[i].integral;
                    info.chi2 = match[i] ? match[i]->chi2 : 0.0;
                    info.npoints = match[i] ? match[i]->n_points : 0;
                    info.saturated = sat.saturated ? 1 : 0;
                    peak_sink_->push_back(info);
                }
            }

            for (size_t i = 0; i < peaks.size(); ++i) {
                const Peak& peak = peaks[i];
                const FittedPeak* f = match[i];
                const double tb = peak.centroid + offset;
                CloudPoint p;
                p.pad = trace.pad;
                p.tb = tb;
                p.x = x;
                p.y = y;
                p.z = (cfg_.windows_edge - tb) /
                      (cfg_.windows_edge - cfg_.micromegas_edge) * cfg_.length_m * 1000.0;
                p.dither_index = first_draw + peak.dither_index;
                p.amplitude = f ? f->amplitude : peak.amplitude;
                p.integral = f ? f->charge : peak.integral;
                p.scale = pads_->Scale(trace.pad);
                cloud.points.push_back(p);
                ++total_peaks_;
            }
            continue;
        }

        if (peaks.size() > 5) continue;

        const double offset = pads_->TimeOffset(trace.pad);
        const double x = pads_->X(trace.pad);
        const double y = pads_->Y(trace.pad);

        if (peak_sink_ != nullptr) {
            const double off = pads_->TimeOffset(trace.pad);
            for (size_t i = 0; i < peaks.size(); ++i) {
                TracePeakInfo info;
                info.pad = trace.pad;
                info.dither_index = first_draw + peaks[i].dither_index;
                info.index = static_cast<int>(i);
                info.npeaks = static_cast<int>(peaks.size());
                info.time = peaks[i].centroid + off;
                info.amplitude = peaks[i].amplitude;
                info.charge = peaks[i].integral;
                info.saturated = 0;
                peak_sink_->push_back(info);
            }
        }

        for (const auto& peak : peaks) {
            const double tb = peak.centroid + offset;
            CloudPoint p;
            p.dither_index = first_draw + peak.dither_index;
            p.pad = trace.pad;
            p.tb = tb;
            p.x = x;
            p.y = y;
            p.z = (cfg_.windows_edge - tb) /
                  (cfg_.windows_edge - cfg_.micromegas_edge) * cfg_.length_m * 1000.0;
            p.amplitude = peak.amplitude;
            p.integral = peak.integral;
            p.scale = pads_->Scale(trace.pad);
            cloud.points.push_back(p);
            ++total_peaks_;
        }
    }

    if (!cloud.pending_dither) {
        std::vector<double> z;
        for (const auto& p : cloud.points) z.push_back(p.z);
        const auto order = NumpyArgsort(std::move(z));
        const auto unsorted = cloud.points;
        for (size_t i=0;i<order.size();++i) cloud.points[i] = unsorted[order[i]];
    }

    return true;
}

void PointCloudBuilder::FinalizeDither(PointCloud& cloud, NumpyPcg64& rng,
                                      const DetectorConfig& cfg, std::vector<TracePeakInfo>* peaks) {
    if (!cloud.pending_dither) return;
    std::vector<double> draws(cloud.dither_draws);
    for (double& value : draws) value = rng.Random();
    if (peaks) for (auto& peak : *peaks) peak.time += draws.at(peak.dither_index);
    std::vector<double> z;
    for (auto& p : cloud.points) {
        p.tb += draws.at(p.dither_index);
        p.z = (cfg.windows_edge - p.tb) / (cfg.windows_edge - cfg.micromegas_edge)
              * cfg.length_m * 1000.0;
        z.push_back(p.z);
    }
    const auto order = NumpyArgsort(std::move(z));
    const auto unsorted = cloud.points;
    for (size_t i=0;i<order.size();++i) cloud.points[i] = unsorted[order[i]];
    cloud.pending_dither = false;
}

} // namespace rk4
