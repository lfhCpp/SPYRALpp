// Independent stage engine. Python launcher only handles configuration/Parquet.
#include "RecoConfig.h"
#include "io/WorkspaceIO.h"
#include "io/PidStateIO.h"
#include "reco/Clusterer.h"
#include "reco/KalmanFitter.h"
#include "reco/TrackFitter.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "ROOT/RConfig.hxx"
#include "TROOT.h"
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <thread>
#include <stdexcept>
#include <utility>
using namespace rk4;
namespace {
using Clock=std::chrono::steady_clock;
void Need(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
std::string Json(const std::string& s){std::string r="\"";for(unsigned char c:s){if(c=='"'||c=='\\'){r+='\\';r+=c;}else if(c=='\n')r+="\\n";else if(c=='\r')r+="\\r";else if(c=='\t')r+="\\t";else if(c<32)throw std::runtime_error("Control character in config");else r+=c;}return r+'"';}
Config Load(const std::string& p){Config c;Need(c.Load(p),"Cannot read config "+p);return c;}
class StageProgress{
public:
    StageProgress(std::string stage,size_t total,std::string unit)
        :stage_(std::move(stage)),total_(total),unit_(std::move(unit)),
         interval_(std::max<size_t>(1,total/50+(total%50!=0))){}
    void Step(){
        if(total_==0)return;
        std::lock_guard<std::mutex> lock(mutex_);
        ++done_;
        if(done_%interval_!=0&&done_!=total_)return;
        std::cout<<stage_<<" progress: "<<done_<<"/"<<total_<<' '<<unit_
                 <<" ("<<(done_*100/total_)<<"%)\n"<<std::flush;
    }
private:
    std::string stage_;
    size_t total_=0,done_=0;
    std::string unit_;
    size_t interval_=1;
    std::mutex mutex_;
};
struct ProgressGuard{StageProgress& progress;~ProgressGuard(){progress.Step();}};
std::string LowerTrim(std::string value){
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    value=value.substr(first,last-first+1);
    for(char& c:value)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}
std::string IsotopeLabel(int A,int Z){
    static const char* symbols[]={"","H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca"};
    if(Z>0&&Z<static_cast<int>(sizeof(symbols)/sizeof(symbols[0])))
        return std::to_string(A)+symbols[Z];
    return std::to_string(A)+"Z"+std::to_string(Z);
}
std::vector<std::pair<int,int>> ParseKalmanHypotheses(const std::string& value,int defaultA,int defaultZ){
    if(LowerTrim(value).empty())return {{defaultA,defaultZ}};
    std::vector<std::pair<int,int>> result;
    std::set<std::pair<int,int>> seen;
    size_t begin=0;
    while(begin<=value.size()){
        const size_t end=value.find(',',begin);
        const std::string token=LowerTrim(value.substr(begin,end==std::string::npos?std::string::npos:end-begin));
        Need(!token.empty(),"solver.kalman_hypotheses contains an empty entry");
        int A=0,Z=0;
        if(token=="p"||token=="proton"||token=="1h"){A=1;Z=1;}
        else if(token=="d"||token=="deuteron"||token=="2h"){A=2;Z=1;}
        else if(token=="t"||token=="triton"||token=="3h"){A=3;Z=1;}
        else if(token=="he3"||token=="3he"){A=3;Z=2;}
        else if(token=="alpha"||token=="4he"){A=4;Z=2;}
        else{
            const size_t colon=token.find(':');
            Need(colon!=std::string::npos,"Unknown Kalman hypothesis '"+token+"' (use p,d,t,he3,alpha or A:Z)");
            try{A=std::stoi(token.substr(0,colon));Z=std::stoi(token.substr(colon+1));}
            catch(...){throw std::runtime_error("Invalid A:Z Kalman hypothesis '"+token+"'");}
        }
        Need(A>0&&Z>0&&Z<=A,"Invalid Kalman isotope A:Z="+std::to_string(A)+":"+std::to_string(Z));
        if(seen.insert({A,Z}).second)result.emplace_back(A,Z);
        if(end==std::string::npos)break;
        begin=end+1;
    }
    Need(!result.empty(),"solver.kalman_hypotheses resolved to no hypotheses");
    return result;
}
template<class Factory>void Parallel(size_t count,const Config& cfg,Factory factory){
    int n=cfg.GetInt("run.threads",1);Need(n>=1,"run.threads must be positive");n=std::min<size_t>(n,std::max<size_t>(count,1));
    std::atomic<size_t> next{0};std::atomic<bool> stop{false};std::exception_ptr error;std::mutex mutex;
    auto worker=[&]{try{auto process=factory();while(!stop){size_t i=next.fetch_add(1);if(i>=count)break;process(i);}}catch(...){std::lock_guard<std::mutex> lock(mutex);if(!error)error=std::current_exception();stop=true;}};
    std::vector<std::thread> pool;for(int i=1;i<n;++i)pool.emplace_back(worker);worker();for(auto& t:pool)t.join();if(error)std::rethrow_exception(error);
}
std::vector<Cell> IDs(long long event,const Cluster& c,const EventMetadata& m){return {event,static_cast<long long>(c.label),static_cast<long long>(c.original_label<0?c.label:c.original_label),m.orig_run,m.orig_event<0?event:m.orig_event};}
void ClusterStage(const std::string& input,const std::string& output,const Config& cfg){
    auto clouds=ReadPointcloudStage(input);ClusterData data;data.info=clouds.info;data.events.resize(clouds.events.size());const auto params=ReadClusterParams(cfg);
    StageProgress progress("Cluster",data.events.size(),"events");
    std::cout<<"Cluster: processing "<<data.events.size()<<" events\n"<<std::flush;
    Parallel(data.events.size(),cfg,[&]{auto clusterer=std::make_shared<Clusterer>();clusterer->SetParams(params);return [&,clusterer](size_t i){auto& dst=data.events[i];const auto& src=clouds.events[i];dst.event=src.event;dst.metadata=src.metadata;clusterer->Run(src,dst.clusters);progress.Step();};});
    size_t clusters=0,points=0;for(auto& e:data.events)for(auto& c:e.clusters){++clusters;points+=c.points.size();}
    WriteClusterStage(output,data);std::cout<<"Cluster: "<<data.events.size()<<" events, "<<clusters<<" clusters, "<<points<<" points\n";
}
void EstimateStage(const std::string& input,const std::string& output,const Config& cfg){
    auto data=ReadClusterStage(input);ScalarTable table;table.info=data.info;table.columns=EstimationColumns();std::vector<std::vector<std::vector<Cell>>> rows(data.events.size());
    const double field=cfg.GetDouble("detector.bfield_t",2.85);
    StageProgress progress("Estimation",data.events.size(),"events");
    std::cout<<"Estimation: processing "<<data.events.size()<<" events\n"<<std::flush;
    Parallel(data.events.size(),cfg,[&]{auto estimator=std::make_shared<Estimator>();ConfigureEstimator(*estimator,cfg);return [&,estimator](size_t i){const auto& e=data.events[i];for(const auto& c:e.clusters){EstimateResult est;if(!estimator->Estimate(c,field,est))continue;
        auto row=IDs(e.event,c,e.metadata);const auto& m=e.metadata;
        for(double value:{m.ic_amplitude,m.ic_centroid,m.ic_integral,m.ic_multiplicity,est.vertex_x,est.vertex_y,est.vertex_z,est.center_x,est.center_y,est.vertex_z,est.polar,est.azimuthal,est.brho,est.dEdx,std::sqrt(est.dEdx),est.charge,est.arclength})row.emplace_back(value);
        row.emplace_back(static_cast<long long>(est.direction));rows[i].push_back(std::move(row));}progress.Step();};});
    for(auto& event:rows)for(auto& row:event)table.rows.push_back(std::move(row));WriteScalarStage(output,"Estimation",table);std::cout<<"Estimation: "<<table.rows.size()<<" accepted trajectories\n";
}
bool Inside(double x,double y,const std::vector<std::pair<double,double>>& polygon){
    if(!std::isfinite(x)||!std::isfinite(y))return false;bool inside=false;
    for(size_t i=0,j=polygon.size()-1;i<polygon.size();j=i++){
        auto [xi,yi]=polygon[i];auto [xj,yj]=polygon[j];
        if((yi>y)!=(yj>y)&&x<(xj-xi)*(y-yi)/(yj-yi)+xi)inside=!inside;
    }return inside;
}
void SolveStage(const std::string& cluster_path,const std::string& estimates_path,const std::string& output,const Config& cfg){
    auto clusters=ReadClusterStage(cluster_path);auto estimates=ReadScalarStage(estimates_path,"Estimation");
    std::map<long long,const ClusterEvent*> events;for(auto& e:clusters.events)Need(events.emplace(e.event,&e).second,"Duplicate cluster event");
    std::vector<std::pair<double,double>> polygon;const int vertices=cfg.GetInt("pid.vertices",0);Need(vertices==0||vertices>=3,"Invalid PID polygon");
    for(int i=0;i<vertices;++i){std::string k="pid.vertex"+std::to_string(i);polygon.emplace_back(cfg.GetDouble(k+".x"),cfg.GetDouble(k+".y"));}
    const std::string pid_state_path=cfg.GetString("output.pid_state_csv","");
    const bool collect_pid_audit=!pid_state_path.empty();
    std::vector<PidCandidateAuditRow> pid_audit;
    std::map<std::pair<long long,long long>,size_t> audit_by_key;
    if(collect_pid_audit)for(const auto& event:clusters.events)for(size_t ci=0;ci<event.clusters.size();++ci){
        const auto& cluster=event.clusters[ci];PidCandidateAuditRow audit;
        audit.event=event.event;audit.cluster_index=static_cast<int>(ci);
        audit.cluster_label=cluster.original_label<0?cluster.label:cluster.original_label;
        audit.n_points=static_cast<int>(cluster.points.size());audit.ic_amplitude=event.metadata.ic_amplitude;
        audit.ic_state=ClassifyIcObservation(audit.ic_amplitude);audit.fit_rejection_stage="estimate_missing";
        audit.pid_state=PidState::PID_UNKNOWN;audit.pid_unknown_reason=PidUnknownReason::NOT_CONFIGURED;
        const auto key=std::make_pair(event.event,static_cast<long long>(ci));
        Need(audit_by_key.emplace(key,pid_audit.size()).second,"Duplicate cluster candidate");
        pid_audit.push_back(std::move(audit));
    }
    std::vector<size_t> estimate_audit_index(collect_pid_audit?estimates.rows.size():0,0);
    const auto write_pid_audit=[&](){if(!pid_state_path.empty())Need(
        WritePidCandidateAuditCsv(pid_state_path,pid_audit),"Cannot write PID state audit: "+pid_state_path);};
    const double ic_min=cfg.GetDouble("solver.ic_min_val",-INFINITY);
    const double ic_max=cfg.GetDouble("solver.ic_max_val",INFINITY);
    const bool ic_gate_configured=cfg.Has("solver.ic_min_val")||cfg.Has("solver.ic_max_val");
    std::vector<size_t> selected;std::map<std::pair<long long,long long>,bool> keys;
    for(size_t i=0;i<estimates.rows.size();++i){
        auto key=std::make_pair(estimates.Integer(i,"event"),estimates.Integer(i,"cluster_index"));Need(keys.emplace(key,true).second,"Duplicate estimate identifier");
        auto found=events.find(key.first);Need(found!=events.end()&&key.second>=0&&static_cast<size_t>(key.second)<found->second->clusters.size(),"Estimation/Cluster identifiers do not match");
        const auto& c=found->second->clusters[key.second];Need(estimates.Integer(i,"cluster_label")==static_cast<long long>(c.original_label<0?c.label:c.original_label),"Estimation/Cluster labels do not match");
        const auto& metadata=found->second->metadata;
        Need((metadata.orig_run<0||estimates.Integer(i,"orig_run")==metadata.orig_run)&&(metadata.orig_event<0||estimates.Integer(i,"orig_event")==metadata.orig_event),"Estimation/Cluster original identifiers do not match");
        PidCandidateAuditRow* audit=nullptr;
        if(collect_pid_audit){const auto audit_it=audit_by_key.find(key);Need(audit_it!=audit_by_key.end(),"Missing PID audit candidate");
            estimate_audit_index[i]=audit_it->second;audit=&pid_audit[audit_it->second];
            audit->event=key.first;audit->cluster_index=static_cast<int>(key.second);
            audit->cluster_label=static_cast<int>(c.original_label<0?c.label:c.original_label);
            audit->n_points=static_cast<int>(c.points.size());audit->estimate_ok=1;
            audit->ic_amplitude=estimates.Number(i,"ic_amplitude");
            audit->ic_state=ClassifyIcObservation(audit->ic_amplitude);
            audit->ic_gate_configured=ic_gate_configured?1:0;
            audit->pid_gate_configured=polygon.empty()?0:1;
            audit->pid_gate_pass=polygon.empty()?-1:0;
            audit->pid_state=PidState::PID_UNKNOWN;
            // The workspace polygon is a selection cut, not a calibrated species classifier.
            audit->pid_unknown_reason=PidUnknownReason::NOT_CONFIGURED;}
        const double ic=estimates.Number(i,"ic_amplitude");
        const bool ic_pass=ic>ic_min&&ic<ic_max;
        if(audit){audit->ic_gate_pass=ic_pass?1:0;audit->fit_rejection_stage="ic_gate";}
        if(!ic_pass)continue;
        if(!polygon.empty()){
            const bool inside=Inside(estimates.Number(i,cfg.GetString("pid.xaxis","sqrt_dEdx")),
                                     estimates.Number(i,cfg.GetString("pid.yaxis","brho")),polygon);
            if(audit)audit->pid_gate_pass=inside?1:0;
            if(!inside){if(audit)audit->fit_rejection_stage="pid_selection_gate";continue;}
        }
        if(audit)audit->fit_rejection_stage="fit_not_run";
        selected.push_back(i);
    }
    ScalarTable table;table.info=clusters.info;table.columns=SolverColumns();
    std::cout<<"InterpSolver: "<<selected.size()<<" / "<<estimates.rows.size()<<" pass IC/PID gates\n";
    if(selected.empty()){write_pid_audit();WriteScalarStage(output,"InterpSolver",table);return;}
    StageProgress fit_progress("InterpSolver",selected.size(),"tracks");
    std::cout<<"InterpSolver: fitting "<<selected.size()<<" tracks\n"<<std::flush;
    std::vector<Element> comp;int n=cfg.GetInt("gas.components",0);
    if(n){for(int i=0;i<n;++i){std::string k="gas.component"+std::to_string(i);Element e;e.a=cfg.GetInt(k+".a");e.z=cfg.GetInt(k+".z");e.s=cfg.GetDouble(k+".s");comp.push_back(e);}}
    else{Element e;e.a=cfg.GetInt("gas.a",1);e.z=cfg.GetInt("gas.z",1);e.s=cfg.GetDouble("gas.stoich",2);comp.push_back(e);}
    GasTarget gas(comp,cfg.GetDouble("gas.pressure_torr",600));const int A=cfg.GetInt("solver.a",1),Z=cfg.GetInt("solver.z",1);
    const double mass=MassData::GetNuclearMassNuclear(A,Z)*Units::u_MeV;Need(mass>0&&std::isfinite(mass)&&Z>0,"Invalid fitted isotope");
    double b=cfg.GetDouble("detector.bfield_t",2.85),e=cfg.GetDouble("detector.efield_vm",55000);
    const auto backend=cfg.GetString("solver.fitter","minuit");Need(backend=="minuit"||backend=="kalman"||backend=="python-objective","solver.fitter must be minuit, python-objective, or kalman");const bool kalman=backend=="kalman",python_objective=backend=="python-objective";
    struct HypothesisSpec{int A;int Z;double mass;std::string label;};
    std::vector<HypothesisSpec> hypothesis_specs;
    if(kalman){
        const auto pairs=ParseKalmanHypotheses(cfg.GetString("solver.kalman_hypotheses",""),A,Z);
        for(const auto& pair:pairs){const double m=MassData::GetNuclearMassNuclear(pair.first,pair.second)*Units::u_MeV;
            Need(m>0&&std::isfinite(m),"No nuclear mass for Kalman hypothesis "+std::to_string(pair.first)+":"+std::to_string(pair.second));
            hypothesis_specs.push_back({pair.first,pair.second,m,IsotopeLabel(pair.first,pair.second)});}
    }
    TrackMesh mesh;MeshConfig mc;
    mc.ke_min=cfg.GetDouble("solver.ke_min",.1);mc.ke_max=cfg.GetDouble("solver.ke_max",70);mc.ke_bins=cfg.GetInt("solver.ke_bins",256);mc.log_ke=cfg.GetBool("solver.log_ke",true);
    mc.polar_min=cfg.GetDouble("solver.polar_min",2);mc.polar_max=cfg.GetDouble("solver.polar_max",88);mc.polar_bins=cfg.GetInt("solver.polar_bins",129);mc.n_time=cfg.GetInt("solver.n_time",400);mc.dt=cfg.GetDouble("solver.dt",3e-10);mc.z_ref=cfg.GetDouble("solver.z_ref",.5);
    const std::string mesh_diagnostics_path=cfg.GetString("solver.mesh_diagnostics_csv","");
    double mesh_load_seconds=0.0,mesh_build_seconds=0.0,mesh_save_seconds=0.0;
    bool mesh_cache_hit=false;
    if(!kalman){
        for(auto key:{"solver.fit_vertex_rho","solver.fit_vertex_phi","solver.fit_azimuthal"})Need(cfg.GetBool(key,true),std::string(key)+"=false is not supported by the native interpolator");
        const auto python_mesh=cfg.GetString("solver.python_mesh","");
        if(python_objective || !python_mesh.empty()){
            const auto mesh_load_start=Clock::now();
            Need(!python_mesh.empty(),"The selected fitting path requires solver.python_mesh from the Python mesh");
            Need(mesh.LoadPythonNpy(python_mesh,mc,A,Z,b,e),"Cannot map the Python NumPy trajectory mesh (missing or mismatched metadata/configuration)");
            mesh_load_seconds=std::chrono::duration<double>(Clock::now()-mesh_load_start).count();
        }else{
            const auto cache=cfg.GetString("solver.mesh_cache","");
            bool loaded=false;
            if(!cache.empty()){
                const auto mesh_load_start=Clock::now();
                loaded=mesh.LoadShared(cache);
                mesh_load_seconds=std::chrono::duration<double>(Clock::now()-mesh_load_start).count();
            }
            const bool ready=loaded&&mesh.MatchesConfig(mc,A,Z,gas,b,e);
            mesh_cache_hit=ready;
            if(!ready){
                const auto mesh_build_start=Clock::now();
                Need(mesh.Build(mc,A,Z,gas,b,e),"Mesh generation failed");
                mesh_build_seconds=std::chrono::duration<double>(Clock::now()-mesh_build_start).count();
                if(!cache.empty()){
                    const auto mesh_save_start=Clock::now();
                    Need(mesh.Save(cache),"Cannot save mesh cache");
                    mesh_save_seconds=std::chrono::duration<double>(Clock::now()-mesh_save_start).count();
                }
            }
        }
    }
    if(!mesh_diagnostics_path.empty()){
        std::ofstream csv(mesh_diagnostics_path);
        Need(static_cast<bool>(csv),"Cannot write mesh diagnostics CSV: "+mesh_diagnostics_path);
        csv<<"backend,cache_hit,cache_load_seconds,build_seconds,save_seconds,total_seconds,n_time,dt,ke_bins,polar_bins\n";
        csv<<backend<<','<<(mesh_cache_hit?1:0)<<','<<std::setprecision(12)<<mesh_load_seconds<<','
           <<mesh_build_seconds<<','<<mesh_save_seconds<<','
           <<(mesh_load_seconds+mesh_build_seconds+mesh_save_seconds)<<','
           <<mc.n_time<<','<<mc.dt<<','<<mc.ke_bins<<','<<mc.polar_bins<<'\n';
        Need(static_cast<bool>(csv),"Failed while writing mesh diagnostics CSV: "+mesh_diagnostics_path);
    }
    KalmanConfig kc;
    kc.meas_sigma_mm=cfg.GetDouble("solver.kalman_meas_sigma",kc.meas_sigma_mm);
    kc.per_cluster_covariance=cfg.GetBool("solver.kalman_per_cluster_covariance",kc.per_cluster_covariance);
    kc.cluster_coef_t_cm2_us=cfg.GetDouble("solver.kalman_cluster_coef_t",kc.cluster_coef_t_cm2_us);
    kc.cluster_coef_l_cm2_us=cfg.GetDouble("solver.kalman_cluster_coef_l",kc.cluster_coef_l_cm2_us);
    kc.cluster_drift_velocity_cm_us=cfg.GetDouble("solver.kalman_cluster_drift_velocity",kc.cluster_drift_velocity_cm_us);
    kc.cluster_tb_time_us=cfg.GetDouble("solver.kalman_cluster_tb_time",kc.cluster_tb_time_us);
    kc.cluster_pad_resolution_xy_mm=cfg.GetDouble("solver.kalman_cluster_pad_resolution_xy",kc.cluster_pad_resolution_xy_mm);
    kc.mom_sigma_frac=cfg.GetDouble("solver.kalman_mom_sigma_frac",kc.mom_sigma_frac);
    kc.pos_sigma_mm=cfg.GetDouble("solver.kalman_pos_sigma",kc.pos_sigma_mm);kc.ang_sigma_deg=cfg.GetDouble("solver.kalman_ang_sigma",kc.ang_sigma_deg);
    kc.min_clusters=cfg.GetInt("solver.kalman_min_clusters",kc.min_clusters);
    kc.min_ke_mev=cfg.GetDouble("solver.kalman_min_ke_mev",kc.min_ke_mev);
    kc.max_steps=cfg.GetInt("solver.kalman_max_steps",kc.max_steps);
    kc.gap_max_steps=cfg.GetInt("solver.kalman_gap_max_steps",kc.gap_max_steps);
    kc.alpha=cfg.GetDouble("solver.kalman_alpha",kc.alpha);kc.beta=cfg.GetDouble("solver.kalman_beta",kc.beta);kc.kappa=cfg.GetDouble("solver.kalman_kappa",kc.kappa);
    kc.smooth=cfg.GetBool("solver.kalman_smooth",kc.smooth);kc.verbose=cfg.GetBool("solver.kalman_verbose",kc.verbose);kc.straggling=cfg.GetBool("solver.kalman_straggling",kc.straggling);kc.mult_scattering=cfg.GetBool("solver.kalman_mult_scattering",kc.mult_scattering);
    kc.back_extrapolate=cfg.GetBool("solver.kalman_back_extrap",kc.back_extrapolate);kc.back_extrap_max_mm=cfg.GetDouble("solver.kalman_back_extrap_max",kc.back_extrap_max_mm);kc.step_mm=cfg.GetDouble("solver.kalman_step_mm",kc.step_mm);kc.chi2_max=cfg.GetDouble("solver.kalman_chi2_max",kc.chi2_max);
    kc.innovation_chi2_max=cfg.GetDouble("solver.kalman_innovation_chi2_max",kc.innovation_chi2_max);
    kc.min_measurement_fraction=cfg.GetDouble("solver.kalman_min_measurement_fraction",kc.min_measurement_fraction);
    kc.iterations=cfg.GetInt("solver.kalman_iterations",kc.iterations);
    kc.ref_covariance_inflation=cfg.GetDouble("solver.kalman_ref_covariance_inflation",kc.ref_covariance_inflation);
    kc.z_pad_plane_mm=cfg.GetDouble("solver.kalman_z_pad_plane_mm",kc.z_pad_plane_mm);
    std::vector<std::vector<Cell>> rows(selected.size());
    const std::string fit_diagnostics_path=cfg.GetString("solver.fit_diagnostics_csv","");
    const double experimental_vertex_radius_limit_mm=
        cfg.GetDouble("solver.experimental_minuit_vertex_radius_mm",0.0);
    const double experimental_vertex_phi_step_mm=
        cfg.GetDouble("solver.experimental_minuit_vertex_phi_step_mm",0.0);
    const bool experimental_centered_vertex_phi_bounds=
        cfg.GetBool("solver.experimental_minuit_centered_vertex_phi_bounds",false);
    const bool centered_azimuth_bounds=cfg.GetBool(
        "solver.minuit_centered_azimuth_bounds",
        cfg.GetBool("solver.experimental_minuit_centered_azimuth_bounds",false));
    const double experimental_azimuth_bounds_center_shift_rad=
        cfg.GetDouble("solver.experimental_minuit_azimuth_bounds_center_shift_rad",0.0);
    const bool experimental_unbounded_azimuth=
        cfg.GetBool("solver.experimental_minuit_unbounded_azimuth",false);
    const bool experimental_fresh_free_vertex_minuit=
        cfg.GetBool("solver.experimental_minuit_fresh_free_vertex",false);
    const bool experimental_exact_minuit_tail_candidates=
        cfg.GetBool("solver.experimental_minuit_exact_tail_candidates",false);
    const bool experimental_mean_squared_distance_objective=cfg.GetBool(
        "solver.experimental_minuit_mean_squared_distance_objective",false);
    std::vector<std::string> fit_diagnostic_rows(fit_diagnostics_path.empty()?0:selected.size());
    struct WorkerHypothesis{
        HypothesisSpec spec;
        std::shared_ptr<KalmanFitter> fitter;
        std::shared_ptr<KalmanFitter> coverage_recovery;
    };
    // Production key takes precedence. Keep the old experimental key as a
    // compatibility override for archived paired controls.
    const bool sort_minuit_candidates = cfg.GetBool(
        "solver.minuit_sort_candidates",
        cfg.GetBool("solver.experimental_minuit_sort_candidates", true));
    Parallel(selected.size(),cfg,[&]{auto fitter=std::make_shared<TrackFitter>();std::vector<WorkerHypothesis> kalman_fitters;
        // Detector config stores positive field magnitudes. Spyral 1.0's
        // propagation convention applies both fields with negative signs.
        const double kalman_bfield=b*cfg.GetDouble("solver.kalman_bfield_sign",-1.0);
        const double kalman_efield=cfg.GetDouble("solver.kalman_efield_vm",e)*
                                   cfg.GetDouble("solver.kalman_efield_sign",-1.0);
        if(kalman){for(const auto& spec:hypothesis_specs){auto kf=std::make_shared<KalmanFitter>();
                Need(kf->Init(gas,spec.A,spec.Z,spec.mass,kalman_bfield,kalman_efield,kc),"Kalman initialization failed for "+spec.label);
                std::shared_ptr<KalmanFitter> recovery;
                if(kc.innovation_chi2_max>0.0&&kc.chi2_max>0.0){
                    KalmanConfig recovery_cfg=kc;
                    recovery_cfg.innovation_chi2_max=0.0;
                    recovery=std::make_shared<KalmanFitter>();
                    Need(recovery->Init(gas,spec.A,spec.Z,spec.mass,kalman_bfield,kalman_efield,recovery_cfg),
                         "Kalman coverage-recovery initialization failed for "+spec.label);
                }
                kalman_fitters.push_back({spec,std::move(kf),std::move(recovery)});}}
        else {Need(fitter->Init(&mesh,mass,Z),"Fitter initialization failed");
            fitter->SetSortCandidates(sort_minuit_candidates);
            fitter->SetExperimentalExactTailCandidates(experimental_exact_minuit_tail_candidates);
            fitter->SetExperimentalVertexRadiusLimit(experimental_vertex_radius_limit_mm);
            fitter->SetExperimentalVertexPhiStepScale(experimental_vertex_phi_step_mm);
            fitter->SetExperimentalCenteredVertexPhiBounds(experimental_centered_vertex_phi_bounds);
            fitter->SetCenteredAzimuthBounds(centered_azimuth_bounds);
            fitter->SetExperimentalAzimuthBoundsCenterShift(experimental_azimuth_bounds_center_shift_rad);
            fitter->SetExperimentalUnboundedAzimuth(experimental_unbounded_azimuth);
            fitter->SetExperimentalFreshFreeVertexMinuit(experimental_fresh_free_vertex_minuit);
            fitter->SetExperimentalMeanSquaredDistanceObjective(experimental_mean_squared_distance_objective);}
        return [&,fitter,kalman_fitters=std::move(kalman_fitters)](size_t j) mutable {ProgressGuard progress_guard{fit_progress};size_t i=selected[j];long long event=estimates.Integer(i,"event"),index=estimates.Integer(i,"cluster_index");const auto& c=events.at(event)->clusters[index];
            PidCandidateAuditRow* audit=collect_pid_audit?&pid_audit[estimate_audit_index[i]]:nullptr;
            if(audit){audit->fit_attempted=1;audit->fit_rejection_stage="fit";}
            EstimateResult guess;guess.ok=true;guess.brho=estimates.Number(i,"brho");guess.polar=estimates.Number(i,"polar");guess.azimuthal=estimates.Number(i,"azimuthal");guess.vertex_x=estimates.Number(i,"vertex_x");guess.vertex_y=estimates.Number(i,"vertex_y");guess.vertex_z=estimates.Number(i,"vertex_z");guess.center_x=estimates.Number(i,"center_x");guess.center_y=estimates.Number(i,"center_y");guess.direction=estimates.Integer(i,"direction");
            FitResult fit;double objective;const double nan=std::numeric_limits<double>::quiet_NaN();
            std::string row_backend=backend,row_objective;
            int fitZ=Z;double fitMass=mass;std::string particle_hypothesis=IsotopeLabel(A,Z);
            if(kalman){
                bool found=false;KalmanFitResult best;double best_score=std::numeric_limits<double>::infinity();
                const HypothesisSpec* best_spec=nullptr;
                for(auto& hypothesis:kalman_fitters){KalmanFitResult candidate;
                    bool ok=hypothesis.fitter->Fit(c.points,guess,candidate)&&candidate.ok&&candidate.ndf>0&&std::isfinite(candidate.chi2);
                    // ATTPCROOT updates every measured cluster instead of
                    // dropping an innovation-gated hit. Use that behavior as
                    // a bounded recovery only when the conservative pass fails
                    // its measurement-coverage requirement. Keep the final
                    // smoothed chi2 gate active, so noisy tracks still fail.
                    if(!ok&&candidate.failure_code==9&&hypothesis.coverage_recovery){
                        KalmanFitResult recovered;
                        const bool recovered_ok=hypothesis.coverage_recovery->Fit(c.points,guess,recovered)&&
                            recovered.ok&&recovered.ndf>0&&std::isfinite(recovered.chi2);
                        if(recovered_ok){candidate=recovered;ok=true;
                            if(kc.verbose)std::fprintf(stderr,
                                "[kalman] event=%lld cluster=%lld hypothesis=%s recovered after coverage failure (chi2/ndf=%.3f)\n",
                                event,index,hypothesis.spec.label.c_str(),candidate.chi2/candidate.ndf);
                        }
                    }
                    if(ok){const double score=candidate.chi2/candidate.ndf;
                        if(!found||score<best_score){found=true;best_score=score;best=candidate;best_spec=&hypothesis.spec;}}
                    else if(kc.verbose)std::fprintf(stderr,"[kalman] event=%lld cluster=%lld hypothesis=%s rejected failure_code=%d n_hits=%d\n",
                        event,index,hypothesis.spec.label.c_str(),candidate.failure_code,candidate.n_hits);
                }
                if(!found||best_spec==nullptr){
                    if(kc.verbose)std::fprintf(stderr,"[kalman] event=%lld cluster=%lld rejected: all %zu particle hypotheses failed\n",event,index,kalman_fitters.size());
                    if(audit)audit->fit_rejection_stage="kalman_hypotheses_rejected";
                    return;
                }
                fitZ=best_spec->Z;fitMass=best_spec->mass;particle_hypothesis=best_spec->label;
                fit.brho=best.brho;fit.sigma_brho=best.sigma_brho;fit.polar=best.polar;fit.sigma_polar=nan;
                fit.azimuthal=best.azimuthal;fit.vx=best.vx;fit.vy=best.vy;fit.vz=best.vz;
                objective=best_score;row_objective="chi2_per_ndf";
            }
            else if(python_objective){
                // Config stores detector length in metres; keep z_error in metres like Spyral.
                const double z_error=(cfg.GetDouble("detector.length_m",1.0)/
                    (cfg.GetDouble("detector.windows_edge",336)-cfg.GetDouble("detector.micromegas_edge",26)))*.5;
                const double beam_radius=cfg.GetDouble("detector.beam_region_radius",25)*.001;
                if(!fitter->FitPythonEquivalent(c.points,guess,z_error,beam_radius,fit)){
                    if(audit){audit->fit_rejection_stage="python_objective_rejected";audit->fit_cost_mm=fit.cost;}return;
                }
                objective=fit.cost;
            }
            else {
                if (fit_diagnostics_path.empty()) {
                    const bool fit_return=fitter->Fit(c.points,guess,guess.azimuthal,fit);
                    if(audit){audit->fit_cost_mm=fit.cost;audit->migrad_status=fit.status;audit->migrad_edm=fit.edm;}
                    if(!fit_return){if(audit)audit->fit_rejection_stage="minuit_acceptance_gate";return;}
                } else {
                    MinuitFitDiagnostics diagnostics;
                    const bool fit_return=fitter->Fit(c.points,guess,guess.azimuthal,fit,&diagnostics);
                    const auto csv_double=[](std::ostringstream& out,double value){
                        if(std::isfinite(value))out<<std::setprecision(17)<<value;
                        else out<<"nan";
                    };
                    std::ostringstream lines;
                    const auto write_row=[&](const MinuitStartDiagnostics* start){
                        lines<<event<<','<<index<<','<<c.label<<','<<diagnostics.point_count<<','
                             <<fit_return<<','<<diagnostics.accepted<<','<<diagnostics.rejection_reason<<',';
                        csv_double(lines,diagnostics.final_cost);
                        lines<<',';csv_double(lines,diagnostics.final_optimizer_objective);
                        lines<<','<<diagnostics.best_start_index;
                        if (!start) {
                            for(int field=0;field<48;++field)lines<<',';
                            lines<<',';csv_double(lines,experimental_vertex_radius_limit_mm);
                            for(int field=0;field<30;++field)lines<<',';
                            for(int field=0;field<12;++field)lines<<',';
                            for(int field=0;field<19;++field)lines<<',';
                            lines<<'\n';
                            return;
                        }
                        lines<<','<<start->start_index<<',';csv_double(lines,start->azimuth_offset);lines<<',';
                        csv_double(lines,start->initial_azimuth);lines<<','
                             <<start->fixed_vertex.minimize_return<<','<<start->fixed_vertex.status<<','
                             <<start->fixed_vertex.covariance_status<<','<<start->fixed_vertex.function_calls<<',';
                        csv_double(lines,start->fixed_vertex.objective);lines<<',';csv_double(lines,start->fixed_vertex.edm);
                        lines<<','<<start->free_vertex.minimize_return<<','<<start->free_vertex.status<<','
                             <<start->free_vertex.covariance_status<<','<<start->free_vertex.function_calls<<',';
                        csv_double(lines,start->free_vertex.objective);lines<<',';csv_double(lines,start->free_vertex.edm);
                        lines<<','<<start->below_existing_cost_gate<<','<<start->selected_as_best<<',';
                        csv_double(lines,start->mean_distance_cost);
                        for(double value:start->initial_values){lines<<',';csv_double(lines,value);}
                        for(double value:start->lower_bounds){lines<<',';csv_double(lines,value);}
                        for(double value:start->upper_bounds){lines<<',';csv_double(lines,value);}
                        for(double value:start->final_values){lines<<',';csv_double(lines,value);}
                        for(double value:start->final_errors){lines<<',';csv_double(lines,value);}
                        lines<<',';csv_double(lines,experimental_vertex_radius_limit_mm);
                        const auto write_search_counters=[&](const CandidateSearchCounters& counters){
                            lines<<','<<counters.objective_evaluations_with_curve
                                 <<','<<counters.objective_evaluations_with_flat_tail
                                 <<','<<counters.flat_tail_segments_total
                                 <<','<<counters.flat_tail_segments_max
                                 <<','<<counters.evaluated_hits
                                 <<','<<counters.hits_with_tail_candidate
                                 <<','<<counters.hits_with_multiple_tail_candidates
                                 <<','<<counters.tail_candidate_slots
                                 <<','<<counters.multiple_tail_hits_checked_exactly
                                 <<','<<counters.multiple_tail_hits_missed_exact
                                 <<','<<std::setprecision(17)<<counters.multiple_tail_miss_distance_sum_mm
                                 <<','<<std::setprecision(17)<<counters.multiple_tail_miss_distance_max_mm
                                 <<','<<counters.objective_evaluations_with_multiple_tail_miss
                                 <<','<<std::setprecision(17)<<counters.multiple_tail_miss_objective_delta_sum_mm
                                 <<','<<std::setprecision(17)<<counters.multiple_tail_miss_objective_delta_max_mm;
                        };
                        write_search_counters(start->fixed_vertex.candidate_search);
                        write_search_counters(start->free_vertex.candidate_search);
                        for(double value:start->fixed_vertex.final_values){lines<<',';csv_double(lines,value);}
                        for(double value:start->fixed_vertex.final_errors){lines<<',';csv_double(lines,value);}
                        lines<<','<<start->fresh_free_vertex_ran<<','<<start->fresh_free_vertex.minimize_return
                             <<','<<start->fresh_free_vertex.status<<','<<start->fresh_free_vertex.covariance_status
                             <<','<<start->fresh_free_vertex.function_calls<<',';
                        csv_double(lines,start->fresh_free_vertex.objective);lines<<',';
                        csv_double(lines,start->fresh_free_vertex.edm);
                        for(double value:start->fresh_free_vertex.final_values){lines<<',';csv_double(lines,value);}
                        for(double value:start->fresh_free_vertex.final_errors){lines<<',';csv_double(lines,value);}
                        lines<<'\n';
                    };
                    if (diagnostics.starts.empty()) write_row(nullptr);
                    else for (const auto& start:diagnostics.starts) write_row(&start);
                    fit_diagnostic_rows[j]=lines.str();
                    if(audit){audit->fit_cost_mm=fit.cost;audit->migrad_status=fit.status;audit->migrad_edm=fit.edm;}
                    if(!fit_return){if(audit)audit->fit_rejection_stage="minuit_acceptance_gate";return;}
                }
                objective=fit.cost;
                if(experimental_mean_squared_distance_objective)
                    row_objective="mean_distance_mm_gate;minuit_mean_squared_distance_mm2";
            }
            if(audit){audit->fit_state=FitState::FIT_ACCEPTED;
                audit->fit_rejection_stage="accepted";audit->fit_cost_mm=fit.cost;
                if(kalman){audit->migrad_status=-1;audit->migrad_edm=-1.0;}}
            double p=299.792458*fitZ*fit.brho,energy=std::hypot(p,fitMass),ke=p*p/(energy+fitMass),sigma_ke=std::abs(p/energy*299.792458*fitZ*fit.sigma_brho);
            std::vector<Cell> row;for(auto key:{"event","cluster_index","cluster_label","orig_run","orig_event"})row.emplace_back(estimates.Integer(i,key));
            // Position m matches Python SolverResult; native unavailable covariance is NaN.
            for(double v:{fit.vx*.001,nan,fit.vy*.001,nan,fit.vz*.001,nan,fit.brho,fit.sigma_brho,ke,sigma_ke,fit.polar,fit.sigma_polar,fit.azimuthal,nan,objective})row.emplace_back(v);
            row.emplace_back(python_objective?"minuit2_python_objective":row_backend);
            row.emplace_back(!row_objective.empty()?row_objective:python_objective?"python_weighted_nearest_distance_m":"mean_distance_mm");
            row.emplace_back(particle_hypothesis);rows[j]=std::move(row);
        };});
    if(!fit_diagnostics_path.empty()){
        std::ofstream csv(fit_diagnostics_path);
        Need(static_cast<bool>(csv),"Cannot write Minuit diagnostics CSV: "+fit_diagnostics_path);
        csv<<"event,cluster_index,cluster_label,n_points,fit_return,fit_accepted,rejection_reason,best_cost_mm,best_optimizer_objective_value,best_start_index,start_index,azimuth_offset_rad,initial_azimuth_rad,fixed_minimize_return,fixed_status,fixed_covariance_status,fixed_ncalls,fixed_objective_mm,fixed_edm,free_minimize_return,free_status,free_covariance_status,free_ncalls,free_objective_mm,free_edm,below_20mm_gate,selected_as_best,mean_distance_cost_mm"
              ",initial_brho_Tm,initial_polar_rad,initial_vx_mm,initial_vy_mm,initial_vz_mm,initial_azimuth_rad"
              ",lower_brho_Tm,lower_polar_rad,lower_vx_mm,lower_vy_mm,lower_vz_mm,lower_azimuth_rad"
              ",upper_brho_Tm,upper_polar_rad,upper_vx_mm,upper_vy_mm,upper_vz_mm,upper_azimuth_rad"
              ",final_brho_Tm,final_polar_rad,final_vx_mm,final_vy_mm,final_vz_mm,final_azimuth_rad"
              ",error_brho,error_polar,error_vx,error_vy,error_vz,error_azimuth,experimental_vertex_radius_limit_mm"
              ",fixed_eval_curve_count,fixed_eval_flat_tail_count,fixed_flat_tail_segments_sum,fixed_flat_tail_segments_max,fixed_eval_hit_count,fixed_hits_with_tail_candidate,fixed_hits_with_multiple_tail_candidates,fixed_tail_candidate_slots"
              ",fixed_multiple_tail_exact_checks,fixed_multiple_tail_missed_exact,fixed_multiple_tail_miss_sum_mm,fixed_multiple_tail_miss_max_mm,fixed_objective_evaluations_with_multiple_tail_miss,fixed_multiple_tail_objective_delta_sum_mm,fixed_multiple_tail_objective_delta_max_mm"
              ",free_eval_curve_count,free_eval_flat_tail_count,free_flat_tail_segments_sum,free_flat_tail_segments_max,free_eval_hit_count,free_hits_with_tail_candidate,free_hits_with_multiple_tail_candidates,free_tail_candidate_slots"
              ",free_multiple_tail_exact_checks,free_multiple_tail_missed_exact,free_multiple_tail_miss_sum_mm,free_multiple_tail_miss_max_mm,free_objective_evaluations_with_multiple_tail_miss,free_multiple_tail_objective_delta_sum_mm,free_multiple_tail_objective_delta_max_mm"
              ",fixed_final_brho_Tm,fixed_final_polar_rad,fixed_final_vx_mm,fixed_final_vy_mm,fixed_final_vz_mm,fixed_final_azimuth_rad"
              ",fixed_error_brho,fixed_error_polar,fixed_error_vx,fixed_error_vy,fixed_error_vz,fixed_error_azimuth"
              ",fresh_free_ran,fresh_free_minimize_return,fresh_free_status,fresh_free_covariance_status,fresh_free_ncalls,fresh_free_objective_mm,fresh_free_edm"
              ",fresh_free_brho_Tm,fresh_free_polar_rad,fresh_free_vx_mm,fresh_free_vy_mm,fresh_free_vz_mm,fresh_free_azimuth_rad"
              ",fresh_free_error_brho,fresh_free_error_polar,fresh_free_error_vx,fresh_free_error_vy,fresh_free_error_vz,fresh_free_error_azimuth\n";
        for(const auto& lines:fit_diagnostic_rows)csv<<lines;
        Need(static_cast<bool>(csv),"Failed while writing Minuit diagnostics CSV: "+fit_diagnostics_path);
    }
    write_pid_audit();
    for(auto& row:rows)if(!row.empty())table.rows.push_back(std::move(row));WriteScalarStage(output,"InterpSolver",table);std::cout<<"InterpSolver: "<<table.rows.size()<<" accepted fits ("<<backend<<")\n";
}
}
int main(int argc,char** argv){try{
    ROOT::EnableThreadSafety();Need(argc>=2,"Usage: rk4_workspace config|schema|convert|cluster|estimate|solve ...");std::string cmd=argv[1];auto start=Clock::now();
    if(cmd=="config"){Need(argc==3,"config CONFIG");auto cfg=Load(argv[2]);std::cout<<'{';bool first=true;for(auto& [k,v]:cfg.Values()){if(!first)std::cout<<',';first=false;std::cout<<Json(k)<<':'<<Json(v);}std::cout<<"}\n";return 0;}
    if(cmd=="schema"){Need(argc==3,"schema STAGE");const auto& s=std::string(argv[2])=="Estimation"?EstimationColumns():SolverColumns();std::cout<<'[';for(size_t i=0;i<s.size();++i){if(i)std::cout<<',';std::cout<<'['<<Json(s[i].name)<<','<<Json(s[i].type==ColumnType::Integer?"int64":s[i].type==ColumnType::Real?"float64":"string")<<']';}std::cout<<"]\n";return 0;}
    if(cmd=="convert"){Need(argc==5,"convert STAGE INPUT OUTPUT");std::string stage=argv[2];if(stage=="Pointcloud")WritePointcloudStage(argv[4],ReadPointcloudStage(argv[3]));else if(stage=="Cluster")WriteClusterStage(argv[4],ReadClusterStage(argv[3]));else{Need(stage=="Estimation"||stage=="InterpSolver","Unknown stage");WriteScalarStage(argv[4],stage,ReadScalarStage(argv[3],stage));}}
    else if(cmd=="cluster"){Need(argc==5,"cluster INPUT OUTPUT CONFIG");ClusterStage(argv[2],argv[3],Load(argv[4]));}
    else if(cmd=="estimate"){Need(argc==5,"estimate INPUT OUTPUT CONFIG");EstimateStage(argv[2],argv[3],Load(argv[4]));}
    else if(cmd=="solve"){Need(argc==6,"solve CLUSTERS ESTIMATES OUTPUT CONFIG");SolveStage(argv[2],argv[3],argv[4],Load(argv[5]));}
    else throw std::runtime_error("Unknown command "+cmd);
    std::cout<<std::setprecision(9)<<"Native stage elapsed: "<<std::chrono::duration<double>(Clock::now()-start).count()<<" s\n";return 0;
}catch(const std::exception& e){std::cerr<<"SpyralPP: "<<e.what()<<'\n';return 1;}}
