#ifndef SPYRALPP_EVENT_METADATA_H
#define SPYRALPP_EVENT_METADATA_H
namespace rk4 {
struct EventMetadata {
    long long orig_run = -1, orig_event = -1;
    double ic_amplitude = -1, ic_centroid = -1, ic_integral = -1, ic_multiplicity = -1;
};
}
#endif
