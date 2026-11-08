#ifndef SPYRALPP_EVENT_META_IO_H
#define SPYRALPP_EVENT_META_IO_H
#include "../core/EventMetadata.h"
class TTree;
namespace rk4 {
void BindEventMetadata(TTree& tree,EventMetadata& metadata);
EventMetadata ReadEventMetadata(TTree& tree,long long event);
}
#endif
