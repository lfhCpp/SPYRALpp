#include "EventMetaIO.h"
#include "TTree.h"
#include "TLeaf.h"
namespace rk4 {
void BindEventMetadata(TTree& tree,EventMetadata& m){
    tree.Branch("orig_run",&m.orig_run,"orig_run/L");tree.Branch("orig_event",&m.orig_event,"orig_event/L");
    tree.Branch("ic_amplitude",&m.ic_amplitude,"ic_amplitude/D");
    tree.Branch("ic_centroid",&m.ic_centroid,"ic_centroid/D");
    tree.Branch("ic_integral",&m.ic_integral,"ic_integral/D");
    tree.Branch("ic_multiplicity",&m.ic_multiplicity,"ic_multiplicity/D");
}
EventMetadata ReadEventMetadata(TTree& tree,long long event){
    EventMetadata m;m.orig_event=event;
    if(auto* leaf=tree.GetLeaf("orig_run"))m.orig_run=leaf->GetValueLong64();
    if(auto* leaf=tree.GetLeaf("orig_event"))m.orig_event=leaf->GetValueLong64();
    if(auto* leaf=tree.GetLeaf("ic_amplitude"))m.ic_amplitude=leaf->GetValue();
    if(auto* leaf=tree.GetLeaf("ic_centroid"))m.ic_centroid=leaf->GetValue();
    if(auto* leaf=tree.GetLeaf("ic_integral"))m.ic_integral=leaf->GetValue();
    if(auto* leaf=tree.GetLeaf("ic_multiplicity"))m.ic_multiplicity=leaf->GetValue();
    return m;
}
}
