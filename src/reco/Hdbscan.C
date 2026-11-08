// Adapted from scikit-learn 1.7.2 cluster/_hdbscan/_linkage.pyx and _tree.pyx.
// Authors: The scikit-learn developers. SPDX-License-Identifier: BSD-3-Clause
// Copyright conditions and disclaimer: notices/scikit-learn-hdbscan.txt.
#include "Hdbscan.h"
#include "NumpyArgsort.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace rk4 {
namespace {
constexpr double Infinity=std::numeric_limits<double>::infinity();
struct Edge {int a,b;double distance;};
struct Node {int left=-1,right=-1,size=1;double distance=0;};
struct Condensed {int parent,child;double lambda;int size;};
double Distance(const CloudPoint&a,const CloudPoint&b,double z_scale){
    const double x=a.x-b.x,y=a.y-b.y,z=a.z*z_scale-b.z*z_scale;
    return std::sqrt(x*x+y*y+z*z);
}
int Find(std::vector<int>& parent,int node){
    int root=node;while(parent[root]!=root)root=parent[root];
    while(parent[node]!=node){const int next=parent[node];parent[node]=root;node=next;}
    return root;
}
std::vector<int> Breadth(const std::vector<Node>& nodes,int root){
    std::vector<int> order{root};
    for(size_t i=0;i<order.size();++i){
        const auto& node=nodes[order[i]];
        if(node.left>=0){order.push_back(node.left);order.push_back(node.right);}
    }
    return order;
}
std::vector<Node> Hierarchy(const PointCloud& cloud,int samples,double z_scale){
    const int n=static_cast<int>(cloud.points.size());
    std::vector<double> core(n),row(n);
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j)row[j]=Distance(cloud.points[i],cloud.points[j],z_scale);
        std::nth_element(row.begin(),row.begin()+samples-1,row.end());
        core[i]=row[samples-1];
    }
    // The source tracking and strict tie comparisons follow sklearn's Prim
    // implementation. Keep distances in distance units throughout.
    std::vector<char> visited(n,0);
    std::vector<double> reach(n,Infinity);
    std::vector<int> sources(n,1);
    std::vector<Edge> edges;edges.reserve(n-1);
    int current=0;
    for(int iteration=0;iteration<n-1;++iteration){
        visited[current]=1;
        double next_distance=Infinity;int source=0,next=0;
        for(int j=0;j<n;++j){
            if(visited[j])continue;
            const double previous=reach[j];
            const double distance=std::max({core[current],core[j],Distance(cloud.points[current],cloud.points[j],z_scale)});
            if(distance<previous){
                reach[j]=distance;sources[j]=current;
                if(distance<next_distance){next_distance=distance;source=current;next=j;}
            }else if(previous<next_distance){next_distance=previous;source=sources[j];next=j;}
        }
        edges.push_back({source,next,next_distance});current=next;
    }
    std::vector<double> weights;weights.reserve(edges.size());
    for(const auto& edge:edges)weights.push_back(edge.distance);
    const auto order=NumpyArgsort(weights);
    const auto unordered=edges;
    for(size_t i=0;i<order.size();++i)edges[i]=unordered[order[i]];
    std::vector<Node> nodes(2*n-1);
    std::vector<int> parent(2*n-1);std::iota(parent.begin(),parent.end(),0);
    for(int i=0;i<n-1;++i){
        const int a=Find(parent,edges[i].a),b=Find(parent,edges[i].b),id=n+i;
        nodes[id]={a,b,nodes[a].size+nodes[b].size,edges[i].distance};
        parent[a]=id;parent[b]=id;
    }
    return nodes;
}
std::vector<Condensed> Condense(const std::vector<Node>& nodes,int n,int min_size,int& next_label){
    const int root=2*n-2;
    std::vector<int> relabel(nodes.size(),-1);relabel[root]=n;next_label=n+1;
    std::vector<char> ignore(nodes.size(),0);
    std::vector<Condensed> result;
    auto prune=[&](int node,int label,double lambda){
        for(int child:Breadth(nodes,node)){
            if(child<n)result.push_back({label,child,lambda,1});
            ignore[child]=1;
        }
    };
    for(int index:Breadth(nodes,root)){
        if(ignore[index]||index<n)continue;
        const auto& node=nodes[index];const int left=node.left,right=node.right;
        const double lambda=node.distance>0?1/node.distance:Infinity;
        const bool keep_left=nodes[left].size>=min_size,keep_right=nodes[right].size>=min_size;
        if(keep_left&&keep_right){
            relabel[left]=next_label++;relabel[right]=next_label++;
            result.push_back({relabel[index],relabel[left],lambda,nodes[left].size});
            result.push_back({relabel[index],relabel[right],lambda,nodes[right].size});
        }else if(!keep_left&&!keep_right){
            prune(left,relabel[index],lambda);prune(right,relabel[index],lambda);
        }else if(!keep_left){relabel[right]=relabel[index];prune(left,relabel[index],lambda);}
        else{relabel[left]=relabel[index];prune(right,relabel[index],lambda);}
    }
    return result;
}
std::vector<int> Select(const std::vector<Condensed>& tree,int n,int next_label,bool single,double epsilon){
    std::vector<double> births(next_label,0),stability(next_label,0);
    std::vector<int> parent(next_label,-1),point_parent(n,-1);
    std::vector<double> point_lambda(n,0);
    std::vector<std::vector<int>> children(next_label);
    for(const auto& edge:tree){
        if(edge.size>1){births[edge.child]=edge.lambda;parent[edge.child]=edge.parent;children[edge.parent].push_back(edge.child);}
        else{point_parent[edge.child]=edge.parent;point_lambda[edge.child]=edge.lambda;}
    }
    for(const auto& edge:tree)stability[edge.parent]+=(edge.lambda-births[edge.parent])*edge.size;
    std::vector<char> selected(next_label,0);
    for(int node=n+(single?0:1);node<next_label;++node)selected[node]=1;
    for(int node=next_label-1;node>=n+(single?0:1);--node){
        double subtree=0;for(int child:children[node])subtree+=stability[child];
        if(subtree>stability[node]){selected[node]=0;stability[node]=subtree;}
        else{
            std::vector<int> descendants=children[node];
            for(size_t i=0;i<descendants.size();++i){
                const int child=descendants[i];selected[child]=0;
                descendants.insert(descendants.end(),children[child].begin(),children[child].end());
            }
        }
    }
    if(epsilon>0){
        std::vector<char> merged(next_label,0);
        for(int node=n;node<next_label;++node){
            if(!selected[node])continue;
            int target=node;
            if(node!=n && 1/births[node]<epsilon){
                while(true){
                    const int p=parent[target];
                    if(p==n){if(single)target=n;break;}
                    target=p;
                    if(1/births[target]>epsilon)break;
                }
            }
            merged[target]=1;
        }
        for(int node=n;node<next_label;++node){
            if(!merged[node])continue;
            for(int p=parent[node];p>=n;p=parent[p])if(merged[p]){merged[node]=0;break;}
        }
        selected.swap(merged);
    }
    std::vector<int> map(next_label,-1);int count=0;
    for(int node=n;node<next_label;++node)if(selected[node])map[node]=count++;
    double root_threshold=epsilon>0?1/epsilon:0;
    if(epsilon==0)for(const auto& edge:tree)if(edge.parent==n)root_threshold=std::max(root_threshold,edge.lambda);
    std::vector<int> labels(n,-1);
    for(int point=0;point<n;++point){
        int node=point_parent[point];
        while(node>n&&!selected[node])node=parent[node];
        if(node>n&&selected[node])labels[point]=map[node];
        else if(node==n&&selected[n]&&count==1&&single&&point_lambda[point]>=root_threshold)labels[point]=map[n];
    }
    return labels;
}
} // namespace

int HdbscanClusterer::Run(const PointCloud& cloud,int samples,int min_size,bool single,double epsilon,double z_scale){
    clusters_.clear();const int n=static_cast<int>(cloud.points.size());point_labels_.assign(n,-1);
    if(samples<1||min_size<2||epsilon<0||!std::isfinite(epsilon)||!(z_scale>0)||!std::isfinite(z_scale))
        throw std::invalid_argument("Invalid HDBSCAN parameters");
    if(n<min_size||n<samples||n<2)return 0;
    for(const auto& p:cloud.points)if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z))
        throw std::invalid_argument("HDBSCAN coordinates must be finite");
    const auto hierarchy=Hierarchy(cloud,samples,z_scale);int next_label=0;
    const auto tree=Condense(hierarchy,n,min_size,next_label);
    point_labels_=Select(tree,n,next_label,single,epsilon);
    int count=0;for(int label:point_labels_)count=std::max(count,label+1);
    clusters_.resize(count);
    for(int i=0;i<count;++i)clusters_[i].label=i;
    for(int i=0;i<n;++i)if(point_labels_[i]>=0)clusters_[point_labels_[i]].points.push_back(cloud.points[i]);
    return count;
}
} // namespace rk4
