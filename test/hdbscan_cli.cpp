#include "reco/Hdbscan.h"
#include <iostream>
int main(){int n,k,size,single;double epsilon,zscale;
while(std::cin>>n>>k>>size>>single>>epsilon>>zscale){
    rk4::PointCloud cloud;cloud.points.resize(n);
    for(auto& p:cloud.points)if(!(std::cin>>p.x>>p.y>>p.z))return 2;
    rk4::HdbscanClusterer hdb;hdb.Run(cloud,k,size,single!=0,epsilon,zscale);
    for(int label:hdb.point_labels_)std::cout<<label<<' ';std::cout<<'\n';
}return 0;}
