#include "io/WorkspaceIO.h"
#include "reco/H5TraceReader.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <unistd.h>
using namespace rk4;
void Check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){
 auto dir=std::filesystem::temp_directory_path()/("spyral_workspace_io_"+std::to_string(getpid()));std::filesystem::create_directory(dir);
 try{
  PointcloudData points;points.info.run=151;points.info.first_event=2;points.info.last_event=5;points.info.raw_events=4;points.info.excluded_ic=1;
  PointCloud cloud;cloud.event=2;cloud.metadata.orig_run=151;cloud.metadata.orig_event=9007199254740993LL;cloud.metadata.ic_amplitude=650;cloud.metadata.ic_centroid=1150;cloud.metadata.ic_integral=9876;cloud.metadata.ic_multiplicity=1;
  CloudPoint p;p.x=1.25;p.y=-2.5;p.z=700;p.amplitude=123;p.integral=456;p.pad=25;p.tb=123.125;p.scale=2;cloud.points.push_back(p);points.events.push_back(cloud);cloud.event=5;cloud.points.clear();cloud.metadata.ic_amplitude=-1;cloud.metadata.ic_centroid=-1;cloud.metadata.ic_integral=-1;cloud.metadata.ic_multiplicity=-1;points.events.push_back(cloud);
  ClusterData clusters;clusters.info=points.info;ClusterEvent e;e.event=2;e.metadata=points.events[0].metadata;Cluster c;c.label=0;c.original_label=9;c.direction=Direction::Backward;c.points.push_back(p);e.clusters.push_back(c);clusters.events.push_back(e);e.event=5;e.clusters.clear();e.metadata=cloud.metadata;clusters.events.push_back(e);
  std::vector<std::string> formats={".root"};if(H5TraceReader::Available())formats.push_back(".h5");
  for(auto extension:formats){
   auto a=(dir/("points"+extension)).string(),b=(dir/("clusters"+extension)).string();WritePointcloudStage(a,points);WriteClusterStage(b,clusters);
   auto pa=ReadPointcloudStage(a);auto cb=ReadClusterStage(b);
   Check(pa.events.size()==2&&pa.events[1].points.empty(),"Empty cloud lost");Check(pa.events[0].metadata.orig_event==9007199254740993LL,"Metadata ID lost precision");
   Check(pa.events[0].metadata.ic_amplitude==650&&pa.events[1].metadata.ic_amplitude==-1,"IC metadata changed");Check(pa.events[0].points[0].tb==123.125&&pa.events[0].points[0].scale==2,"Point values changed");
   Check(cb.events.size()==2&&cb.events[1].clusters.empty(),"Empty cluster event lost");Check(cb.events[0].clusters[0].original_label==9,"Original HDBSCAN label lost");Check(cb.events[0].clusters[0].direction==Direction::Backward,"Cluster direction lost");
   Check(pa.info.raw_events==4&&cb.info.excluded_ic==1&&cb.info.first_event==2&&cb.info.last_event==5,"Stage provenance changed");
  }
  for(auto stage:{"Estimation","InterpSolver"}){
   ScalarTable table;table.info=points.info;table.columns=std::string(stage)=="Estimation"?EstimationColumns():SolverColumns();std::vector<Cell> row;
   for(const auto& col:table.columns){if(col.type==ColumnType::Integer)row.emplace_back(9007199254740993LL);else if(col.type==ColumnType::Real)row.emplace_back(col.name=="sigma_vx"?NAN:1.2345678901234567);else row.emplace_back(std::string("diagnostic \"quoted\", text"));}
   table.rows.push_back(row);
   for(auto ext:{".root",".csv"}){
    auto path=(dir/(std::string(stage)+ext)).string();WriteScalarStage(path,stage,table);auto back=ReadScalarStage(path,stage);
    Check(back.rows.size()==1&&back.Integer(0,"event")==9007199254740993LL,"Scalar Int64 lost precision");
    for(size_t i=0;i<row.size();++i){if(std::holds_alternative<double>(row[i])&&std::isnan(std::get<double>(row[i])))Check(std::isnan(std::get<double>(back.rows[0][i])),"NaN changed");else Check(row[i]==back.rows[0][i],"Scalar roundtrip changed value");}
   }
   table.rows.clear();auto path=(dir/(std::string(stage)+"_empty.root")).string();WriteScalarStage(path,stage,table);auto back=ReadScalarStage(path,stage);Check(back.rows.empty()&&back.columns.size()==table.columns.size(),"Empty scalar schema lost");
  }
  std::filesystem::remove_all(dir);std::cout<<"workspace I/O: both formats, empty events/tables, metadata, labels, Int64, NaN and text pass\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';std::filesystem::remove_all(dir);return 1;}
}
