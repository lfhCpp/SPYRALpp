// Exercise the actual independent solver executable using a known mesh trajectory.
#include "io/WorkspaceIO.h"
#include "reco/TrackMesh.h"
#include "core/MassData.h"
#include "core/Units.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
using namespace rk4;
void Check(bool good,const char* message){if(!good)throw std::runtime_error(message);}
int main(int argc,char** argv){
 auto dir=std::filesystem::temp_directory_path()/("spyral_workspace_solver_"+std::to_string(getpid()));std::filesystem::create_directory(dir);
 try{
  Check(argc==2,"Need stage engine path");Element h;h.a=1;h.z=1;h.s=2;GasTarget gas({h},600);TrackMesh mesh;MeshConfig mc;mc.ke_bins=12;mc.ke_min=2;mc.ke_max=20;mc.polar_bins=12;mc.polar_min=20;mc.polar_max=80;mc.n_time=300;
  Check(mesh.Build(mc,1,1,gas,2.85,55000),"Cannot build test mesh");std::vector<double> xyz;Check(mesh.Interpolate(8,45*Units::DEG2RAD,xyz),"Cannot interpolate known trajectory");
  ClusterData data;data.info.run=151;data.info.first_event=7;data.info.last_event=7;ClusterEvent e;e.event=7;e.metadata.orig_run=151;e.metadata.orig_event=42;e.metadata.ic_amplitude=650;Cluster c;c.label=0;c.original_label=9;c.direction=Direction::Forward;
  for(int i=0;i<120;i+=2){CloudPoint p;p.x=xyz[3*i];p.y=xyz[3*i+1];p.z=xyz[3*i+2];p.integral=100;p.scale=1;c.points.push_back(p);}e.clusters.push_back(c);data.events.push_back(e);
  const auto clusters=(dir/"clusters.root").string(),estimates=(dir/"estimates.root").string(),output=(dir/"fit.root").string(),config=(dir/"config.cfg").string();WriteClusterStage(clusters,data);
  const double mass=MassData::GetNuclearMassNuclear(1,1)*Units::u_MeV,brho=std::sqrt(8*(8+2*mass))/299.792458;
  ScalarTable table;table.info=data.info;table.columns=EstimationColumns();std::vector<Cell> row;
  for(auto& col:table.columns){if(col.type==ColumnType::Integer){long long v=0;if(col.name=="event")v=7;if(col.name=="cluster_label")v=9;if(col.name=="orig_event")v=42;if(col.name=="orig_run")v=151;row.emplace_back(v);}else{double v=1;if(col.name=="vertex_x"||col.name=="vertex_y"||col.name=="azimuthal")v=0;if(col.name=="vertex_z"||col.name=="center_z")v=500;if(col.name=="brho")v=brho;if(col.name=="polar")v=45*Units::DEG2RAD;if(col.name=="ic_amplitude")v=650;row.emplace_back(v);}}table.rows.push_back(row);WriteScalarStage(estimates,"Estimation",table);
  std::ofstream cfg(config);cfg<<"run.threads=1\nsolver.fitter=minuit\nsolver.ke_bins=12\nsolver.ke_min=2\nsolver.ke_max=20\nsolver.polar_bins=12\nsolver.polar_min=20\nsolver.polar_max=80\nsolver.n_time=300\nsolver.ic_min_val=500\nsolver.ic_max_val=800\n";cfg.close();
  auto pid=fork();Check(pid>=0,"Cannot fork stage engine");if(pid==0){execl(argv[1],argv[1],"solve",clusters.c_str(),estimates.c_str(),output.c_str(),config.c_str(),nullptr);_exit(127);}int status=0;Check(waitpid(pid,&status,0)==pid&&WIFEXITED(status)&&WEXITSTATUS(status)==0,"Stage engine failed");
  auto fit=ReadScalarStage(output,"InterpSolver");Check(fit.rows.size()==1,"Known good independent solver input was rejected");Check(fit.Integer(0,"orig_event")==42&&fit.Integer(0,"cluster_label")==9,"Original IDs changed");Check(std::abs(fit.Number(0,"brho")/brho-1)<.005,"Independent solver brho error exceeds .5 percent");Check(std::abs(fit.Number(0,"vertex_z")-.5)<.003,"Solver position is not in meters");Check(std::abs(fit.Number(0,"ke")-8)<.08,"Solver kinetic energy changed");Check(std::get<std::string>(fit.rows[0][fit.Index("objective_kind")])=="mean_distance_mm","Solver objective units missing");
  std::filesystem::remove_all(dir);std::cout<<"Independent Minuit2 solver: known trajectory, IDs, meters, energy and objective pass\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';std::filesystem::remove_all(dir);return 1;}
}
