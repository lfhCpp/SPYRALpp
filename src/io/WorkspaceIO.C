#include "WorkspaceIO.h"
#include "CloudIO.h"
#include "EventMetaIO.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TNamed.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>
namespace rk4 {
namespace {
void Require(bool value,const std::string& text){if(!value)throw std::runtime_error(text);}
bool Root(const std::string& path){return path.size()>=5&&path.substr(path.size()-5)==".root";}
bool Csv(const std::string& path){return path.size()>=4&&path.substr(path.size()-4)==".csv";}
long long Parameter(TFile& file,const char* name,long long fallback){auto* p=dynamic_cast<TParameter<Long64_t>*>(file.Get(name));return p?p->GetVal():fallback;}
std::string Named(TFile& file,const char* name,const std::string& fallback){auto* p=dynamic_cast<TNamed*>(file.Get(name));return p?p->GetTitle():fallback;}
StageFileInfo ReadInfo(TFile& file,long long first,long long last){
    StageFileInfo i;i.first_event=Parameter(file,"min_event",first);i.last_event=Parameter(file,"max_event",last);
    i.run=Parameter(file,"run_number",-1);i.raw_events=Parameter(file,"raw_events",0);i.excluded_ic=Parameter(file,"excluded_ic",0);i.good_downscale=Parameter(file,"good_ic_downscale",0);
    i.seed=std::stoull(Named(file,"dither_seed","20261026"));i.waveform_method=Named(file,"waveform_method","baseline");return i;
}
void WriteInfo(TFile& file,const StageFileInfo& i,const std::string& stage){
    file.cd();
    for(const auto& item:std::vector<std::pair<std::string,long long>>{{"min_event",i.first_event},{"max_event",i.last_event},{"run_number",i.run},{"raw_events",i.raw_events},{"excluded_ic",i.excluded_ic},{"good_ic_downscale",i.good_downscale}}){
        TParameter<Long64_t> p(item.first.c_str(),item.second);Require(p.Write(item.first.c_str(),TObject::kOverwrite)>0,"Cannot write ROOT file metadata");
    }
    for(const auto& item:std::vector<std::pair<std::string,std::string>>{{"dither_seed",std::to_string(i.seed)},{"waveform_method",i.waveform_method},{"spyralpp_stage",stage},{"spyralpp_schema","workspace-1"},{"position_unit",stage=="InterpSolver"?"m":"mm"},{"angle_unit","rad"}}){
        TNamed p(item.first.c_str(),item.second.c_str());Require(p.Write(item.first.c_str(),TObject::kOverwrite)>0,"Cannot write ROOT text metadata");
    }
}
TLeaf* Leaf(TTree& tree,const std::string& name,const std::string& alias=""){
    auto* leaf=tree.GetLeaf(name.c_str());return leaf||alias.empty()?leaf:tree.GetLeaf(alias.c_str());
}
long long Integer(TTree& tree,const std::string& name,long long fallback,const std::string& alias=""){
    auto* leaf=Leaf(tree,name,alias);return leaf?leaf->GetValueLong64():fallback;
}
double Real(TTree& tree,const std::string& name,double fallback,const std::string& alias=""){
    auto* leaf=Leaf(tree,name,alias);return leaf?leaf->GetValue():fallback;
}
void BindVector(TTree& tree,const char* name,std::vector<double>** address,bool required){
    if(!tree.GetBranch(name)){Require(!required,"Missing ROOT vector "+std::string(name));return;}
    Require(tree.SetBranchAddress(name,address)>=0,"Incompatible ROOT vector "+std::string(name));
}
std::vector<std::string> ParseCsv(const std::string& line){
    std::vector<std::string> cells;std::string cell;bool quoted=false;
    for(size_t i=0;i<line.size();++i){const char c=line[i];
        if(c=='"'){if(quoted&&i+1<line.size()&&line[i+1]=='"'){cell+='"';++i;}else quoted=!quoted;}
        else if(c==','&&!quoted){cells.push_back(cell);cell.clear();}else if(c!='\r')cell+=c;
    }
    Require(!quoted,"Unterminated CSV field");cells.push_back(cell);return cells;
}
void CsvText(std::ostream& out,const std::string& text){out<<'"';for(char c:text){if(c=='"')out<<'"';out<<c;}out<<'"';}
const std::vector<Column>& Columns(const std::string& stage){return stage=="Estimation"?EstimationColumns():SolverColumns();}
std::string TreeName(const std::string& stage){return stage=="Estimation"?"estimation":"fit";}
ScalarTable ReadCsv(const std::string& path,const std::string& stage){
    std::ifstream in(path);Require(bool(in),"Cannot open stage CSV: "+path);std::string line;Require(bool(std::getline(in,line)),"Missing stage CSV header");
    const auto header=ParseCsv(line);ScalarTable table;table.columns=Columns(stage);std::vector<size_t> indices;
    for(const auto& col:table.columns){const auto it=std::find(header.begin(),header.end(),col.name);
        if(it==header.end()&&stage=="InterpSolver"&&col.name=="particle_hypothesis"){indices.push_back(header.size());continue;}
        Require(it!=header.end(),"Missing scalar column "+col.name);indices.push_back(it-header.begin());}
    while(std::getline(in,line)){if(line.empty())continue;const auto cells=ParseCsv(line);Require(cells.size()==header.size(),"Invalid stage CSV row width");std::vector<Cell> row;
        for(size_t i=0;i<table.columns.size();++i){const std::string value=indices[i]<cells.size()?cells[indices[i]]:"unknown";
            if(table.columns[i].type==ColumnType::Integer)row.emplace_back(std::stoll(value));
            else if(table.columns[i].type==ColumnType::Real)row.emplace_back(value.empty()?std::numeric_limits<double>::quiet_NaN():std::stod(value));
            else row.emplace_back(indices[i]<cells.size()?value:"unknown");
        }
        table.rows.push_back(std::move(row));
    }
    if(!table.rows.empty()){table.info.first_event=table.Integer(0,"event");table.info.last_event=table.info.first_event;table.info.run=table.Integer(0,"orig_run");
        for(size_t i=0;i<table.rows.size();++i){const auto e=table.Integer(i,"event");table.info.first_event=std::min(table.info.first_event,e);table.info.last_event=std::max(table.info.last_event,e);}}
    return table;
}
void WriteCsv(const std::string& path,const ScalarTable& table){
    std::ofstream out(path);Require(bool(out),"Cannot write stage CSV: "+path);out<<std::setprecision(17);
    for(size_t i=0;i<table.columns.size();++i)out<<(i?",":"")<<table.columns[i].name;out<<'\n';
    for(const auto& row:table.rows){Require(row.size()==table.columns.size(),"Scalar row width mismatch");
        for(size_t i=0;i<row.size();++i){if(i)out<<',';
            if(const auto* value=std::get_if<long long>(&row[i]))out<<*value;
            else if(const auto* value=std::get_if<double>(&row[i]))out<<*value;
            else CsvText(out,std::get<std::string>(row[i]));
        }
        out<<'\n';
    }
    out.close();Require(bool(out),"Incomplete stage CSV output");
}
}
size_t ScalarTable::Index(const std::string& name)const{for(size_t i=0;i<columns.size();++i)if(columns[i].name==name)return i;throw std::runtime_error("Missing scalar column "+name);}
double ScalarTable::Number(size_t row,const std::string& name)const{
    const auto& cell=rows.at(row).at(Index(name));if(const auto* value=std::get_if<double>(&cell))return *value;if(const auto* value=std::get_if<long long>(&cell))return static_cast<double>(*value);throw std::runtime_error("Column is not numeric: "+name);
}
long long ScalarTable::Integer(size_t row,const std::string& name)const{
    const auto& cell=rows.at(row).at(Index(name));if(const auto* value=std::get_if<long long>(&cell))return *value;
    throw std::runtime_error("Column is not integer: "+name);
}
const std::vector<Column>& EstimationColumns(){
    static const std::vector<Column> schema=[](){std::vector<Column> s;
        for(const auto* name:{"event","cluster_index","cluster_label","orig_run","orig_event"})s.push_back({name,ColumnType::Integer});
        for(const auto* name:{"ic_amplitude","ic_centroid","ic_integral","ic_multiplicity","vertex_x","vertex_y","vertex_z","center_x","center_y","center_z","polar","azimuthal","brho","dEdx","sqrt_dEdx","dE","arclength"})s.push_back({name,ColumnType::Real});
        s.push_back({"direction",ColumnType::Integer});return s;}();return schema;
}
const std::vector<Column>& SolverColumns(){
    static const std::vector<Column> schema=[](){std::vector<Column> s;
        for(const auto* name:{"event","cluster_index","cluster_label","orig_run","orig_event"})s.push_back({name,ColumnType::Integer});
        for(const auto* name:{"vertex_x","sigma_vx","vertex_y","sigma_vy","vertex_z","sigma_vz","brho","sigma_brho","ke","sigma_ke","polar","sigma_polar","azimuthal","sigma_azimuthal","redchisq"})s.push_back({name,ColumnType::Real});
        s.push_back({"solver_method",ColumnType::Text});s.push_back({"objective_kind",ColumnType::Text});
        s.push_back({"particle_hypothesis",ColumnType::Text});return s;}();return schema;
}
PointcloudData ReadPointcloudStage(const std::string& path){
    if(!Root(path))return ReadPointcloudH5(path);PointcloudData output;Require(ReadPointClouds(path,output.events),"Cannot read Pointcloud ROOT: "+path);
    std::sort(output.events.begin(),output.events.end(),[](const auto&a,const auto&b){return a.event<b.event;});
    for(size_t i=1;i<output.events.size();++i)Require(output.events[i-1].event!=output.events[i].event,"Duplicate pointcloud event");
    TFile file(path.c_str(),"READ");Require(!file.IsZombie(),"Cannot open ROOT file");output.info=ReadInfo(file,output.events.empty()?0:output.events.front().event,output.events.empty()?-1:output.events.back().event);return output;
}
void WritePointcloudStage(const std::string& path,const PointcloudData& data){
    if(!Root(path)){WritePointcloudH5(path,data);return;}Require(WritePointClouds(path,data.events),"Cannot write Pointcloud ROOT");
    TFile file(path.c_str(),"UPDATE");Require(!file.IsZombie(),"Cannot update ROOT metadata");WriteInfo(file,data.info,"Pointcloud");
}
ClusterData ReadClusterStage(const std::string& path){
    if(!Root(path))return ReadClusterH5(path);TFile file(path.c_str(),"READ");Require(!file.IsZombie(),"Cannot open Cluster ROOT");
    auto* tree=dynamic_cast<TTree*>(file.Get("cluster"));Require(tree,"Cluster ROOT lacks cluster tree");
    std::map<long long,ClusterEvent> events;std::map<long long,size_t> declared;
    if(auto* index=dynamic_cast<TTree*>(file.Get("events"))){
        Require(index->GetLeaf("event")&&index->GetLeaf("nclusters"),"Invalid ROOT events tree");
        for(Long64_t i=0;i<index->GetEntries();++i){Require(index->GetEntry(i)>0,"Cannot read ROOT event entry");const auto e=Integer(*index,"event",-1),n=Integer(*index,"nclusters",-1);
            Require(e>=0&&n>=0&&n<=std::numeric_limits<int>::max()&&!events.count(e),"Invalid ROOT event index");ClusterEvent item;item.event=e;item.metadata=ReadEventMetadata(*index,e);events[e]=std::move(item);declared[e]=static_cast<size_t>(n);}
    }
    Require(tree->GetLeaf("event")&&Leaf(*tree,"cluster_index","cluster"),"Cluster ROOT lacks identifiers");
    std::vector<double>* x=nullptr,*y=nullptr,*z=nullptr,*q=nullptr,*scale=nullptr;
    BindVector(*tree,"x",&x,true);BindVector(*tree,"y",&y,true);BindVector(*tree,"z",&z,true);BindVector(*tree,"integral",&q,true);BindVector(*tree,"scale",&scale,false);
    for(Long64_t i=0;i<tree->GetEntries();++i){Require(tree->GetEntry(i)>0,"Cannot read ROOT cluster entry");const auto event=Integer(*tree,"event",-1),index=Integer(*tree,"cluster_index",-1,"cluster");
        Require(event>=0&&index>=0&&index<=std::numeric_limits<int>::max(),"Invalid ROOT cluster identifier");
        auto it=events.find(event);if(it==events.end()){ClusterEvent item;item.event=event;item.metadata=ReadEventMetadata(*tree,event);it=events.emplace(event,std::move(item)).first;}
        auto& item=it->second;Require(index==static_cast<long long>(item.clusters.size()),"ROOT cluster indices must be ordered and contiguous");
        Require(x&&y&&z&&q&&y->size()==x->size()&&z->size()==x->size()&&q->size()==x->size()&&(!scale||scale->size()==x->size()),"Inconsistent ROOT cluster columns");
        Cluster cluster;cluster.label=static_cast<int>(index);cluster.original_label=static_cast<int>(Integer(*tree,"cluster_label",index));
        const auto direction=Integer(*tree,"direction",-1);Require(direction>=-1&&direction<=1,"Invalid ROOT cluster direction");cluster.direction=static_cast<Direction>(direction);
        for(size_t j=0;j<x->size();++j){CloudPoint p;p.x=(*x)[j];p.y=(*y)[j];p.z=(*z)[j];p.integral=(*q)[j];p.scale=scale?(*scale)[j]:1.;
            Require(std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)&&std::isfinite(p.integral)&&std::isfinite(p.scale),"Nonfinite ROOT cluster point");cluster.points.push_back(p);}
        item.clusters.push_back(std::move(cluster));
    }
    tree->ResetBranchAddresses();ClusterData output;
    for(auto& entry:events){if(declared.count(entry.first))Require(entry.second.clusters.size()==declared.at(entry.first),"Declared ROOT cluster count differs");output.events.push_back(std::move(entry.second));}
    output.info=ReadInfo(file,output.events.empty()?0:output.events.front().event,output.events.empty()?-1:output.events.back().event);return output;
}
void WriteClusterStage(const std::string& path,const ClusterData& data){
    if(!Root(path)){WriteClusterH5(path,data);return;}TFile file(path.c_str(),"RECREATE");Require(!file.IsZombie(),"Cannot write Cluster ROOT");
    WriteInfo(file,data.info,"Cluster");TTree events("events","event metadata, including events without clusters");TTree tree("cluster","cleaned clusters");
    Long64_t event=-1,nclusters=0,index=-1,label=-1;Int_t npoints=0,direction=-1;EventMetadata metadata;std::vector<double>x,y,z,q,scale;
    events.Branch("event",&event,"event/L");events.Branch("nclusters",&nclusters,"nclusters/L");BindEventMetadata(events,metadata);
    tree.Branch("event",&event,"event/L");tree.Branch("cluster_index",&index,"cluster_index/L");tree.Branch("cluster_label",&label,"cluster_label/L");
    tree.Branch("npoints",&npoints,"npoints/I");tree.Branch("direction",&direction,"direction/I");BindEventMetadata(tree,metadata);
    tree.Branch("x",&x);tree.Branch("y",&y);tree.Branch("z",&z);tree.Branch("integral",&q);tree.Branch("scale",&scale);
    for(const auto& item:data.events){event=item.event;nclusters=item.clusters.size();metadata=item.metadata;if(metadata.orig_event<0)metadata.orig_event=event;Require(events.Fill()>=0,"Cannot fill ROOT event metadata");
        for(size_t i=0;i<item.clusters.size();++i){const auto& c=item.clusters[i];index=i;label=c.original_label<0?c.label:c.original_label;direction=static_cast<int>(c.direction);npoints=static_cast<int>(c.points.size());
            x.clear();y.clear();z.clear();q.clear();scale.clear();for(const auto& p:c.points){x.push_back(p.x);y.push_back(p.y);z.push_back(p.z);q.push_back(p.integral);scale.push_back(p.scale);}
            Require(tree.Fill()>=0,"Cannot fill ROOT cluster");
        }
    }
    Require(events.Write()>0&&tree.Write()>0,"Cannot write ROOT cluster trees");
}
ScalarTable ReadScalarStage(const std::string& path,const std::string& stage){
    if(Csv(path))return ReadCsv(path,stage);Require(Root(path),"Parquet is handled by the workspace entry point");
    TFile file(path.c_str(),"READ");Require(!file.IsZombie(),"Cannot open scalar ROOT");
    auto* tree=dynamic_cast<TTree*>(file.Get(TreeName(stage).c_str()));if(!tree&&stage=="Estimation")tree=dynamic_cast<TTree*>(file.Get("estimate"));Require(tree,"Missing scalar stage ROOT tree");
    Require(tree->GetLeaf("event")&&Leaf(*tree,"cluster_index","cluster"),"Scalar ROOT lacks identifiers");
    ScalarTable table;table.columns=Columns(stage);std::vector<std::string> text(table.columns.size());std::vector<std::string*> pointers(table.columns.size());
    for(size_t i=0;i<table.columns.size();++i)if(table.columns[i].type==ColumnType::Text){
        if(tree->GetBranch(table.columns[i].name.c_str())){pointers[i]=&text[i];Require(tree->SetBranchAddress(table.columns[i].name.c_str(),&pointers[i])>=0,"Invalid ROOT string column");}
        else if(stage=="InterpSolver"&&table.columns[i].name=="particle_hypothesis")text[i]="unknown";}
    for(Long64_t entry=0;entry<tree->GetEntries();++entry){Require(tree->GetEntry(entry)>0,"Cannot read scalar ROOT row");if(tree->GetLeaf("ok")&&Integer(*tree,"ok",0)!=1)continue;
        const auto event=Integer(*tree,"event",-1),index=Integer(*tree,"cluster_index",-1,"cluster");const auto meta=ReadEventMetadata(*tree,event);std::vector<Cell> row;
        for(size_t i=0;i<table.columns.size();++i){const auto& col=table.columns[i];
            if(col.type==ColumnType::Integer){long long fallback=-1;if(col.name=="cluster_label")fallback=index;else if(col.name=="orig_event")fallback=meta.orig_event;else if(col.name=="orig_run")fallback=meta.orig_run;
                row.emplace_back(Integer(*tree,col.name,fallback,col.name=="cluster_index"?"cluster":""));}
            else if(col.type==ColumnType::Text)row.emplace_back(text[i]);
            else{
                std::string alias;if(col.name=="vertex_x")alias="vx";if(col.name=="vertex_y")alias="vy";if(col.name=="vertex_z")alias="vz";if(col.name=="dEdx")alias="dedx";if(col.name=="dE")alias="charge";if(col.name=="redchisq")alias="cost";
                double value=Real(*tree,col.name,std::numeric_limits<double>::quiet_NaN(),alias);
                if(stage=="InterpSolver"&&col.name.compare(0,7,"vertex_")==0&&!tree->GetLeaf(col.name.c_str())&&Leaf(*tree,alias))value*=.001;
                if(col.name=="sqrt_dEdx"&&!std::isfinite(value))value=std::sqrt(Real(*tree,"dEdx",0,"dedx"));
                if(col.name=="center_z"&&!std::isfinite(value))value=Real(*tree,"vertex_z",0,"vz");
                row.emplace_back(value);
            }
        }
        table.rows.push_back(std::move(row));
    }
    tree->ResetBranchAddresses();table.info=ReadInfo(file,table.rows.empty()?0:table.Integer(0,"event"),table.rows.empty()?-1:table.Integer(table.rows.size()-1,"event"));return table;
}
void WriteScalarStage(const std::string& path,const std::string& stage,const ScalarTable& table){
    if(Csv(path)){WriteCsv(path,table);return;}Require(Root(path),"Parquet is handled by the workspace entry point");
    TFile file(path.c_str(),"RECREATE");Require(!file.IsZombie(),"Cannot write scalar ROOT");WriteInfo(file,table.info,stage);TTree tree(TreeName(stage).c_str(),stage.c_str());
    std::vector<Long64_t> integers(table.columns.size());std::vector<Double_t> reals(table.columns.size());std::vector<std::string> text(table.columns.size());
    for(size_t i=0;i<table.columns.size();++i){const auto& c=table.columns[i];if(c.type==ColumnType::Integer)tree.Branch(c.name.c_str(),&integers[i],(c.name+"/L").c_str());else if(c.type==ColumnType::Real)tree.Branch(c.name.c_str(),&reals[i],(c.name+"/D").c_str());else tree.Branch(c.name.c_str(),&text[i]);}
    for(const auto& row:table.rows){Require(row.size()==table.columns.size(),"Scalar output row width differs");for(size_t i=0;i<row.size();++i){if(table.columns[i].type==ColumnType::Integer)integers[i]=std::get<long long>(row[i]);else if(table.columns[i].type==ColumnType::Real)reals[i]=std::get<double>(row[i]);else text[i]=std::get<std::string>(row[i]);}Require(tree.Fill()>=0,"Cannot fill scalar ROOT");}
    Require(tree.Write()>0,"Cannot write scalar ROOT tree");
}
}
