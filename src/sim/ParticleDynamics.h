#ifndef ParticleDynamics_h
#define ParticleDynamics_h

//----------------------------------------------------------------------
// 粒子在 AT-TPC 中的运动学与轨迹积分
//
// 本类只负责：
//   1. 反应/衰变的两体运动学（InelasticScattering / BreakupReaction 等）
//   2. 运动方程 equationMotion 与 RK4 积分的初态、停止条件
//   3. 靶中能损（catima）
//
// 实验相关的截面、DWBA、拟合等分析函数暂留在这里，
// 后续 M6 会拆分到 analysis/ 目录。
//
// 头文件自洽：本文件自己包含所需的全部头，
// 既能被 ROOT 宏 #include（解释执行），也能被 g++ 编译。
//----------------------------------------------------------------------

// ---- ROOT ----
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TRotation.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

// ---- 物理 ----
#include "catima/catima.h"

// ---- 本工程 ----
#include "../core/MassData.h"
#include "../core/Units.h"
#include "../../extern/rk4_solver/include/rk4_solver.hpp"
#include "TargetDensity.h"

// ---- 标准库 ----
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


//------------------轨迹输出-----------------------//
// 轨迹上的一个点（实验室系）
struct TrackPoint {
    double x = 0.0, y = 0.0, z = 0.0; // m
    double ke = 0.0;                  // 动能 MeV
};

//------------------class define-----------------------//
class ParticleDynamics{

public :

   ParticleDynamics();// 无参构造函数

   void  SetNuclearData();// 设置核反应数据

   double  GetNuclearMass(Int_t A,Int_t Z);// 获取原子核质量 (u，已扣除电子)

   TString   GetNuclearName(Int_t A,Int_t Z);// 获取原子核名称

   void equationMotion(double t, const double (&state)[6], double (&results)[6]);// 微分方程

   void initialState(
      double kinetic_energy_mev,      // 动能 (MeV)
      double x0, double y0, double z0, // 初始位置 (m)
      double theta, double phi,       // theta/phi 均为弧度；phi 是 XY 平面方位角
      double state[6]                 // 输出: [x,y,z,gvx,gvy,gvz]
   );// 设置初始状态

   bool StopEvent(const double t, const double (&x)[6], double (&x_plus)[6]);// 事件检测函数

   // 重置单次积分的状态（每次新事件前必须调用）
   void SetInitialConditions(){this->StopCondition = false;};
   void Reset(){this->StopCondition = false;};

   // ---- 停止条件（几何/能量，单位 m 与 MeV）----
   void SetStopRho(double rho_m){ this->rho_max = rho_m; }          // 默认 0.28 m
   void SetStopZ(double z_min_m, double z_max_m){ this->z_min = z_min_m; this->z_max = z_max_m; }
   void SetStopEnergy(double ke_mev){ this->ke_min_stop = ke_mev; } // 默认 1e-5 MeV 约 10 eV

   // ---- 探测器模拟用接口 ----
   // 设置被积分的粒子种类（会重建 catima 投影物）
   void SetParticle(int A, int Z);

   // 从顶点出发积分一条轨迹；theta/phi 输入单位 rad，dt 单位 s
   void Propagate(double ke_mev, const double vtx[3], double theta, double phi,
                  std::vector<TrackPoint>& track, double dt = 1.0e-10,
                  int max_steps = 200000);

   // 按靶材配置重建材料（stoich 为化学计量数，density 单位 g/cm^3）
   void ConfigureTarget(int A, int Z, double stoich, double density_gcm3, double length_cm = 100.0);
   void ConfigureTarget(const std::vector<Element>& compound, double density_gcm3,
                        double length_cm = 100.0);
   // 直接使用配置层构造的靶对象；GasTarget/GasMixTarget 会负责由压力计算密度。
   void ConfigureTarget(const Target& target, double length_cm = 100.0);

   // 设置场内场强
   void ConfigureFields(double bfield_t, double efield_vm) { Bfield = bfield_t; Efield = efield_vm; }

   // dE/dx 缩放因子（默认 1）。给 Kalman 滤波用：
   //   > 0：正常能损，σ 点按能损噪声缩放；
   //   < 0：反向回推（能量增加），用于把径迹外推到顶点。
   void SetELossScale(double s) { eloss_scale = s; }
   double GetELossScale() const { return eloss_scale; }
   double GetDensity() const { return density; }        // g/cm^3
   double GetMaterialDensity() const { return density; }

   void dataAnalysis();

   void loopTree(TTree* tree , double fitLeft=-2., double fitRight=2., bool reset=true);

   void test();

   void SRIMtest();

   TGraphErrors* GetCrossSection(TH1D *hTheta ,TGraphErrors* grCrossSection);////input: Be-12 hTheta  output TH1D *hCrossSection

   void GetEfficiency(TH1D *hTheta , TString filename="Efficiency.txt");// 计算重建效率

   void SetEfficiency(TGraphErrors* grCrossSection, TString filename="Efficiency.txt");// 计算重建效率

   TGraphErrors* CutCrossSectionData(TGraphErrors* grCrossSection, double Left=20., double Right=180.);// Cut CrossSectionData

   double GetEnergyLoss(double distance);//input: cm  output: MeV

   // ke->MeV theta->radian
   double GetExciteEnergy(double projectile_ke, double ejectile_ke, double theta);

   // 质心系转换 ke->MeV theta->radian  out: theta->degree
   double AngleTrans(double projectile_ke,  double ejectile_ke, double theta, double target_ke=0.);// 质心系转换

   TGraph* GetRefLine(TGraph* T, double exciteEnergy);// calculate Excite Energy Reference Line

   TGraph* GetRefLine(TGraph* T, TGraph* Tcm, double exciteEnergy);// Function overload

   //目前只适用于弹散/非弹散
   TGraph* GetRefLineRe(TGraph* T, TGraph* Tcm, double exciteEnergy);// Function overload 相对论下质心系转换, mising mass 轻核

   //目前只适用于弹散/非弹散
   TGraph* GetRefLineS(TGraph* T, TGraph* T2, double exciteEnergy);// Function overload 相对论下质心系转换  两个核

   // 破裂反应 运动学曲线
   TGraph* GetBreakupRefLineS(TGraph* T, TGraph* T2,double Ek, double exciteEnergy);// Function overload 相对论下质心系转换  两个核

   TGraph* computeDiff(const TGraph* gr);// 计算微分

   TGraph* ReadDWBA(TGraph* grDwBA, TString filename="DWBA.Xsec.txt");//

   void WriteCrossSection(TGraphErrors* grCrossSection, TString filename="CrossSection.gs.txt");// 输出微分截面文件

   double Normalize(TGraph* theory, TGraph* experiment,Int_t start=20,Int_t end=40);// 拟合 获取归一化刻度系数

   double GausIntegral(double A, double mu, double sigma, double x1, double x2, double binWidth=1.);// 计算给定高斯函数参数的区间积分

   void minsGaus(TH1D *hTheta , double x1, double x2, TString filename="minsGaus.txt");// 将参数减去给定高斯函数参数的区间积分

   void SetGausPar(TF1* fgaus);// 对双高斯函数的参数设定拟合初值

   void AnaThetaCM(TTree* tree, double Emean=2.1,TString filename="AnaThetaCM");// 分质心系角度的数据拟合、激发态参考能量Emean MeV

   TH1D* AnaThetaCrossSection(TTree* tree,double Emean=2.1, TString filename="AnaThetaCM");// 分质心系角度计算微分截面、激发态参考能量Emean MeV

   // 考虑分角度后的计数还有另一个激发态导致的pileup，计算考虑pileup的微分截面----->与 minsGaus()函数功能相近
   TH1D* minsPileUp(TH1D *hTheta , double Emean=0.,TString Emeanfilename="AnaThetaCM", TString filename="AnaThetaCM");// 此处的filename应为另一个峰的拟合文件

   void  minsPileUp(TH1D *hexPerBin,double bin_center_cm, double &minValue, Int_t binValue, Int_t peak);// 函数重载

   static double Response(double *x, double *par);// DWBA数据的响应函数

   static double twoGaus(double *x, double *par);// 双高斯拟合函数

   static double bgOnly(double *x, double *par); // 本底函数拟合

   static double combinedFit(double *x, double *par); //复合拟合函数

   void SetBgPar(TF1* fbg);// 对本底函数的参数设定拟合初值

   bool OpenFile(std::ifstream &inf, TString filename);// 打开文件

   void angleResolution();// 角度分辨引入的能谱展宽 mising mass

   void targetDeposit();// 引入的能谱展宽

   // 获取计算DWBA数据的概率分布函数,gr运动学参考线,thetaL,thetaR为实验室系计算角度范围
   double GetDWBApdf(TGraph* grDwBA, TGraph* gr, double powerTheta, double thetaL=58., double thetaR=82);

   TTree* GetTree(){return this->tree;}// 获取tree的指针

   double GetProjectileMass(){return this->projectile_mass;}// nucleus mass : u

   double GetEjectileMass(){return this->ejectile_mass;}// nucleus mass : u

   double GetTargetMass(){return this->target_mass;}// nucleus mass : u

   double GetResidualMass(){return this->residual_mass;}// nucleus mass : u

   double GetProjectileMassMeV(){return this->projectile_mass*u_MeV;}// nucleus mass : u

   double GetEjectileMassMeV(){return this->ejectile_mass*u_MeV;}// nucleus mass : u

   double GetTargetMassMeV(){return this->target_mass*u_MeV;}// nucleus mass : u

   double GetResidualMassMeV(){return this->residual_mass*u_MeV;}// nucleus mass : u

   double GetInitialBeamEnergy(){return this->proj_energy_start;}// MeV

   double GetBeamEnergy(double distance=0);//  [in]:cm [out]:MeV

   TVector3 ParticleToWorld(const TVector3& vec_particle, double theta, double phi);// 旋转矩阵函数：将粒子坐标系中的向量转换到世界坐标系

   void CalclatePID(Int_t m_A=12 , Int_t m_Z=4);// 计算能损与粒子鉴别

   catima::Material GetTargetMaterial(){return this->material;}// 建议引用接收使用

   catima::Projectile GetProjectile(){return this->projectile;}// 建议引用接收使用

    void InelasticScattering(
      double &exit_energy,
      double &recoil_energy,
      double &exit_angle, //Lab frame
      double &recoil_angle, //Lab frame
      double beam_energy,
      double excited_energy,
      double theta
   );// simulate process of inelastic scattering


   void BreakupReaction(
      double parent_mass,
      double fragment1_mass,
      double fragment2_mass,
      double parent_Ek,
      double parent_Ex,
      double theta,
      double &fragment1_energy,
      double &fragment2_energy,
      double &fragment1_angle,
      double &fragment2_angle
   );// simulate process of breakup reaction

public :

   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   TF1* gaus;

   static TGraph* grDWBA;// DWBA计算结果
   static TGraph* grDWBA2;// DWBA计算结果
   static TGraph* grDWBA3;// DWBA计算结果
   static Int_t RangeStart;// DWBA数据拟合区间
   static Int_t RangeEnd;

   TRandom3* r;

   // 定义一个 map 存储角度区间对应的多个系数
   static std::map<double, std::vector<double>> bgCoeffMap;


private :

   // Stopping data loaded from external file (energies in MeV, dedx in MeV/mm)
   std::vector<double> stop_energies_mev;
   std::vector<double> stop_dedx_mev_per_mm;

   // Load stopping data file (returns true on success)
   bool LoadStoppingFile(const char* filename);

   // Interpolate dedx from loaded stopping data. Input energy in MeV.
   // Returns dedx in MeV/mm. If no data loaded, returns 0.
   double DedxFromFile(double energy_mev) const;



   TTree* tree;

   // Specify your nuclei  定义变量  内联初始化：C++11

   // The nucleus we observe (the one we fitted)      m_b
   Int_t ejectile_z = 1;
   Int_t ejectile_a = 1;
   double ejectile_mass;// nucleus mass : u

   // The incoming nucleus (the beam)                 m_a
   Int_t projectile_z= 4;
   Int_t projectile_a= 12;
   Int_t particle_in=45280*100*10;// IC Counts    IC1==45888、IC2==51639 IC版本18=45280
   double projectile_mass;// nucleus mass : u

   // The target nucleus                              m_A
   Int_t target_z = 1;
   Int_t target_a = 1;
   double density=6.6168e-05;// material density  in g/cm3 600Torr-》6.6168e-05
   double target_mass;// nucleus mass : u

   // The residual nucleus is equal to projectile     m_B                             m_A
   Int_t residual_z = 4;
   Int_t residual_a = 12;
   double residual_mass;// nucleus mass : u

   // The fragment1 nucleus    破裂反应碎片1  重核
   Int_t fragment1_z = 6;
   Int_t fragment1_a = 10;
   double fragment1_mass;// nucleus mass : u

   // The fragment2 nucleus     破裂反应碎片2
   Int_t fragment2_z = 2;
   Int_t fragment2_a = 4;
   double fragment2_mass;// nucleus mass : u

   // The fragment3 nucleus     破裂反应碎片3
   Int_t fragment3_z = 5;
   Int_t fragment3_a = 10;
   double fragment3_mass;// nucleus mass : u

   // catima information (calculate energy loss program)
   catima::Material material;      // define target material
   catima::Projectile projectile; // define projectile

   // EoM 专用缓存：避免在积分内层循环里反复构造 catima 对象
   catima::Projectile proj_eom;

   // Initial beam energy
   double proj_energy=12.8;    // MeV/u
   double proj_energy_start; // MeV
   double proj_energy_stop;  // MeV

   double Bfield=2.85;           // Tesla
   double Efield=60000;//45000;       // V/m

   bool StopCondition=false; // 当变量为true时，积分器将停止积分

   // ---- 停止条件（可配置）----
   double rho_max = Units::at_tpc_rho_max;   // m
   double z_max   = Units::at_tpc_length;    // m
   double z_min   = 0.0;                     // m
   double ke_min_stop = 1.0e-5;              // MeV (约 10 eV)
   double eloss_scale = 1.0;                 // dE/dx 缩放（Kalman 用）

   // 单位常数：统一指向 core/Units.h，避免多处重复定义
   static constexpr double me           = Units::me_u;
   static constexpr double u_MeV        = Units::u_MeV;
   static constexpr double C            = Units::C;
   static constexpr double E_CHARGE     = Units::E_CHARGE;
   static constexpr double MEV_2_KG     = Units::MEV_2_KG;
   static constexpr double MEV_2_JOULE  = Units::MEV_2_JOULE;

      // Declaration of leaf types   ->  read root file
   Long64_t        event;
   Long64_t        cluster_index;
   Long64_t        cluster_label;
   Double_t        vertex_x;
   Double_t        sigma_vx;
   Double_t        vertex_y;
   Double_t        sigma_vy;
   Double_t        vertex_z;
   Double_t        sigma_vz;
   Double_t        brho;
   Double_t        sigma_brho;
   Double_t        ke;
   Double_t        sigma_ke;
   Double_t        polar;
   Double_t        sigma_polar;
   Double_t        azimuthal;
   Double_t        sigma_azimuthal;
   Double_t        redchisq;
   Double_t        excitation;
   Double_t        theta_cm;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_cluster_index;   //!
   TBranch        *b_cluster_label;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_sigma_vx;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_sigma_vy;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_sigma_vz;   //!
   TBranch        *b_brho;   //!
   TBranch        *b_sigma_brho;   //!
   TBranch        *b_ke;   //!
   TBranch        *b_sigma_ke;   //!
   TBranch        *b_polar;   //!
   TBranch        *b_sigma_polar;   //!
   TBranch        *b_azimuthal;   //!
   TBranch        *b_sigma_azimuthal;   //!
   TBranch        *b_redchisq;   //!
   TBranch        *b_excitation;   //!
   TBranch        *b_theta_cm;   //!


};// Class end here!

#endif

