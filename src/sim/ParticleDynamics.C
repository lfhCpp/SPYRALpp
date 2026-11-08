#include "ParticleDynamics.h"

#include "../core/DataPath.h"

#include <mutex>

// 本文件既可被 ROOT 宏 #include（解释执行），也可单独编译。
// using 只作用于本翻译单元；所有头文件仍显式写 std::。
using namespace std;

// ROOT 的对象注册表（gROOT 的列表）不是线程安全的：多线程同时 new TRandom3
// 会踩堆。这里只把"创建 ROOT 对象"这一步串行化。
static std::mutex& ROOTObjectMutex() {
    static std::mutex m;
    return m;
}

// catima 的 Material/Projectile 构造会初始化它的全局存储（核素表、spline 缓存），
// 这部分不是线程安全的。构造/配置串行化，积分（equationMotion）保持并行。
static std::mutex& CatimaSetupMutex() {
    static std::mutex m;
    return m;
}


//---------------------------- static define-----------------------//

 TGraph* ParticleDynamics::grDWBA = nullptr; 
 TGraph* ParticleDynamics::grDWBA2 = nullptr; 
 TGraph* ParticleDynamics::grDWBA3 = nullptr; 
 Int_t ParticleDynamics::RangeStart=10;// DWBA数据拟合区间 10-50
 Int_t ParticleDynamics::RangeEnd=140;
 std::map<double, std::vector<double>> ParticleDynamics::bgCoeffMap;

//--------------------->main function start<-----------------------//
ParticleDynamics::ParticleDynamics(){ 
    // 整个构造过程串行化（catima 全局初始化 + ROOT 对象创建都不是线程安全的）
    std::lock_guard<std::mutex> catima_lock(CatimaSetupMutex());
    fChain=nullptr;
    
    tree=nullptr;

    {
        std::lock_guard<std::mutex> lock(ROOTObjectMutex());
        r = new TRandom3(0);
    }

    // TF1 / TGraph / gStyle 会注册到 ROOT 的全局列表，构造过程不是线程安全的。
    // 这些对象只服务于遗留分析代码，放到守卫里，保证模拟路径可以多线程调用。
#ifdef RK4_WITH_LEGACY_ANALYSIS
    ParticleDynamics::grDWBA = new TGraph;
    ParticleDynamics::grDWBA2 = new TGraph;
    ParticleDynamics::grDWBA3 = new TGraph;

    gStyle->SetOptFit(1);
    gStyle->SetTitleFontSize(0.05); // 设置标题字体大小
    gStyle->SetTitleAlign(23); // 居中对齐

    gaus=new TF1("gaus","gaus",-10.,20.);
#endif
    
    //Declaration();

    SetNuclearData();
    // try load external stopping data for comparison
    // 数据文件统一放在 data/ 下，用 DataPath 定位（兼容源码树与安装后的 share 目录）
    if (!this->LoadStoppingFile(rk4::DataFile("Hydrogen in Hydrogen (gas).txt").c_str())) {
        if (!this->LoadStoppingFile("data/Hydrogen in Hydrogen (gas).txt")) {
            this->LoadStoppingFile("../data/Hydrogen in Hydrogen (gas).txt");
        }
    }
    
#ifdef RK4_WITH_LEGACY_ANALYSIS
    test();
#endif

    //angleResolution();

    //targetDeposit();

    //dataAnalysis();

    //CalclatePID();

  }

#ifdef RK4_WITH_LEGACY_ANALYSIS
void  ParticleDynamics::test(){

   
    
}

void  ParticleDynamics::SRIMtest(){


    catima::Config cfg_bare;
    cfg_bare.corrections |= catima::corrections::no_barkas;
    cfg_bare.corrections |= catima::corrections::no_lindhard;
    cfg_bare.corrections |= catima::corrections::no_shell_correction;
    

    double mass_kg = this->ejectile_mass*u_MeV * MEV_2_KG;
    // 获取能量损失并计算减速度
    catima::Material D2(ejectile_a, ejectile_z, 2); //19.2 eV
    D2.density(this->density).thickness(100*this->density);//  g/cm3 and g/cm2
    catima::Projectile ejectile(ejectile_a, ejectile_z);
    ejectile.T =11./this->ejectile_mass; // setting energy to  MeV/u;
    D2.I(50);
    //cout<<"D2.I():"<<D2.I()<<endl;

    double dedx = catima::dedx(ejectile, D2);; // MeV/g/cm2
    // eloss_scale 默认 1（乘 1.0 是精确的，不影响仿真结果）；
    // Kalman 滤波用 >0 的缩放描述能损噪声、用 <0 做反向回推。
    double deceleration =
        (dedx * MEV_2_JOULE * this->density * 100.0 * this->eloss_scale) / mass_kg;//F = dE/dx  dE=Fdx  density:g/cm^3

    cout<<"T:"<<ejectile.T*ejectile_mass<<"\tdedx:"<<dedx*this->density/10.<<endl;

    // Also get dedx from external SRIM file (if loaded). File provides dE/dx in MeV/mm.
    double dedx_file_mev_per_mm = this->DedxFromFile(ejectile.T*ejectile_mass); // MeV/mm
    double deceleration_file = 0.0;
    if(dedx_file_mev_per_mm > 0.0){
        // convert MeV/mm to MeV/m (x1000), then to J/m and to acceleration
        deceleration_file = (dedx_file_mev_per_mm * 1000.0 * MEV_2_JOULE) / mass_kg; // m/s^2
    }

    // Print a single comparison at t near zero for quick check
    
    std::cout << "[Energy] " << ejectile.T*ejectile_mass << "MeV "<<std::endl<<"catima dedx(MeV/mm):\t" 
                  << dedx*this->density/10.<<std::endl<< "file dedx(MeV/mm):\t" << dedx_file_mev_per_mm
                  << std::endl;
   
    
}

// Load stopping-file produced by SRIM (expects energy then unit then dE/dx columns)
#endif // RK4_WITH_LEGACY_ANALYSIS
bool ParticleDynamics::LoadStoppingFile(const char* filename) {
    stop_energies_mev.clear();
    stop_dedx_mev_per_mm.clear();

    std::ifstream inf(filename);
    if(!inf.is_open()) return false;

    std::string line;
    while(std::getline(inf, line)){
        // skip empty or header lines that don't start with a number
        std::istringstream iss(line);
        double energyVal;
        if(!(iss >> energyVal)) continue;
        std::string unit;
        if(!(iss >> unit)) continue;

        // read two dedx columns if present (elec and nuclear), else single
        double d1=0.0, d2=0.0;
        if(!(iss >> d1)) continue;
        iss >> d2; // if fails, d2 stays 0

        double energy_mev = energyVal;
        // handle unit like "keV" or "MeV"
        if(unit.find("keV") != std::string::npos) energy_mev = energyVal / 1000.0;
        // total stopping (elec + nuclear) reported in file (units shown as MeV/mm)
        double dedx_mev_per_mm = d1 + d2;

        stop_energies_mev.push_back(energy_mev);
        stop_dedx_mev_per_mm.push_back(dedx_mev_per_mm);
    }

    // ensure monotonic energy ordering
    if(stop_energies_mev.size() > 1){
        std::vector<std::pair<double,double>> pairs;
        for(size_t i=0;i<stop_energies_mev.size();++i) pairs.emplace_back(stop_energies_mev[i], stop_dedx_mev_per_mm[i]);
        std::sort(pairs.begin(), pairs.end());
        for(size_t i=0;i<pairs.size();++i){ stop_energies_mev[i]=pairs[i].first; stop_dedx_mev_per_mm[i]=pairs[i].second; }
    }

    return !stop_energies_mev.empty();
}

// Linear interpolation of dedx (energy in MeV) -> returns dedx in MeV/mm
double ParticleDynamics::DedxFromFile(double energy_mev) const {
    if(stop_energies_mev.empty()) return 0.0;
    // clamp
    if(energy_mev <= stop_energies_mev.front()) return stop_dedx_mev_per_mm.front();
    if(energy_mev >= stop_energies_mev.back()) return stop_dedx_mev_per_mm.back();

    auto it = std::lower_bound(stop_energies_mev.begin(), stop_energies_mev.end(), energy_mev);
    size_t idx = std::distance(stop_energies_mev.begin(), it);
    if(idx == 0) return stop_dedx_mev_per_mm[0];
    double x1 = stop_energies_mev[idx-1];
    double x2 = stop_energies_mev[idx];
    double y1 = stop_dedx_mev_per_mm[idx-1];
    double y2 = stop_dedx_mev_per_mm[idx];
    double t = (energy_mev - x1) / (x2 - x1);
    return y1 + t*(y2 - y1);
}

// 运动方程函数 定义微分方程
void ParticleDynamics::equationMotion(
    double t,
    const double (&state)[6], 
    double (&results)[6]
) {

    // State = [x, y, z, vx, vy, vz] m
    // Derivative = [vx, vy, vz, ax, ay, az] (returns)
    // 状态向量 = [x, y, z, γβ_x, γβ_y, γβ_z]
    // 导数向量 = [v_x, v_y, v_z, d(γβ_x)/dt, d(γβ_y)/dt, d(γβ_z)/dt]

    // 计算 gv = sqrt(gvx² + gvy² + gvz²) ==》γβ
    double gv = std::sqrt(state[3] * state[3] + 
                         state[4] * state[4] + 
                         state[5] * state[5]);
    
    // 计算 beta 和 gamma
    double beta = std::sqrt(gv * gv / (1.0 + gv * gv)); // γ = √(1 + (γβ)²) = √(1 + gv²)
    double gamma = gv / beta;                           // β = (γβ)/γ = gv/γ
    
    // 计算单位向量与速度 (m/s)
    // 用定长数组，避免在积分内层循环里做堆分配
    double unit_vector[3];
    double velo[3];
    for (int i = 0; i < 3; ++i) {
        unit_vector[i] = state[3 + i] / gv;
        velo[i] = unit_vector[i] * beta * C;
    }
    
    // 计算动能 (MeV)
    double kinetic_energy = this->ejectile_mass*u_MeV * (gamma - 1.0);// ejectile.mass：粒子的静止质量（MeV/c²）
    
    // 低能停止统一交给 StopEvent() 判断，本函数保持无副作用

    // 计算电荷和质量
    double charge_c = ejectile_z * E_CHARGE;
    double mass_kg = this->ejectile_mass*u_MeV * MEV_2_KG;
    double q_m = charge_c / mass_kg;
    
    // 获取能量损失并计算减速度
    // 注意：不再每次构造 catima 对象。RK4 每步调用 4 次，
    // 原来的写法会在 10^5 步的轨迹上产生 4x10^5 次堆分配。
    proj_eom.T = kinetic_energy / this->ejectile_mass; // MeV/u
    double dedx = catima::dedx(proj_eom, this->material); // MeV/(g/cm2)
    //double dedx =catima::calculate(ejectile,D2,kinetic_energy/this->ejectile_mass).dEdxi; // MeV/g/cm2
    //cout<<"kinetic_energy:"<<kinetic_energy<<"  dedx:"<<dedx<<endl;
    double deceleration = (dedx * MEV_2_JOULE * this->density * 100.0) / mass_kg;//F = dE/dx  dE=Fdx  density:g/cm^3

    // // Also get dedx from external SRIM file (if loaded). File provides dE/dx in MeV/mm.
    // double dedx_file_mev_per_mm = this->DedxFromFile(kinetic_energy); // MeV/mm
    // double deceleration_file = 0.0;
    // if(dedx_file_mev_per_mm > 0.0){
    //     // convert MeV/mm to MeV/m (x1000), then to J/m and to acceleration
    //     deceleration_file = (dedx_file_mev_per_mm * 1000.0 * MEV_2_JOULE) / mass_kg; // m/s^2
    // }

    // // Print a single comparison at t near zero for quick check
    // if(t < 1e-9){
    //     std::cout << "[Energy] " << kinetic_energy << " MeV | catima dedx(MeV/g/cm2): " << dedx
    //               << " -> decel: " << deceleration
    //               << " m/s^2 | file dedx(MeV/mm): " << dedx_file_mev_per_mm
    //               << " -> decel_file: " << deceleration_file << " m/s^2" << std::endl;
    // }
  
    
    // 计算导数:
    
    // 位置导数 = 速度
    results[0] = velo[0];
    results[1] = velo[1];
    results[2] = velo[2];
    
    //  速度导数 = 加速度（考虑洛伦兹力和能量损失） 相对论力学下
    //  F = d(γ m0 v)/dt = m0 d(γ v)/dt
    //  导数需要转换为d(γβ)/dt ==》γβ = (γv)/c
    results[3] = (q_m * velo[1] * Bfield - deceleration * unit_vector[0]) / C;
    results[4] = (q_m * (-1.0 * velo[0] * Bfield) - deceleration * unit_vector[1]) / C;
    results[5] = (q_m * Efield - deceleration * unit_vector[2]) / C;// 注意探测器的摆放决定了E场方向
    
}

// 计算gv格式的初始状态 (gamma * v)
void ParticleDynamics::initialState(
    double kinetic_energy_mev,      // 动能 (MeV)
    double x0, double y0, double z0, // 初始位置 (m)
    double theta, double phi,       // theta/phi 均为弧度
    double state[6]                 // 输出: [x,y,z,gvx,gvy,gvz]
) {
    // 1. 设置初始位置
    state[0] = x0;
    state[1] = y0;
    state[2] = z0;
    
    // 2. 计算总能量和gamma因子
    double mass_mev = ejectile_mass*u_MeV;  // 将amu转换为MeV/c²
    double total_energy = kinetic_energy_mev + mass_mev;  // 总能量 (MeV)
    double gamma = total_energy / mass_mev;      // 洛伦兹因子
    
    // 3. 计算速度大小
    double momentum = sqrt(total_energy * total_energy - mass_mev * mass_mev);
    double beta = momentum / total_energy;  // v/c
    double speed = beta * C;                // 速度大小 (m/s)
    
    // 4. 计算速度分量

    double vx = speed * sin(theta) * cos(phi);
    double vy = speed * sin(theta) * sin(phi);
    double vz = speed * cos(theta);
    
    // 5. 转换为gv格式
    state[3] = gamma * vx / C;  // gvx
    state[4] = gamma * vy / C;  // gvy
    state[5] = gamma * vz / C;  // gvz
    //cout<<"Initial gv: "<<state[3]<<" , "<<state[4]<<" , "<<state[5]<<endl;
}

bool ParticleDynamics::StopEvent(const double t, const double (&x)[6], double (&x_plus)[6]) {

    // 位置边界：灵敏体积是半径 rho_max 的圆柱，z 属于 [z_min, z_max]
    const double rho = std::sqrt(x[0]*x[0] + x[1]*x[1]);
    if (rho >= this->rho_max || x[2] >= this->z_max || x[2] < this->z_min) {
        this->StopCondition = true;
    }

    // 能量边界：动能低于阈值即认为粒子已停止
    // （原实现写在 equationMotion 里，让运动方程带副作用，也无法单独测试）
    if (!this->StopCondition) {
        const double gv = std::sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
        if (gv <= 0.0) {
            this->StopCondition = true;
        } else {
            const double beta  = std::sqrt(gv*gv / (1.0 + gv*gv));
            const double gamma = gv / beta;
            const double kinetic_energy = this->ejectile_mass * u_MeV * (gamma - 1.0);
            if (kinetic_energy < this->ke_min_stop) {
                this->StopCondition = true;
            }
        }
    }

    if (this->StopCondition) {
        for (int i = 0; i < 6; ++i) {
            x_plus[i] = x[i];
        }
    }

    return this->StopCondition;
}

#ifdef RK4_WITH_LEGACY_ANALYSIS
void ParticleDynamics::dataAnalysis(){// fun
    
    // step 1：卡Z轴坐标，排除异常坐标事件
    TCut cutZ = "vertex_z > 0.005 && vertex_z < 0.995";
    TCut cutZ2 = "vertex_z > 0.005 && vertex_z < 0.5";
    this->tree= fChain->CopyTree(cutZ);  // 应用第一个 cut 条件  单位：m
    this->tree= tree->CopyTree(cutZ2);  // 应用第一个 cut 条件  单位：m
    tree= tree->CopyTree("!(3*0.1775 < excitation && excitation < 2.1-3*0.1775 && ke<8 )"); // sigma 基态：0.1775

    TH1D *hKe=new TH1D("Ke","",500,135,160);// projectile_ke : MeV
    TH1D *hdistance=new TH1D("hdistance","",500,0,105);    
    TH1D *hex=new TH1D("hex","",150,-1.,8);// excitation : MeV
    TH2D *hKL=new TH2D("hKL","",1000,0,180,1000,0,50);
    TH2D *hKcm=new TH2D("hKcm","",1000,0,180,1000,0,50);
    TH2D *hEZ=new TH2D("hEZ","",150,-1.,8,100,0,1.);//激发能与Z位置
    TH1D *hBepolar=new TH1D("hBepolar","",36,0,180);// 设置截面bin的大小
    TH1D *hBepolarAdd=new TH1D("hBepolarAdd","",36,0,180);// 对分角度得到的截面进行加和
    TH1D *hBepolar_bgsub   = new TH1D("hBepolar_bgsub","",18,0,180);// 
    
    // 运动学参考线
    TGraph *gr = new TGraph;      // GetRefLine
    TGraph *grex1 = new TGraph;   // GetRefLine
    TGraph *grex2 = new TGraph;   // GetRefLine
    TGraph *grex3 = new TGraph;   // GetRefLine
    TGraph *grex4 = new TGraph;   // GetRefLine
    TGraph *grex5 = new TGraph;   // GetRefLine
    TGraph *grCM = new TGraph;    // GetRefLine
    TGraph *grCMex1 = new TGraph; // GetRefLine
    TGraph *grCMex2 = new TGraph; // GetRefLine
    TGraph *grCMex3 = new TGraph; // GetRefLine
    TGraph *grCMex4 = new TGraph; // GetRefLine
    TGraph *grCMex5 = new TGraph; // GetRefLine
    this->GetRefLine(gr,grCM,0.);
    this->GetRefLine(grex1,grCMex1,2.109);
    this->GetRefLine(grex2,grCMex2,2.715);
    this->GetRefLine(grex3,grCMex3,4.412);
    this->GetRefLine(grex4,grCMex4,5.724);
    this->GetRefLine(grex5,grCMex5,6.6);

    // step 2：循环事件，得到反应点分布图、出射粒子角度-动能图、Be-12激发能谱图
    TCanvas* can=new TCanvas("can", "EnergyLoss", 1000, 600);
     
    
    this->loopTree(tree);
    can->Print("12BeAnalysis.pdf(");

    // step 3：选择CUTG、再次循环事件
    TCanvas* c1=new TCanvas("c1", "", 1000, 600);
    c1->Divide(2,1);
    c1->cd(1);  
    tree->Draw("ke:180*polar/TMath::Pi()>>h1(500,0,90,500,0,50)","","colz"); 

    // 获取直方图对象
    TH2D *h1 = (TH2D*)gPad->GetPrimitive("h1");

    // 设置坐标轴标题h 
    h1->GetXaxis()->SetTitle("#theta / degree");
    h1->GetYaxis()->SetTitle("Kinetic Energy / MeV");
    h1->GetXaxis()->CenterTitle(); 
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->SetNdivisions(5,5,0); 
    h1->GetYaxis()->SetNdivisions(5,5,0); 
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleOffset(0.90);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(0.9);
    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);
    // 更新画布
    gPad->Modified();
    gPad->Update();
    
    
    //gROOT->LoadMacro("CUTG0.C");// Load and execute the ROOT macro
    gr->Draw("L"); // "L" 表示画线 
    grex1->Draw("L"); // "L" 表示画线 
    grex2->Draw("L"); // "L" 表示画线 
    grex3->Draw("L"); // "L" 表示画线 
    grex4->Draw("L"); // "L" 表示画线 
    grex5->Draw("L"); // "L" 表示画线 
    c1->cd(2);
    hKcm->Draw("colz");
    grCM->Draw("L"); // "L" 表示画线
    grCMex1->Draw("L"); // "L" 表示画线
    grCMex2->Draw("L"); // "L" 表示画线
    grCMex3->Draw("L"); // "L" 表示画线
    grCMex4->Draw("L"); // "L" 表示画线
    grCMex5->Draw("L"); // "L" 表示画线
    c1->Print("12BeAnalysis.pdf");


    //*********************************分角度分析数据************************************************//
    
    // 分角度拟合  获取数据 
    //00  bin5: 2.109、2.715、4.412、5.724、6.6 ？---》2.012 2.6、4.532、5.75、6.6
    
    //this->AnaThetaCM(tree,2.012,"AnaThetaCMs_test");// 生成分角度分析文本文件 

    // 分角度分析截面
    //hBepolar=this->AnaThetaCrossSection(tree,2.012,"AnaThetaCM21");// 生成分角度分析截面，需注意loopTree的设置
    //hBepolar=this->AnaThetaCrossSection(tree,2.6,"AnaThetaCMsub27");
    //hBepolar=this->AnaThetaCrossSection(tree,4.532,"AnaThetaCM45");
    //hBepolar=this->AnaThetaCrossSection(tree,5.75,"AnaThetaCM57");
    //hBepolar=this->AnaThetaCrossSection(tree,6.6,"AnaThetaCM66");

    // 分角度减去pileup影响
    //hBepolar=this->minsPileUp(hBepolar,2.012,"AnaThetaCM21","AnaThetaCMsub27");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,2.6,"AnaThetaCMsub27","AnaThetaCMsub21");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,4.532,"AnaThetaCM45","AnaThetaCM57");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,4.532,"AnaThetaCM45","AnaThetaCMsub27");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,5.75,"AnaThetaCM57","AnaThetaCM45");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,5.75,"AnaThetaCM57","AnaThetaCM66");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,6.6,"AnaThetaCM66","AnaThetaCM45");// 文件名：分析的激发态  减去pileup的激发态
    //hBepolar=this->minsPileUp(hBepolar,6.6,"AnaThetaCM66","AnaThetaCM57");// 文件名：分析的激发态  减去pileup的激发态
    //*********************************分角度分析数据************************************************//


    //******************************************自由CUT分析**************************************//
    //应用第二个 cut 条件 ---->使用运动学、激发能、质心系角度 cut  .q
    double peak_low  = 4.532-3*0.1775;
    double peak_high = 4.532+3*0.1775;
    //tree= tree->CopyTree("CUTG");  // 运动学 cut
    //tree= tree->CopyTree("ke<18");  
    //tree= tree->CopyTree("-6*0.1775 < excitation && excitation < 6*0.1775");  // sigma 基态：0.1775
    //tree= tree->CopyTree("2.012-6*0.1775  < excitation && excitation .q<2.012+6*0.1775");  // 第一激发态 sigma 激发态：2.1 
    //tree= tree->CopyTree("2.012-0.1775  < excitation && excitation <2.012+0.1775");  // 第一激发态 sigma 激发态：2.1 
    
    //tree= tree->CopyTree("2.6-3*0.1775  < excitation && excitation <2.6+3*0.1775");  // sigma 激发态：2.75
    tree= tree->CopyTree("4.532-3*0.1775 < excitation && excitation < 4.532+3*0.1775");  //  sigma ：0.1775
    //tree= tree->CopyTree("5.75-3*0.1775 < excitation && excitation < 5.75 +3*0.1775");  //  sigma ：0.1775
   
    //tree= tree->CopyTree("6.6-3*0.1775 < excitation && excitation < 6.6 +3*0.1775");  //  sigma ：0.1775
    
    // tree= tree->CopyTree("-3*0.1775 < excitation && excitation < 3*0.1775");  //  sigma ：0.1775
    // tree= tree->CopyTree("0 < theta_cm && theta_cm < 180");  // degree 卡质子的质心系角度  
    this->loopTree(tree); 
    can->Print("12BeAnalysis.pdf");
    //hBepolar=this->minsPileUp(hBepolar,2.012,"AnaThetaCM21","AnaThetaCMsub27");// 文件名：分析的激发态  减去pileup的激发态 
    // hBepolar=this->minsPileUp(hBepolar,5.75,"AnaThetaCM57","AnaThetaCM66");// 文件名：分析的激发态  减去pileup的激发态
    //*****************************************自由CUT分析***********************************//
    
    
    delete can; can=nullptr;
    
    //********************************刻度效率、截面、DWBA、标准化********************************//

    // step 4：Be-12 角度分布
    TCanvas* c5=new TCanvas("c5", "c5", 1000, 600);
    c5->Divide(2,1);
    c5->cd(1);  
    hex->Fit(gaus,"MW","",-1,6);
    //hex->Draw();
    c5->cd(2);  
    hEZ->Draw("colz");// 激发态与位置关系
    //hBepolar->Draw();
    c5->Draw(); c5->Print("12BeAnalysis.pdf");

    
    //********************************* 新增部分：分角度绘制激发能谱***********************************************
    TCanvas* cAngleBins = new TCanvas("cAngleBins", "Excitation per Angle Bin", 1200, 900);
    cAngleBins->Print("ExcitationPerAngleBin.pdf("); // 开始PDF文件
    hBepolar->Reset(); // 重置hBepolar以用于新的角度分布分析

    // 将0-180度分成36个bin，每个5度
    const int nAngleBins = 18;
    const double binWidth = 180.0 / nAngleBins;
    
    TH1D* hexPerBin[nAngleBins];
    TString histName, histTitle;

    // bgCoeffMap 存储角度区间对应的多个系数, 成员变量
    ifstream bgFile("background_coefficients.txt");
    double angle, coeff, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3;// 读取系数
    while (bgFile >> angle >> coeff>> A1 >> mu1 >> sigma1 >> A2 >> mu2 >> sigma2 >> A3 >> mu3 >> sigma3) {
        //bgCoeffMap[angle] = coeff; // 存储到 map
        bgCoeffMap[angle].push_back(coeff);
        bgCoeffMap[angle].push_back(A1);
        bgCoeffMap[angle].push_back(mu1);
        bgCoeffMap[angle].push_back(sigma1);                
        bgCoeffMap[angle].push_back(A2);
        bgCoeffMap[angle].push_back(mu2);
        bgCoeffMap[angle].push_back(sigma2);
        bgCoeffMap[angle].push_back(A3);
        bgCoeffMap[angle].push_back(mu3);
        bgCoeffMap[angle].push_back(sigma3);
    }
    bgFile.close();

    for (int i = 0; i < nAngleBins; i++) {
        // 定义当前角度范围
        double minAngle = i * binWidth;
        double maxAngle = (i + 1) * binWidth;
        double bin_center = (minAngle + maxAngle) / 2 ;
        double bin_center_cm = 180 - bin_center ;
        
        // 创建只包含本底的函数
        TF1 *fBgOnly = new TF1("fBgOnly", ParticleDynamics::bgOnly, 4, 12, 1);
        double bgCoeff = bgCoeffMap.count(bin_center_cm) ? bgCoeffMap[bin_center_cm][0] : 0.0;// 从 map 中获取容器第一个值[0]
        fBgOnly->SetParameter(0,bgCoeff);
        fBgOnly->SetLineColor(kBlue);
        fBgOnly->SetLineStyle(2);
        fBgOnly->SetLineWidth(2);

        // 创建直方图
        histName = Form("hex_bin_%d", i);
        histTitle = Form("Excitation Energy (%.1f#circ - %.1f#circ);Excitation Energy (MeV);Counts", 
                         (180-minAngle), (180-maxAngle));
        hexPerBin[i] = new TH1D(histName, histTitle, 91, -2, 12); // 91 bins from -2 to 10 MeV
        
        // 应用角度筛选条件
        TCut angleCut = Form("theta_cm >= %.2f && theta_cm < %.2f", minAngle, maxAngle);
        TTree* treeAngleBin = tree->CopyTree(angleCut);
        
        // 填充直方图
        treeAngleBin->Draw(Form("excitation>>%s", histName.Data()));
        
        // 绘制直方图
        cAngleBins->cd();
        hexPerBin[i]->SetLineColor(kBlue);
        hexPerBin[i]->SetLineWidth(2);
        hexPerBin[i]->Draw();
        fBgOnly->Draw("same");

        TF1 *fCombined = new TF1("fCombined", ParticleDynamics::combinedFit, 4, 12, 10);
        this->SetBgPar(fCombined); // 参数命名与设定拟合初值

        if((180-maxAngle)==20)
            fCombined->SetParLimits(0,  0, 0.00005);  // 对mu值限制
        else 
            fCombined->SetParLimits(0,  0, 0.00008);  // 对mu值限制

        // 第一个峰 (MeV)
        fCombined->SetParameter(1, 50);  // P1_A
        fCombined->SetParameter(2, 4.412);  // P1_mu 固定在0
        fCombined->SetParLimits(2,  4.412-0.2, 4.412+0.2);  // 对mu值限制

        // 第二个峰 (MeV)
        fCombined->SetParameter(4, 30);  // P2_A
        fCombined->SetParameter(5, 5.724); // P2_mu 固定在1.78
        fCombined->SetParLimits(5,5.724-0.2, 5.724+0.2);  //对mu值限制

        // // 第三个峰 (MeV)
        fCombined->SetParameter(7, 30);  // P3_A
        fCombined->SetParameter(8, 6.6); // P3_mu 固定在2.65
        fCombined->SetParLimits(8,6.6-0.2, 6.6+0.2);  // 对mu值限制

     
        // //执行拟合
        hexPerBin[i]->Fit(fCombined, "R+", "same", 4, 8);
        // //保存本底系数（par[0]）到文件，格式：角度中心值 本底系数
        // ofstream bgFile("background_coefficients.txt", ios::app); // 追加模式  一行本底+一行峰系数
        // bgFile << bin_center_cm << "\t" << fCombined->GetParameter(0) <<"\t\t";// 保存本底系数
        // bgFile << fCombined->GetParameter(1) <<"\t"<<fCombined->GetParameter(2) << "\t"<<fCombined->GetParameter(3) <<"\t";// 保存第一个峰系数
        // bgFile << fCombined->GetParameter(4) <<"\t"<<fCombined->GetParameter(5) << "\t"<<fCombined->GetParameter(6) <<"\t";// 保存第二个峰系数
        // bgFile << fCombined->GetParameter(7) <<"\t"<<fCombined->GetParameter(8) << "\t"<<fCombined->GetParameter(9) <<endl;// 保存第三个峰系数
        // bgFile.close();

        // 计算峰区域总积分
        Int_t bin_low = hexPerBin[i]->GetXaxis()->FindBin(peak_low);
        Int_t bin_high = hexPerBin[i]->GetXaxis()->FindBin(peak_high);
        double total_integral = hexPerBin[i]->Integral(bin_low, bin_high);

        // 初始化净信号和本底加和
        double net_signal = 0.0;
        double total_bg = 0.0;
        double minValue = 0.0;
        // 逐个 bin 处理
        for (Int_t b = bin_low; b <= bin_high; b++) {
             double data_count = hexPerBin[i]->GetBinContent(b);
              
             double bin_center = hexPerBin[i]->GetXaxis()->GetBinCenter(b);
                         
            // 获取背景函数在当前 bin 中心的值
             double bg_value = fBgOnly->Eval(bin_center);
                                      
             total_bg += bg_value;
                                                 
            // 比较实际计数与背景值、堆叠值，计算净信号
             if (data_count > bg_value){
                // 如果实际计数大于背景值，认为有信号
                minValue=data_count - bg_value; 

                //this->minsPileUp(hexPerBin[i], bin_center_cm, minValue, b, 2);// 减去第二个峰堆叠的影响, start from 1
                //this->minsPileUp(hexPerBin[i], bin_center_cm, minValue, b, 3);// 减去第三个峰堆叠的影响

                net_signal += minValue;
            }else {
                // 如果实际计数小于等于背景值，全部认为是本底
                // 不需要额外操作，因为信号为0
            }
     
        }// mins bin counts loop!
        //cout<<"total_intefral = "<<total_integral<<", bg_integral = "<<total_bg<<", net_signal = "<<net_signal<<endl;

        // 填充直方图
        hBepolar->Fill(bin_center_cm, total_integral);// 指定bin的填充值
        hBepolar_bgsub->Fill(bin_center_cm, net_signal);

        TLatex* result = new TLatex(0.3, 0.75, Form("Total: %.1f\n  Real_cts: %.1f", total_integral, net_signal));
        result->SetNDC();
        result->SetTextSize(0.05);
        result->Draw("same");
        
        // 添加角度范围标签
        TLatex* latex = new TLatex(0.3, 0.85, Form("%.1f#circ - %.1f#circ", (180-minAngle), (180-maxAngle)));
        latex->SetNDC();
        latex->SetTextSize(0.06);
        latex->Draw("same");
        
        // 保存当前页到PDF
        hexPerBin[i]->GetYaxis()->SetRangeUser(0., 10.);
        cAngleBins->Modified();
        cAngleBins->Update();
        cAngleBins->Print("ExcitationPerAngleBin.pdf");
        
        // 清理
        delete treeAngleBin;
        delete latex;
    }
    cAngleBins->Print("ExcitationPerAngleBin.pdf)"); // 结束PDF文件
    delete cAngleBins;



    // 获取粒子重建效率
    //this->GetEfficiency(hBepolar,"./Efficiency/Efficiency_test.txt");

    // step 5：得到微分截面
    TCanvas* cSigma=new TCanvas("cSigma", "dSigma", 1000, 600);
    TGraphErrors* grCrossSection=new TGraphErrors;
    
    // 该函数暂时已弃用
    //this->minsGaus(hBepolar,2.6-3*0.1775,2.6);//扣除gaus拟合的本底  
    //this->minsGaus(hBepolar,2.1-3*0.1775,2.1+3*0.1775);//扣除gaus拟合的本底

    // 设置重建效率、计算截面
    this->GetCrossSection(hBepolar_bgsub,grCrossSection);
    //->GetCrossSection(hBepolar,grCrossSection);
    //this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency00.txt");
    //this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency21.txt");
    //this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency27.txt");
    this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency45.txt");
    //this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency57.txt");
    //this->SetEfficiency(grCrossSection,"./Efficiency/Efficiency66.txt");

    // DWBA计算截面结果
    TString path="./DWBA/"; 
    //this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC.txt");
    //this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC21.txt");//2.1 Mev
    //this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC27.txt");//2.7 Mev
    this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC45.txt");//2.7 Mev
    //this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC57.txt");//2.7 Mev
    //this->ReadDWBA(grDWBA,path+"DWBA.Xsec.VVC66.txt");//2.7 Mev


    // this->ReadDWBA(grDWBA2,"DWBA.XsecVV572.txt");//2.7 Mev 
    // this->ReadDWBA(grDWBA3,"DWBA.XsecVV573.txt");//2.7 Mev
    // grDWBA2->SetMarkerColor(3);
    // grDWBA3->SetMarkerColor(12);

    // 标准化基态微分截面--
    // 推荐缩放：Pulse Height Normalization Factor :0.447514      Deviation:-0.0586111
    grCrossSection->Scale(1./0.370);//激发态标准化 
    grCrossSection->Scale(1./0.5);  //vertex_z 效率标准化
    double scaleFactor=this->Normalize(grDWBA,grCrossSection,28,70);// 设定放缩因子IC2拟合区间 == 22 38  、28-70
    if(scaleFactor<0.3){// IC1== 0.400316\0.368813   IC2==0.447514  IC1 18=0.370 (250718)
        scaleFactor= 0.400316 ;
    }  scaleFactor=  0.290;
   // grCrossSection->Scale(1./scaleFactor);

    // 输出截面数据至文本文件
    grCrossSection=this->CutCrossSectionData(grCrossSection,28.,140.);//截取数据28-100
    this->WriteCrossSection(grCrossSection,"./CrossSection/CrossSectionbg.extestMeV.txt");
    
    grCrossSection->GetXaxis()->CenterTitle(); 
	grCrossSection->GetYaxis()->CenterTitle();
    grCrossSection->GetXaxis()->SetTitleSize(0.05);
    grCrossSection->GetXaxis()->SetTitleOffset(0.90);
    grCrossSection->GetYaxis()->SetTitleSize(0.05);
    grCrossSection->GetYaxis()->SetTitleOffset(0.9);
    grCrossSection->GetYaxis()->SetRangeUser(0.01, 1e2);
    grCrossSection->Draw("AP Same"); //AP
    grDWBA->Draw("L same");
    //grDWBA2->Draw("P same");
    //grDWBA3->Draw("P same");

    // 创建图例对象
    TLegend* legend = new TLegend(0.77, 0.77, 0.87, 0.89);// 图例的位置 (x1, y1, x2, y2)
    legend->AddEntry(grCrossSection, "Data", "p"); // 添加图例条目 (图形对象, 标签, 图例标志)  
    legend->AddEntry(grDWBA, "DWBA", "L");
    //legend->AddEntry(grDWBA2, "DWBA 3-", "P");
    //legend->AddEntry(grDWBA3, "DWBA 4+", "P");
    legend->SetBorderSize(0); 
    legend->Draw();// 绘制图例
   
    // cSigma->SetGridx();
    // cSigma->SetGridy();  
    // cSigma->SetTicks(1,1);
    cSigma->SetLogy();
    cSigma->SetFrameLineWidth(2);
    cSigma->SetFrameBorderMode(0);
    cSigma->SetFrameLineWidth(2);
    cSigma->SetFrameBorderMode(0);
    cSigma->Update();
    cSigma->Print("12BeAnalysis.pdf)");


    cSigma->Update();

}// ParticleDynamics::dataAnalysis end here! 

void ParticleDynamics::SetBgPar(TF1* fbg){// 对函数的参数设定拟合初值

    fbg->SetParName(0, "BgAmp");      // 本底权重
    fbg->SetParName(1, "A_{1}");      // 第一个高斯幅度
    fbg->SetParName(2, "#mu_{1}");    // 第一个高斯幅度
    fbg->SetParName(3, "#sigma_{1}"); // 第一个高斯标准差
    fbg->SetParName(4, "A_{2}");      // 第二个高斯幅度
    fbg->SetParName(5, "#mu_{2}");    // 第二个高斯均值
    fbg->SetParName(6, "#sigma_{2}"); // 第二个高斯标准差
    fbg->SetParName(7, "A_{3}");      // 第三个高斯幅度
    fbg->SetParName(8, "#mu_{3}");    // 第三个高斯均值
    fbg->SetParName(9, "#sigma_{3}"); // 第三个高斯标准差

    fbg->SetParLimits(1, 0., 1e5);   // par限制  使峰幅度不为负值
    fbg->SetParLimits(4, 0., 1e5);   // par限制
    fbg->SetParLimits(7, 0., 1e5);   // par限制
    fbg->SetParLimits(3, 0.1, 0.2);  // par限制
    fbg->SetParLimits(6, 0.1, 0.2);  // par限制
    fbg->SetParLimits(9, 0.1, 0.2);  // par限制
   

}

double ParticleDynamics::bgOnly(double *x, double *par) {
    double xx = x[0];
    double bg = 0.0;
    if (xx > 4 && xx < 12.0){
        bg = par[0] * (3.21377 * pow(xx,5) - 139.852 * pow(xx,4) + 2353.49 * pow(xx,3) -19890.5 * pow(xx,2) + 86793.3 * xx - 140050);
   }
    //cout<<"xx = "<< xx <<", bg = "<<bg<<endl;
    return bg;
}// ParticleDynamics:bgOnly end here! 

double ParticleDynamics::combinedFit(double *x, double *par) {

   double xx = x[0];
                                
   double bg = 0.0;
   if (xx > 4 && xx < 12.0){
        bg = par[0] * (3.21377 * pow(xx,5) - 139.852 * pow(xx,4) + 2353.49 * pow(xx,3) -19890.5 * pow(xx,2) + 86793.3 * xx - 140050);
   }
   double peak1 = 0.0;
   double peak2 = 0.0;
   double peak3 = 0.0;
   double peak4 = 0.0;
   double peak5 = 0.0;
   // 四个高斯峰
   if(xx >= -1 && xx <= 8.0){
        peak1 = par[1] * TMath::Exp(-0.5*TMath::Power((xx - par[2])/par[3], 2));
        peak2 = par[4] * TMath::Exp(-0.5*TMath::Power((xx - par[5])/par[6], 2));
        peak3 = par[7] * TMath::Exp(-0.5*TMath::Power((xx - par[8])/par[9], 2));
        // peak4 = par[10] * TMath::Exp(-0.5*TMath::Power((xx - par[11])/par[12], 2));
        // peak5 = par[13] * TMath::Exp(-0.5*TMath::Power((xx - par[14])/par[15], 2));
   }
   //                                                                 
   return bg + peak1 + peak2 + peak3;// + peak4 + peak5;
}// ParticleDynamics:combinedFit end here! 




void ParticleDynamics::loopTree(TTree* tree, double fitLeft, double fitRight, bool reset){

    TCanvas* can= (TCanvas*) gROOT->FindObject("can");

    TH1D* hKe= (TH1D*) gROOT->FindObject("Ke");
    TH1D* hdistance= (TH1D*) gROOT->FindObject("hdistance");
    TH1D* hex= (TH1D*) gROOT->FindObject("hex");
    TH2D* hKL= (TH2D*) gROOT->FindObject("hKL");
    TH2D* hKcm= (TH2D*) gROOT->FindObject("hKcm");
    TH2D* hEZ= (TH2D*) gROOT->FindObject("hEZ");
    TH1D* hBepolar= (TH1D*) gROOT->FindObject("hBepolar");
    
    if(reset){
        can->Clear();
        hKe->Reset();;
        hdistance->Reset();;
        hex->Reset();
        hKcm->Reset();
        hBepolar->Reset();   
    }
     

    Long64_t nentries = tree->GetEntriesFast();//cout<<"after cut  nentries:"<<nentries<<endl;
    
    double distance , projectile_ke, ex, polarCM;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
       tree->GetEntry(jentry); 

   
       distance=sqrt(vertex_x*vertex_x+vertex_y*vertex_y+vertex_z*vertex_z)*100.;// 单位 cm
       projectile_ke=GetEnergyLoss(distance);
        

       hKe->Fill(projectile_ke);
       hdistance->Fill(distance);

       ex=GetExciteEnergy(projectile_ke, ke, polar);
       hex->Fill(ex);
       hEZ->Fill(ex,vertex_z);

       polarCM=AngleTrans(projectile_ke, ke, polar);

       hKL->Fill(polar,ke);
       hKcm->Fill(polarCM,ke);

       hBepolar->Fill(180.- polarCM);

    }

    can->Divide(2,2);
    can->cd(1);  hKe->Draw();
    can->cd(2);  hdistance->Draw();
    can->cd(3);  hKcm->Draw();   
    can->cd(4);  hex->Draw(); // 激发态

    // ----->获取拟合参数
    //hex ->Fit(gaus,"QN","",fitLeft,fitRight);// 激发态
    //double Sigma=gaus->GetParameter(2);//cout<<"test取值范围par[1]:"<<Sigma<<endl;

    hdistance->GetXaxis()->SetTitle("Range of Track [cm]"); 
	hdistance->GetYaxis()->SetTitle("Counts");
    hdistance->GetXaxis()->CenterTitle(); 
	hdistance->GetYaxis()->CenterTitle();
    
    hBepolar->SetTitle("Be-12");
    hBepolar->GetXaxis()->SetTitle("#theta / degree"); 
	hBepolar->GetYaxis()->SetTitle("Counts");
    hBepolar->GetXaxis()->CenterTitle(); 
	hBepolar->GetYaxis()->CenterTitle();

    hKe->GetXaxis()->SetTitle("#theta / degree"); 
	hKe->GetYaxis()->SetTitle("Kinetic Energy / MeV");
    hKe->GetXaxis()->CenterTitle(); 
	hKe->GetYaxis()->CenterTitle();

    hKcm->GetXaxis()->SetTitle("#theta / degree"); 
	hKcm->GetYaxis()->SetTitle("Kinetic Energy / MeV");
    hKcm->GetXaxis()->CenterTitle(); 
	hKcm->GetYaxis()->CenterTitle();
    
    hex->SetTitle("Be-12 Spectrum");
    hex->GetXaxis()->SetTitle("Excitation Energy [MeV]"); 
	hex->GetYaxis()->SetTitle("Counts");
    hex->GetXaxis()->CenterTitle(); 
	hex->GetYaxis()->CenterTitle();

}// ParticleDynamics::loopTree end here! 


void ParticleDynamics::angleResolution(){// test fun

    // 运动学参考线
      TGraph *gr = new TGraph;      // GetRefLine
   
      TGraph *grex1 = new TGraph;   // GetRefLine
      TGraph *grex2 = new TGraph;   // GetRefLine
      TGraph *grex3 = new TGraph;   // GetRefLine
      TGraph *grex4 = new TGraph;   // GetRefLine
      TGraph *grex5 = new TGraph;   // GetRefLine
      TGraph *grCM = new TGraph;    // GetRefLine
      TGraph *grLise = new TGraph;    // GetRefLine
      TGraph *grCMex1 = new TGraph; // GetRefLine
      TGraph *grCMex2 = new TGraph; // GetRefLine
      TGraph *grCMex3 = new TGraph; // GetRefLine
      TGraph *grCMex4 = new TGraph; // GetRefLine
      TGraph *grCMex5 = new TGraph; // GetRefLine
      this->GetRefLineRe(gr,grCM,0.); 
      this->GetRefLineRe(grex1,grCMex1,0.961);
      this->GetRefLineRe(grex2,grCMex2,1.19);
      this->GetRefLineRe(grex3,grCMex3,1.8);
      this->GetRefLineRe(grex4,grCMex4,2.439);
      this->GetRefLineRe(grex5,grCMex5,3.132);
      
      

      TGraph *dgr = this->computeDiff(gr);  // 差分：GetRefLine

      gr->GetXaxis()->SetTitle("#theta [degree]"); 
      gr->GetYaxis()->SetTitle("Kinetic Energy [MeV]");
      gr->GetXaxis()->CenterTitle(); 
      gr->GetYaxis()->CenterTitle();
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetXaxis()->SetTitleOffset(0.90);
      gr->GetYaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetTitleOffset(0.9);

      TH2D *hTheta=new TH2D("hTheta","",1000,40,90,1000,40,90);
      TH1D *hde=new TH1D("hde","",1000,-5,5);
      TH2D *hde2=new TH2D("hde2","",1000,10,90,1000,-5,5);
      TH1D *h1_theta=new TH1D("h1_dtheta","",1000,-10,10);
      TH1D *hex=new TH1D("hex","",150,-2.,7);// excitation : MeV
     
 
      this->ReadDWBA(grDWBA,"DWBA.Xsec12N.txt");//gs MeV
      this->ReadDWBA(grDWBA2,"DWBA.Xsec12Nex1.txt");//0.961 MeV
      this->ReadDWBA(grDWBA3,"DWBA.Xsec12Nex3.txt");//1.8 MeV
      //   this->ReadDWBA(grDWBA4,"DWBA.Xsec12Nex4.txt");//2.439 MeV



      // 已知均匀分布的方差为 var_uniform
      double var_uniform = 0.268; // 标准差
      double half_width = TMath::Sqrt(3)* var_uniform;
      double power;// 权重

      double ex;
       hex->Fill(ex);

       //  gs
      for(int j=55;j<85;j++){
            for(int i=0;i<1000;i++){// 抽样次数
                //double de=gr->Eval(r->Gaus(j, 0.4128))-gr->Eval(j);// TPC高斯分布
                double Jteheta=r->Uniform(j, j + 1);// 角度均匀分布
                double thetaL=TMath::Pi()*r->Uniform(Jteheta - half_width, Jteheta + half_width)/180;
                double ex=GetExciteEnergy(proj_energy_start, gr->Eval(Jteheta) ,thetaL); // 计算激发能
                power=this->GetDWBApdf(grDWBA,gr,Jteheta);
                hde->Fill(ex,power);
                hde2->Fill(Jteheta,ex);
        }   
      }
    // // gs 考虑能量分辨
    //   for(int j=55;j<85;j++){
    //         for(int i=0;i<1000;i++){// 抽样次数
    //             //double de=gr->Eval(r->Gaus(j, 0.4128))-gr->Eval(j);// TPC高斯分布
    //             double Jteheta=r->Uniform(j, j + 1);// 角度均匀分布
    //             double de=gr->Eval(r->Uniform(Jteheta - half_width, Jteheta + half_width))-gr->Eval(Jteheta); // 硅均匀分布
    //             power=this->GetDWBApdf(grDWBA,gr,Jteheta);
    //             hde->Fill(de,power);
    //             hde2->Fill(Jteheta,de);
    //     }   
    //   }
      
     

      // 0.961 grex1
    //   for(int j=55;j<79;j++){
    //         for(int i=0;i<1000;i++){// 抽样次数
    //             //double de=gr->Eval(r->Gaus(j, 0.4128))-gr->Eval(j);// TPC高斯分布
    //             double Jteheta=r->Uniform(j, j + 1);// 角度均匀分布
    //             double thetaL=TMath::Pi()*r->Uniform(Jteheta - half_width, Jteheta + half_width)/180; // 硅均匀分布
    //             double ex=GetExciteEnergy(proj_energy_start, grex1->Eval(Jteheta) ,thetaL); // 计算激发能
    //             //power=this->GetDWBApdf(grDWBA2,grex1,Jteheta);
    //             power=this->GetDWBApdf(grDWBA2,grex1,Jteheta,55,79);
    //             hde->Fill(ex,power);
    //             hde2->Fill(Jteheta,ex);
    //     }   
    //   }

    //  // 1.8 grex3
    //   for(int j=55;j<70;j++){
    //         for(int i=0;i<1000;i++){// 抽样次数
    //             //double de=gr->Eval(r->Gaus(j, 0.4128))-gr->Eval(j);// TPC高斯分布
    //             double Jteheta=r->Uniform(j, j + 1);// 角度均匀分布
    //             double thetaL=TMath::Pi()*r->Uniform(Jteheta - half_width, Jteheta + half_width)/180; // 硅均匀分布
    //             double ex=GetExciteEnergy(proj_energy_start, grex3->Eval(Jteheta) ,thetaL); // 计算激发能
    //             power=this->GetDWBApdf(grDWBA3,grex3,Jteheta,55,74);
    //             hde->Fill(ex,power);
    //             hde2->Fill(Jteheta,ex);
    //     }   
    //   }

      // 创建画布和3D视图
      TCanvas* c = new TCanvas("c", "Si Detectors", 1200, 800);
      c->Divide(2,2);
      c->cd(1);
     // gr->Draw(); // "L" 表示画线
    //  grex1->Draw("L"); // "L" 表示画线
      //grex2->Draw("L"); // "L" 表示画线
    //  grex3->Draw("L"); // "L" 表示画线
    //  grex4->Draw("L"); // "L" 表示画线
      //grLise->Draw("L"); // "L" 表示画线
     // grCM->Draw("L"); // "L" 表示画线
    //   grCMex1->Draw("L"); // "L" 表示画线
    //   grCMex3->Draw("L"); // "L" 表示画线
    //   grCMex4->Draw("L"); // "L" 表示画线

      c->cd(2);
      dgr->Draw(); // "L" 表示画线
      hde->GetXaxis()->SetTitle("Excitation Energy [MeV]"); 
      hde->GetXaxis()->CenterTitle();
      hde->Draw();
      //hex->Draw();
      c->cd(3);
      hTheta->GetXaxis()->SetTitle("#theta L"); 
      hTheta->GetYaxis()->SetTitle("Cal #theta L");
      TF1 *f1 = new TF1("f1", "[0]+[1]*x", 0, 180);
      //hTheta->Fit(f1, "R");
      hTheta->Draw("colz"); 
      c->cd(4);
 
      hde2->GetXaxis()->SetTitle("#theta degree "); 
      hde2->GetYaxis()->SetTitle("Deviation [MeV]");
      hde2->GetYaxis()->CenterTitle(); 
      hde2->GetXaxis()->CenterTitle();
      hde2->GetXaxis()->SetTitleSize(0.05);
      hde2->GetXaxis()->SetTitleOffset(0.90);
      hde2->GetYaxis()->SetTitleSize(0.05);
      hde2->GetYaxis()->SetTitleOffset(0.9);
      hde2->Draw("colz"); 
    //   h1_theta->GetXaxis()->SetTitle("#theta L deviation"); 
    //   h1_theta->GetXaxis()->CenterTitle();
    //   h1_theta->Draw(); 
    //   h1_theta->Fit("gaus");
      
 

}// ParticleDynamics::test end here!

double ParticleDynamics::GetDWBApdf(TGraph* grDwBA, TGraph* gr,double powerTheta, double thetaL, double thetaR) {

    // 质心系转换 ke->MeV theta->radian  out: theta->degree
    // double AngleTrans(double projectile_ke,  double ejectile_ke, double theta);// 质心系转换

    double Emax=gr->Eval(thetaL);// 轻核实验室系高能小角度，
    double Emin=gr->Eval(thetaR);// 轻核实验室系低能大角度
    double E=gr->Eval(powerTheta);// 目标权重对应的能量
    // if(Emin<0.001) Emin=0.807912;
    // if(E<0.001) E=0.807912;
    double thetaCMmax =180- AngleTrans(proj_energy_start, Emax, thetaL*TMath::Pi()/180.0);// 重核质心系大角度
    double thetaCMmin =180- AngleTrans(proj_energy_start, Emin, thetaR*TMath::Pi()/180.0);// 重核质心系小角度
    double thetaCM =180- AngleTrans(proj_energy_start, E, powerTheta*TMath::Pi()/180.0);//   重核质心系：目标角度
             thetaCMmax=std::round(thetaCMmax);
             thetaCMmin=std::round(thetaCMmin);
             thetaCM=std::round(thetaCM);
    //cout<<"thetaCMmax="<<thetaCMmax<<" thetaCMmin="<<thetaCMmin<<" thetaCM="<<thetaCM<<endl;
    
    // 获取点的数量
    Int_t nPoints = grDwBA->GetN();
    double theta,Y;
    grDwBA->GetPoint(0, theta, Y); // 获取第 i 个点的 x 和 y 值
    double width=(180.-theta)/nPoints;// cout<<"width="<<width<<endl;

    //grDwBA->Integral(thetaCMmin, thetaCMmax);// 计算积分

    // 计算积分
    double integral = 0.0;
    for (Int_t i = thetaCMmin; i < thetaCMmax; i=i+width) {

        grDwBA->GetPoint(i, theta, Y); // 获取第 i 个点的 x 和 y 值: floor向下取整
        integral=integral + Y*width;
    }
    grDwBA->GetPoint(thetaCM, theta, Y); // 获取第 i 个点的 x 和 y 值
    
    //cout<<"powerT="<<Y*width/integral<<endl;
    return Y*width/integral;

}// ParticleDynamics::GetDWBApdf end here!

void ParticleDynamics::targetDeposit(){

    // 运动学参考线
    TGraph *gr = new TGraph;      // GetRefLine
    TGraph *grCMex1 = new TGraph;   // GetRefLine
    this->GetRefLineS(gr,grCMex1,2.4);

    TGraph *grtest = new TGraph;      // GetRefLine
    TGraph *grCMtest = new TGraph;    // GetRefLine
    this->GetBreakupRefLineS(grtest,grCMtest,10*10,10);

    this->ReadDWBA(grDWBA,"DWBA.Xsec12Npt.txt");//gs MeV
    
    TH1D *hex=new TH1D("hex","",100,-2,5.);// 出射能量分布 MeV
    TH1D *hde=new TH1D("hde","",100,-10,10);// 出射能量分布 MeV
    TH1D *hde2=new TH1D("hde2","",100,0,105);// 出射能量分布 MeV
    TH2D *hDepthEnergy = new TH2D("hDepthEnergy", "Reaction Depth vs Fragment Energy", 
                                  100, 0, 50, 500, 293,305); // 深度与能量关系 //100, 200, 400


    TH2D *hE_theta = new TH2D("hE_theta", "E vs theta  / MeV", 
                                    1000, 0, 90,1800, 0, 150); // 深度与能量关系 //100, 200, 400

    TH2D *hde_E = new TH2D("hde_E", "dE2 vs dE1  / MeV", 
                                    2800, 0, 280,2800, 0, 280); // 深度与能量关系 //100, 200, 400
    TH2D *hde_E2 = new TH2D("hde_E2", "dE3 vs dE2  / MeV", 
                                    2800, 0, 280,2800, 0, 280); // 深度与能量关系 //100, 200, 400 
    TH2D *hde_E3 = new TH2D("hde_E3", "dE4 vs dE3  / MeV", 
                                    2800, 0, 280,2800, 0, 280); // 深度与能量关系 //100, 200, 400   
    TH2D *hde_E4 = new TH2D("hde_E4", "dE5 vs dE4  / MeV", 
                                    2200, 0, 220,2200, 0, 220); // 深度与能量关系 //100, 200, 400                                                               
    TH2D *hde_E5 = new TH2D("hde_E5", "dE6 vs dE5  / MeV", 
                                    2200, 0, 220,2200, 0, 220); // 深度与能量关系 //100, 200, 400  

    double thickness_target=100*1e-4; // um -> cm

    double distance, projectile_ke;

    double  exit_energy,recoil_energy,exit_angle,recoil_angle;// degree, MeV

    double  fragment1_energy,fragment2_energy,fragment1_angle,fragment2_angle;// degree MeV
        
    double E_out;

    int SiCount1=0,SiCount2=0,SiCount3=0,SiCount4=0,SiCount5=0, CsICount=0; // Si detector count

    catima::Material T1Si({28,14,1});
    T1Si.density(2.321).thickness(1520.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material SiD1({28,14,1});
    SiD1.density(2.321).thickness(40.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material SiD2({28,14,1});
    SiD2.density(2.321).thickness(1000.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material SiD3({28,14,1});
    SiD3.density(2.321).thickness(1000.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material SiD4({28,14,1});
    SiD4.density(2.321).thickness(1500.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material SiD5({28,14,1});
    SiD5.density(2.321).thickness(1500.0*1e-4*2.321);//  g/cm3 and g/cm2

    catima::Material CsI({{133,55,1},{127,53,1}});
    CsI.density(4.51).thickness(40000.0*1e-4*4.51);//  g/cm3 and g/cm2


    for(int i=0;i<1000;i++){ // 距离抽样循环
       
        distance=r->Uniform(0, thickness_target);// um -> cm
        
        this->material.thickness(distance*this->density);   // thickness in g/cm2
        catima::Projectile p0(this->projectile_a,this->projectile_z); //
        p0.T = proj_energy; // setting energy to  MeV/u;
        this->projectile=p0;
        projectile_ke=catima::energy_out(this->projectile,this->material) * projectile_mass;
        
        for(Int_t j=1;j<1800;j++){ // 角度抽样循环
            
            // 第一步：非弹性散射
            double thetaCM1=r->Uniform(0, 180.)*TMath::Pi()/180.0;// 质心系均匀分布角度 TMath::Pi()/2.;//
            this->InelasticScattering(exit_energy,
                                    recoil_energy,
                                    exit_angle,
                                    recoil_angle, 
                                    projectile_ke,
                                    15, //exciteEnergy
                                    thetaCM1);
            // 调整能量（减去静止质量）
            //recoil_angle=recoil_angle*180.0/TMath::Pi();            // 轻核角度/反冲核角度
            //exit_angle=exit_angle*180.0/TMath::Pi();               // 剩余核出射粒子角度 
            recoil_energy=recoil_energy-ejectile_mass*u_MeV;      // 轻核能量
            exit_energy=exit_energy-residual_mass*u_MeV;         // 剩余核出射粒子能量
            //cout<<"i:"<<i<<"\t thetaCM:"<<thetaCM1<<"\t exit_energy:"<<exit_energy<<"\t recoil_energy:"<<recoil_energy-target_mass*u_MeV<<endl;
            
            
            // 第二步：破裂反应
            double thetaCM2=r->Uniform(0, 180.)*TMath::Pi()/180.0;// 质心系均匀分布角度  TMath::Pi()/2.; //
            this->BreakupReaction(this->residual_mass,
                                this->fragment1_mass,
                                this->fragment2_mass,
                                exit_energy,//Ek
                                15,//exciteEnergy
                                thetaCM2,
                                fragment1_energy,
                                fragment2_energy,
                                fragment1_angle,
                                fragment2_angle);

            // 第三步：计算碎片能量
            fragment1_energy=fragment1_energy-this->fragment1_mass*u_MeV;
            fragment2_energy=fragment2_energy-this->fragment2_mass*u_MeV;


            // 第四步：计算碎片在材料中的能量损失                       
            // 计算碎片1的实验室角度
            TVector3 frag1_vec(0, 0, 1);    // 视为第二个束流坐标系，初始方向沿z轴 实际上是角度相加
            frag1_vec.RotateY(fragment1_angle);         //cout<<"frag1_angle="<<fragment1_angle*180.0/TMath::Pi()<<endl;
            frag1_vec.RotateY(exit_angle);              //cout<<"exit_angle="<<exit_angle*180.0/TMath::Pi()<<endl;
            double frag1_lab_angle = frag1_vec.Theta(); //cout<<"frag1_lab_angle="<<frag1_lab_angle*180.0/TMath::Pi()<<endl;


            
            // 计算碎片2的实验室角度
            TVector3 frag2_vec(0, 0, 1);
            frag2_vec.RotateY(fragment2_angle);
            frag2_vec.RotateY(exit_angle);
            double frag2_lab_angle = frag2_vec.Theta();//cout<<"frag2_lab_angle="<<frag2_lab_angle*180.0/TMath::Pi()<<endl;

            
            // 计算实际穿过的材料厚度（考虑角度）
            double effective_thickness0 = (thickness_target - distance) / TMath::Abs(TMath::Cos(recoil_angle));
            double effective_thickness1 = (thickness_target - distance) / TMath::Abs(TMath::Cos(frag1_lab_angle));
            double effective_thickness2 = (thickness_target - distance) / TMath::Abs(TMath::Cos(frag2_lab_angle));


            // 设置材料属性、计算出射靶后的碎片能量

            // 计算出射轻核  在靶中能量损失
            this->material.thickness(effective_thickness0*this->density);   // thickness in g/cm2
            catima::Projectile particle0(this->ejectile_a,this->ejectile_z); //
            particle0.T = recoil_energy/this->ejectile_mass; // setting energy to  MeV/u;
            double E0_out=catima::energy_out(particle0,this->material) * ejectile_mass; // MeV
            //cout<<"E0_out="<<E0_out<<"\t recoil_angle:"<<recoil_angle*180.0/TMath::Pi()<<endl;

            // 计算出射轻核在硅中的能量损失
            particle0.T = E0_out/this->ejectile_mass; // setting energy to  MeV/u;
            double ET1Si=catima::energy_out(particle0,T1Si) * ejectile_mass; // 硅出射能量MeV    
            double dET1Si=E0_out-ET1Si; // 硅能量沉积MeV  
            hE_theta->Fill(recoil_angle*180.0/TMath::Pi(),E0_out);
            
            double ex=GetExciteEnergy(proj_energy_start, E0_out,recoil_angle); // 计算激发能
            double power=this->GetDWBApdf(grDWBA,gr,recoil_angle);
            //if((recoil_angle*180.0/TMath::Pi())>(70-12)&&(recoil_angle*180.0/TMath::Pi())<(70+12)&&E0_out>0.001)
            hex->Fill(ex,power);


            // 计算碎片1  C-11在靶中能量损失
            this->material.thickness(effective_thickness1*this->density);   // thickness in g/cm2
            catima::Projectile particle1(this->fragment1_a,this->fragment1_z); //
            particle1.T = fragment1_energy/this->fragment1_mass; // setting energy to  MeV/u;
            double E1_out=catima::energy_out(particle1,this->material) * fragment1_mass; // MeV

            // 计算出射碎片1在硅中的能量损失
            particle1.T = E1_out/this->fragment1_mass; // setting energy to  MeV/u;
            double ESiD1=catima::energy_out(particle1,SiD1) * fragment1_mass; // 第一块硅出射能量MeV    
            double dESiD1=E1_out-ESiD1; // 第一块硅能量沉积MeV 
            if(ESiD1==0) SiCount1++;
            
            particle1.T = ESiD1/this->fragment1_mass; // setting energy to  MeV/u;
            double ESiD2=catima::energy_out(particle1,SiD2) * fragment1_mass; // 第二块硅出射能量MeV    
            double dESiD2=ESiD1-ESiD2; // 第二块硅能量沉积MeV 
            if(ESiD2==0&&ESiD1!=0) SiCount2++;
            
            
            particle1.T = ESiD2/this->fragment1_mass; // setting energy to  MeV/u;
            double ESiD3=catima::energy_out(particle1,SiD3) * fragment1_mass; // 第三块硅出射能量MeV  
            double dESiD3=ESiD2-ESiD3; // 第三块硅能量沉积MeV 
            if(ESiD3==0&&ESiD2!=0) SiCount3++;
            
            particle1.T = ESiD3/this->fragment1_mass; // setting energy to  MeV/u; 
            double ESiD4=catima::energy_out(particle1,SiD4) * fragment1_mass; // 第四块硅出射能量MeV  
            double dESiD4=ESiD3-ESiD4; // 第四块硅能量沉积MeV
            if(ESiD4==0&&ESiD3!=0) SiCount4++;

            particle1.T = ESiD4/this->fragment1_mass; // setting energy to  MeV/u; 
            double ESiD5=catima::energy_out(particle1,SiD5) * fragment1_mass; // 第五块硅出射能量MeV  
            double dESiD5=ESiD4-ESiD5; // 第五块硅能量沉积MeV
            if(ESiD5==0&&ESiD4!=0) SiCount5++;

            
            particle1.T = ESiD5/this->fragment1_mass; // setting energy to  MeV/u; 
            double ECsI=catima::energy_out(particle1,CsI) * fragment1_mass; // CsI出射能量MeV  
            double dECsI=ESiD5-ECsI; // CsI能量沉积MeV
            if(ECsI==0) CsICount++;
            // cout<<"dECsI="<<dECsI<<" dESiD5="<<dESiD5<<" dESiD4="<<dESiD4<<" dESiD3="<<dESiD3
            // <<" dESiD2="<<dESiD2<<" dESiD1="<<dESiD1<<" E1_out:"<<E1_out<<endl;
 


            // 计算碎片2  p在靶中能量损失
            this->material.thickness(effective_thickness2*this->density);   // thickness in g/cm2
            catima::Projectile particle2(this->fragment2_a,this->fragment2_z); //
            particle2.T = fragment2_energy/this->fragment2_mass; // setting energy to  MeV/u;
            double E2_out=catima::energy_out(particle2,this->material) * fragment2_mass;// MeV

            // 计算出射碎片2在硅中的能量损失
            particle2.T = E2_out/this->fragment2_mass; // setting energy to  MeV/u;
            double ESiD1_=catima::energy_out(particle2,SiD1) * fragment2_mass; // 第一块硅出射能量MeV    
            double dESiD1_=E2_out-ESiD1_; // 第一块硅能量沉积MeV 
            if(ESiD1_==0) SiCount1++;
            
            particle2.T = ESiD1_/this->fragment2_mass; // setting energy to  MeV/u;
            double ESiD2_=catima::energy_out(particle2,SiD2) * fragment2_mass; // 第二块硅出射能量MeV    
            double dESiD2_=ESiD1_-ESiD2_; // 第二块硅能量沉积MeV 
            if(ESiD2_==0&&ESiD1_!=0) SiCount2++;
            
            particle2.T = ESiD2_/this->fragment2_mass; // setting energy to  MeV/u;
            double ESiD3_=catima::energy_out(particle2,SiD3) * fragment2_mass; // 第三块硅出射能量MeV  
            double dESiD3_=ESiD2_-ESiD3_; // 第三块硅能量沉积MeV 
            if(ESiD3_==0&&ESiD2_!=0) SiCount3++;
            
            particle2.T = ESiD3_/this->fragment2_mass; // setting energy to  MeV/u; 
            double ESiD4_=catima::energy_out(particle2,SiD4) * fragment2_mass; // 第四块硅出射能量MeV  
            double dESiD4_=ESiD3_-ESiD4_; // 第四块硅能量沉积MeV
            if(ESiD4_==0&&ESiD3_!=0) SiCount4++;

            particle2.T = ESiD4_/this->fragment2_mass; // setting energy to  MeV/u; 
            double ESiD5_=catima::energy_out(particle2,SiD5) * fragment2_mass; // 第五块硅出射能量MeV  
            double dESiD5_=ESiD4_-ESiD5_; // 第五块硅能量沉积MeV
            if(ESiD5_==0&&ESiD4_!=0) SiCount5++;

            
            particle2.T = ESiD5_/this->fragment2_mass; // setting energy to  MeV/u; 
            double ECsI_=catima::energy_out(particle2,CsI) * fragment2_mass; // CsI出射能量MeV  
            double dECsI_=ESiD5_-ECsI_; // CsI能量沉积MeV
            if(ECsI_==0) CsICount++;
            // cout<<"dECsI="<<dECsI_<<" dESiD5="<<dESiD5_<<" dESiD4="<<dESiD4_<<" dESiD3="<<dESiD3_
            // <<" dESiD2="<<dESiD2_<<" dESiD1="<<dESiD1_<<" E1_out:"<<E2_out<<endl;


            
            E_out= E1_out+E2_out;//cout<<"E1_out="<<E1_out<<"\t E2_out="<<E2_out<<endl;

            hde2->Fill(E1_out);
            hDepthEnergy->Fill(distance*1e4, E_out);

            hde_E->Fill(dESiD2, dESiD1);
            hde_E->Fill(dESiD2_, dESiD1_);

            hde_E2->Fill(dESiD3, dESiD2);
            hde_E2->Fill(dESiD3_, dESiD2_);

            hde_E3->Fill(dESiD4, dESiD3);
            hde_E3->Fill(dESiD4_, dESiD3_);

            hde_E4->Fill(dESiD5, dESiD4);
            hde_E4->Fill(dESiD5_, dESiD4_);

            hde_E5->Fill(dECsI, dESiD5);
            hde_E5->Fill(dECsI_, dESiD5_);
            
        }

    }
   
    cout<<"SiCount1="<<SiCount1/(2*1.799*1e6)<<"\t SiCount2="<<SiCount2/(2*1.799*1e6)<<"\t SiCount3="
    <<SiCount3/(2*1.799*1e6)<<"\t SiCount4="<<SiCount4/(2*1.799*1e6)<<"\t SiCount5="<<SiCount5/(2*1.799*1e6)
    <<"\t All counts:"<<SiCount1+SiCount2+SiCount3+SiCount4+SiCount5<<endl;



    // 创建画布和3D视图  
    TCanvas* c0 = new TCanvas("c0", "Target Effects Analysis", 1200, 800);
    c0->Divide(2,2);
    c0->cd(1);
    hex->GetXaxis()->SetTitle("Excitation Energy [MeV]"); 
    hex->GetXaxis()->CenterTitle();
    hex->GetXaxis()->SetTitleSize(0.05);
    hex->GetXaxis()->SetTitleOffset(0.90);
    hex->Draw("hist");

    // 绘图1: 碎片总能量分布
    hE_theta->GetXaxis()->SetTitle("#theta [degree]"); 
    hE_theta->GetYaxis()->SetTitle("Kinetic Energy [MeV]");
    hE_theta->GetXaxis()->CenterTitle();
    hE_theta->GetXaxis()->SetTitleSize(0.05);
    hE_theta->GetXaxis()->SetTitleOffset(0.90);
    hE_theta->GetYaxis()->CenterTitle();
    hE_theta->GetYaxis()->SetTitleSize(0.05);
    hE_theta->GetYaxis()->SetTitleOffset(0.90);
    hE_theta->Draw("colz"); // "colz" 表示用颜色填充
    gr->Draw("L"); // "L" 表示画线

    c0->cd(2);
    hex->GetXaxis()->SetTitle("Excitation Energy [MeV]"); 
    hex->GetXaxis()->CenterTitle();
    hex->GetXaxis()->SetTitleSize(0.05);
    hex->GetXaxis()->SetTitleOffset(0.90);
    hex->Draw("hist");

    // 绘图3: 深度对能量损失的影响
    c0->cd(3);
    hde2->GetXaxis()->SetTitle("Total Fragment Energy [MeV]");
    hde2->GetYaxis()->SetTitle("Counts");
    hde2->SetLineColor(kBlue);
    hde2->SetLineWidth(2);
    hde2->Draw();
    
    // 绘图4: 深度与能量关系
    c0->cd(4);
    hDepthEnergy->GetXaxis()->SetTitle("Reaction Depth [um]");
    hDepthEnergy->GetYaxis()->SetTitle("Total Fragment Energy [MeV]");
    hDepthEnergy->GetXaxis()->CenterTitle(); 
    hDepthEnergy->GetYaxis()->CenterTitle();
    hDepthEnergy->GetXaxis()->SetTitleSize(0.05);
    hDepthEnergy->GetXaxis()->SetTitleOffset(0.90);
    hDepthEnergy->GetYaxis()->SetTitleSize(0.05);
    hDepthEnergy->GetYaxis()->SetTitleOffset(1.);
    hDepthEnergy->Draw("COLZ");



    TCanvas* c = new TCanvas("c", "Target Effects Analysis", 1200, 800);
    c->Divide(2,2);
    
    c->cd(1); 
    hde_E2->GetXaxis()->SetTitle("dE3 [MeV]");
    hde_E2->GetYaxis()->SetTitle("dE2 [MeV]");
    hde_E2->GetXaxis()->CenterTitle(); 
    hde_E2->GetYaxis()->CenterTitle();
    hde_E2->GetXaxis()->SetTitleSize(0.05);
    hde_E2->GetXaxis()->SetTitleOffset(0.90);
    hde_E2->GetYaxis()->SetTitleSize(0.05);
    hde_E2->GetYaxis()->SetTitleOffset(1.);
    hde_E2->Draw("colz");

    c->cd(2);
    hde_E3->GetXaxis()->SetTitle("dE4 [MeV]");
    hde_E3->GetYaxis()->SetTitle("dE3 [MeV]");
    hde_E3->GetXaxis()->CenterTitle(); 
    hde_E3->GetYaxis()->CenterTitle();
    hde_E3->GetXaxis()->SetTitleSize(0.05);
    hde_E3->GetXaxis()->SetTitleOffset(0.90);
    hde_E3->GetYaxis()->SetTitleSize(0.05);
    hde_E3->GetYaxis()->SetTitleOffset(1.);
    hde_E3->Draw("colz");

    c->cd(3);
    hde_E4->GetXaxis()->SetTitle("dE5 [MeV]");
    hde_E4->GetYaxis()->SetTitle("dE4 [MeV]");
    hde_E4->GetXaxis()->CenterTitle(); 
    hde_E4->GetYaxis()->CenterTitle();
    hde_E4->GetXaxis()->SetTitleSize(0.05);
    hde_E4->GetXaxis()->SetTitleOffset(0.90);
    hde_E4->GetYaxis()->SetTitleSize(0.05);
    hde_E4->GetYaxis()->SetTitleOffset(1.);
    hde_E4->Draw("colz");

    c->cd(4);

    // hde_E->GetXaxis()->SetTitle("dE1 [MeV]");
    // hde_E->GetYaxis()->SetTitle("dE2 [MeV]");
    // hde_E->GetXaxis()->CenterTitle(); 
    // hde_E->GetYaxis()->CenterTitle();
    // hde_E->GetXaxis()->SetTitleSize(0.05);
    // hde_E->GetXaxis()->SetTitleOffset(0.90);
    // hde_E->GetYaxis()->SetTitleSize(0.05);
    // hde_E->GetYaxis()->SetTitleOffset(1.);
    // hde_E->Draw("colz");

    hde_E5->GetXaxis()->SetTitle("E [MeV]");
    hde_E5->GetYaxis()->SetTitle("dE5 [MeV]");
    hde_E5->GetXaxis()->CenterTitle(); 
    hde_E5->GetYaxis()->CenterTitle();
    hde_E5->GetXaxis()->SetTitleSize(0.05);
    hde_E5->GetXaxis()->SetTitleOffset(0.90);
    hde_E5->GetYaxis()->SetTitleSize(0.05);
    hde_E5->GetYaxis()->SetTitleOffset(1.);
    hde_E5->Draw("colz");
    
    

}// ParticleDynamics::targetDeposit end here!



#endif // RK4_WITH_LEGACY_ANALYSIS

void  ParticleDynamics::SetNuclearData(){// 设置核反应数据


    Int_t residual_z = target_z + projectile_z - ejectile_z;
    Int_t residual_a = target_a + projectile_a - ejectile_a;

    if(residual_z < 0)
        cout<<"Illegal nuclei! Residual Z:"<<residual_z<<endl;
    if (residual_a < 1)
        cout<<"Illegal nuclei! Residual A:"<<residual_a<<endl;

    // from amdc get mass
    this->projectile_mass= GetNuclearMass(projectile_a,projectile_z);
    this->target_mass= GetNuclearMass(target_a,target_z);
    this->ejectile_mass= GetNuclearMass(ejectile_a,ejectile_z);
    this->residual_mass=GetNuclearMass(residual_a,residual_z);

    this->fragment1_mass=GetNuclearMass(fragment1_a,fragment1_z);
    this->fragment2_mass=GetNuclearMass(fragment2_a,fragment2_z);
    this->fragment3_mass=GetNuclearMass(fragment3_a,fragment3_z);

    //this->projectile_energy = projectile_energy * 1e-3; // MeV

    this->proj_energy_start = proj_energy * projectile_mass; // 11.5 MeV  12.29  12.5

    // define target material 
    //this->material.add_element(1,1,4);     // arguments are A,Z, stoichiometric number
    this->material.add_element(1,1,2);     // arguments are A,Z, stoichiometric number
    this->material.density(density);    // density in g/cm3
    double length=100;//cm
    this->material.thickness(length*density);   // thickness in g/cm2

    // define projectile
    catima::Projectile p(this->projectile_a,this->projectile_z,this->projectile_z,proj_energy); 
    p.T = proj_energy; // setting energy to  MeV/u;
    this->projectile=p;
    // EoM 内层循环复用这个对象，避免每次调用都构造 catima::Projectile
    this->proj_eom = catima::Projectile(this->projectile_a, this->projectile_z);


    // Energy at far end of detector
    this->proj_energy_stop = catima::energy_out(this->projectile,this->material) * projectile_mass;

    // cout reaction information 
    cout<<"Reaction:"<<projectile_a<<GetNuclearName(projectile_a,projectile_z)<<"("<<target_a
        <<GetNuclearName(target_a,target_z)<<","<<ejectile_a<<GetNuclearName(ejectile_a,ejectile_z)
        <<")"<<residual_a<<GetNuclearName(residual_a,residual_z)<<endl;;
    cout<<"proj_energy_stop:"<<catima::energy_out(this->projectile,this->material)<<" MeV/u"<<endl;
    cout<<"Beam energy range::"<<this->proj_energy_start<<" to "<<this->proj_energy_stop<<" MeV"<<endl;
    cout<<GetNuclearName(projectile_a,projectile_z);
    printf(" Mass:%0.8lf u\n",projectile_mass);


}// ParticleDynamics::SetNuclearDat end here!

double  ParticleDynamics::GetNuclearMass(Int_t A,Int_t Z){// 获取原子核质量
    // 统一走 MassData 的缓存表；原实现每次调用都重读 2550 行质量表
    return MassData::GetNuclearMassNuclear(A, Z);
}// ParticleDynamics::GetNuclearMassn end here!

TString  ParticleDynamics::GetNuclearName(Int_t A,Int_t Z){// 获取原子核名称
    return TString(MassData::GetNuclearName(A, Z));
}// ParticleDynamics::GetNuclearName end here!

double ParticleDynamics::GetEnergyLoss(double distance){

    this->material.thickness(distance*this->density);   // thickness in g/cm2

    return  catima::energy_out(this->projectile,this->material) * projectile_mass;

}// ParticleDynamics::GetEnergyLoss end here!

double ParticleDynamics::GetBeamEnergy(double distance){

    this->material.thickness(distance*this->density);   // thickness in g/cm2

    return  catima::energy_out(this->projectile,this->material) * projectile_mass;

}// ParticleDynamics::GetBeamEnergy end here!

double ParticleDynamics::GetExciteEnergy(double projectile_ke, double ejectile_ke, double theta){//
	
    const double unit = u_MeV;
	
	double m1 = this->projectile_mass * unit; // incoming nucleus mass 
	double m2 = this->target_mass * unit;     // target nucleus mass
	double m3 = this->ejectile_mass * unit;   // 出射、观测核
	double m4 = this->residual_mass * unit; // residual nucleus 反冲核、剩余核
	
	double E1 = m1 + projectile_ke;
	double E3 = m3 + ejectile_ke;

	double p1 = sqrt(pow(E1,2) - pow(m1,2));
	double p3 = sqrt(pow(E3,2) - pow(m3,2));

	double thetalab = theta;
	
	double s = pow(m1,2) + pow(m2,2) + 2*m2*E1;
	double t = pow(m1,2) + pow(m3,2) + 2*(p1*p3*cos(thetalab) - E1*E3);
	double u = pow(m2,2) + pow(m3,2) - 2*m2*E3;

	double Ex = sqrt(s + t + u - pow(m1,2) - pow(m2,2) - pow(m3,2)) - m4;
	
	return Ex;

}// ParticleDynamics::GetEnergyLoss end here!


double ParticleDynamics::AngleTrans(double projectile_ke,  double ejectile_ke, double theta, double target_ke){
		
	const double unit = u_MeV;
	double m1 = this->projectile_mass * unit; // incoming nucleus mass 
	double m2 = this->target_mass * unit;     // target nucleus mass
	double m3 = this->ejectile_mass * unit;   // 出射、观测核
	double m4 = this->residual_mass * unit;   // residual nucleus 反冲核、剩余核
	
	double E1 = m1 + projectile_ke;
	double E2 = m2 + target_ke;
	double E3 = m3 + ejectile_ke;
	double E4 = E1 + E2 - E3;

	double thetalab = theta;
	
	double s = pow(m1,2)+pow(m2,2)+2*m2*E1;
    double u = pow(m2,2)+pow(m3,2)-2*m2*E3;

    double  omega0=sqrt(s * s + pow(m1, 2) *pow(m1, 2) + pow(m2, 2) * pow(m2, 2) - 2 * s * pow(m1, 2) - 2 * pow(m1, 2) * pow(m2, 2) - 2 * s * pow(m2, 2));
    double  omega1=sqrt(u * u + pow(m2, 2) *pow(m2, 2) + pow(m3, 2) * pow(m3, 2) - 2 * u * pow(m2, 2) - 2 * pow(m2, 2) * pow(m3, 2) - 2 * u * pow(m3, 2));
    

    double m4_ex = sqrt((cos(thetalab) * omega0 * omega1 -(s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) / (2 * pow(m2, 2)) + s + u - pow(m2, 2));

    double   Ex = m4_ex - m4;

    double t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * E4;
	
    double  omega2=sqrt(s * s + pow(m3, 2) *pow(m3, 2) + pow(m4_ex, 2) *pow(m4_ex, 2) - 2 * s * pow(m3, 2) - 2 *pow(m3, 2) * pow(m4_ex, 2) - 2 * s * pow(m4_ex, 2));

	double theta_cm = acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) + (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /(omega0 * omega2));
	
	return  180.*theta_cm/TMath::Pi();

}// ParticleDynamics::TransCM end here!

#ifdef RK4_WITH_LEGACY_ANALYSIS
TGraph* ParticleDynamics::GetRefLine(TGraph* T, double exciteEnergy){

    double theta[1800];// degree

    double term1, term2, term3, term4, Emax, energy[1800];// MeV

    double Q=(projectile_mass + target_mass- ejectile_mass - residual_mass)*u_MeV - exciteEnergy;
 
    theta[0]=1E-3*TMath::Pi()/180.;
    
    term1=sqrt(projectile_mass*ejectile_mass*proj_energy_start)/(ejectile_mass+residual_mass)*cos(theta[0]);
    term2=sqrt(proj_energy_start*((residual_mass - projectile_mass)/(residual_mass + ejectile_mass)) + term1*term1
         + residual_mass*Q/(residual_mass+ejectile_mass)); 
    Emax=(term1+term2)*(term1+term2);
    
    energy[0]=Emax;  //cout<<"energy[0]:"<<energy[0]<<endl;
    T->SetPoint(0, theta[0]/TMath::Pi()*180.,Emax); 

    for(Int_t i=1;i<1800;i++){ //solve Q function
        
        // Lab
        energy[i] = Emax * (1 - i / 1799.); //energy[i] = Emax *(1- exp(-(1800.-i)/30000)); // 逐渐减小的能量值  
        term3=(projectile_mass/residual_mass-1)*proj_energy_start + (ejectile_mass/residual_mass+1)*energy[i] - Q;
        term4=2/residual_mass*sqrt(projectile_mass*ejectile_mass*proj_energy_start*energy[i]);
        if(term3/term4<=1)
            theta[i]=180.*acos(term3/term4)/TMath::Pi();   
        else
            continue;

        if (std::isnan(theta[i])){// theta[i] is equal to nan
            theta[i]=0.;
        }
        T->SetPoint(i, theta[i],energy[i]);   
    }
    T->SetLineColor(kRed);
    T->SetLineWidth(2);

    return T;

}// ParticleDynamics::GetRefLine end here!

TGraph* ParticleDynamics::GetRefLine(TGraph* T, TGraph* Tcm, double exciteEnergy){//Function overload


    double theta[1800], thetaCM[1800] ;// degree

    double term1, term2, term3, term4, Emax, energy[1800];// MeV

    double Q=(projectile_mass + target_mass- ejectile_mass - residual_mass)*u_MeV - exciteEnergy;
 
    theta[0]=1E-3*TMath::Pi()/180.;
    
    term1=sqrt(projectile_mass*ejectile_mass*proj_energy_start)/(ejectile_mass+residual_mass)*cos(theta[0]);
    term2=sqrt(proj_energy_start*((residual_mass - projectile_mass)/(residual_mass + ejectile_mass)) + term1*term1
         + residual_mass*Q/(residual_mass+ejectile_mass)); 
    Emax=(term1+term2)*(term1+term2);
    
    energy[0]=Emax;  //cout<<"energy[0]:"<<energy[0]<<endl;
    thetaCM[0]=AngleTrans(proj_energy_start, Emax, theta[0]);// CM
    
    if (std::isnan(thetaCM[0])){
        thetaCM[0]=0.;
    }

    T->SetPoint(0, theta[0]/TMath::Pi()*180.,Emax); 
    Tcm->SetPoint(0, thetaCM[0],Emax);// CM

    for(Int_t i=1;i<1800;i++){ //solve Q function
        
        // Lab
        energy[i] = Emax * (1 - i / 1799.); //energy[i] = Emax *(1- exp(-(1800.-i)/30000)); // 逐渐减小的能量值  
        term3=(projectile_mass/residual_mass-1)*proj_energy_start + (ejectile_mass/residual_mass+1)*energy[i] - Q;
        term4=2/residual_mass*sqrt(projectile_mass*ejectile_mass*proj_energy_start*energy[i]);
        if(term3/term4<=1)
            theta[i]=180.*acos(term3/term4)/TMath::Pi();   
        else
            continue;

        if (std::isnan(theta[i])){// theta[i] is equal to nan
             continue;
        }
        T->SetPoint(i, theta[i],energy[i]);   
       
        // CM
        thetaCM[i]=AngleTrans(proj_energy_start, energy[i], TMath::Pi()*theta[i]/180.);        
        if (std::isnan(thetaCM[i])){// theta[i] is equal to nan
             continue;
        } 
        Tcm->SetPoint(i, thetaCM[i],energy[i]);     

    }
    T->SetLineColor(kRed);
    T->SetLineWidth(2);
    Tcm->SetLineColor(kRed);
    Tcm->SetLineWidth(2);

    return Tcm;

}// ParticleDynamics::GetRefLine end here!

TGraph* ParticleDynamics::GetRefLineRe(TGraph* T, TGraph* Tcm, double exciteEnergy){// Function overload 相对论下质心系转换

    double theta, thetaCM;// degree

    double  exit_energy,recoil_energy,exit_angle,recoil_angle;// MeV

    for(Int_t i=1;i<1800;i++){ 
        
        InelasticScattering(exit_energy,recoil_energy,exit_angle,recoil_angle, 
                         proj_energy_start,exciteEnergy, TMath::Pi()*i/1800.0);
   
        theta=recoil_angle*180.0/TMath::Pi();      //轻核角度
        thetaCM=180-180*i/1800.0;         //计算的是束流粒子的质心系角度  轻核：180-角度
        //cout<<"i:"<<i<<"\t theta:"<<theta<<"\t thetaCM:"<<thetaCM<<"\t exit_energy:"<<exit_energy<<"\t recoil_energy:"<<recoil_energy-target_mass*u_MeV<<endl;
        T->SetPoint(i, theta,recoil_energy-ejectile_mass*u_MeV); 
        Tcm->SetPoint(i, thetaCM,recoil_energy-ejectile_mass*u_MeV);

    }

    T->SetLineColor(kRed);
    T->SetLineWidth(2);
    Tcm->SetLineColor(kRed);
    Tcm->SetLineWidth(2);

    return Tcm;


}// ParticleDynamics::GetRefLineRe end here!

TGraph* ParticleDynamics::GetRefLineS(TGraph* T, TGraph* T2, double exciteEnergy){// Function overload 相对论下质心系转换

    double theta, theta2;// degree

    double  exit_energy,recoil_energy,exit_angle,recoil_angle;// MeV

    for(Int_t i=1;i<1800;i++){ 
        
        this->InelasticScattering(exit_energy,recoil_energy,exit_angle,recoil_angle, 
                         proj_energy_start,exciteEnergy, TMath::Pi()*i/1800.0);
   
        theta=recoil_angle*180.0/TMath::Pi();      //轻核角度
        theta2=exit_angle*180.0/TMath::Pi();      //束流出射粒子角度 
        //cout<<"i:"<<i<<"\t theta:"<<theta<<"\t thetaCM:"<<thetaCM<<"\t exit_energy:"<<exit_energy<<"\t recoil_energy:"<<recoil_energy-target_mass*u_MeV<<endl;
        T->SetPoint(i-1, theta,recoil_energy-ejectile_mass*u_MeV); 
        T2->SetPoint(i-1, theta2,exit_energy-residual_mass*u_MeV);

    }

    T->SetLineColor(kRed);
    T->SetLineWidth(2);
    T2->SetLineColor(kRed);
    T2->SetLineWidth(2);

    return T;


}// ParticleDynamics::GetRefLineS end here!

TGraph* ParticleDynamics::GetBreakupRefLineS(TGraph* T, TGraph* T2,double Ek, double exciteEnergy){// Function overload 相对论下质心系转换

    double theta, theta2;// degree

    double  fragment1_energy,fragment2_energy,fragment1_angle,fragment2_angle;// MeV

    for(Int_t i=1;i<1800;i++){ 
        
        this->BreakupReaction(this->residual_mass,
                              this->fragment1_mass,
                              this->fragment2_mass,
                              Ek,
                              exciteEnergy,
                              TMath::Pi()*i/1800.0,
                              fragment1_energy,
                              fragment2_energy,
                              fragment1_angle,
                              fragment2_angle
             );

   
        theta=fragment1_angle*180.0/TMath::Pi();      //束流出射粒子角度
        theta2=fragment2_angle*180.0/TMath::Pi();      //轻核角度
        //cout<<"i:"<<i<<"\t theta:"<<theta<<"\t theta2:"<<theta2<<"\t fragment1_energy:"
        //    <<fragment1_energy-fragment1_mass*u_MeV<<"\t fragment2_energy:"<<fragment2_energy-fragment2_mass*u_MeV<<endl;
        T->SetPoint(i-1, theta,fragment1_energy-fragment1_mass*u_MeV); 
        T2->SetPoint(i-1, theta2,fragment2_energy-fragment2_mass*u_MeV);

    }

    T->SetLineColor(kRed);
    T->SetLineWidth(2);
    T2->SetLineColor(kRed);
    T2->SetLineWidth(2);

    return T;

}// ParticleDynamics::GetBreakupRefLineS end here!

TGraphErrors* ParticleDynamics::GetCrossSection(TH1D *hTheta,TGraphErrors* grCrossSection){

    // constant
    static  double const Na=6.02214076e23;

    const double length=100.;//cm
    double thickness=this->density*length; // thickness in g/cm2
    double targetParticle=(thickness/this->target_mass)*Na; // Particle Counts : cm-2
    cout<<"targetParticle:"<<targetParticle<<endl;

    // 获取总 bin 数
    Int_t nBins = hTheta->GetNbinsX();
    cout<<"CrossSection nBins:"<<nBins<<endl;
    
    double BinCenterX,all=0.;
    double binContent,binContentError;// 获取每个bin的计数，并计算误差
    double cirArea=0.,dtheta=TMath::Pi()/nBins;
    for(Int_t i=0;i<nBins;i++){

        binContent = hTheta->GetBinContent(i+1);// bin start -> 1  提取每个bin的计数
        binContentError = (binContent > 0) ? 1. / TMath::Sqrt(binContent) : 0;// 防止除零错误
        BinCenterX= hTheta->GetBinCenter(i+1); //提取每个bin左边界的横坐标x。

        cirArea=2*TMath::Pi()*( cos(i*dtheta) - cos((i+1)*dtheta) );
        all=all+cirArea;
        
        binContent=binContent/(particle_in*targetParticle*cirArea);// cm2/sr
        binContent=binContent*1e27;// mb/sr
        binContentError=binContentError*binContent;
   
        grCrossSection->SetPoint(i,BinCenterX,binContent);     
        grCrossSection->SetPointError(i, 0., binContentError);
        //cout<<"theta:"<<(i+1)*dtheta*180/TMath::Pi()<<"\tCrossSection:"<<binContent<<endl;
       
    }//cout<<"allcirArea:"<<all<<endl;
    grCrossSection->SetMarkerColor(4);//蓝色 方块
    grCrossSection->SetMarkerStyle(21);
    grCrossSection->GetXaxis()->SetTitle("#theta [degree]"); 
	grCrossSection->GetYaxis()->SetTitle("Cross Section [mb/sr]");

        
    return grCrossSection;
}// ParticleDynamics::GetCrossSection end here!

void ParticleDynamics::WriteCrossSection(TGraphErrors* grCrossSection, TString filename){

    // 将截面数据输出到txt文件中
    fstream OutFile;
    TString out_filename=filename;
    OutFile.open(out_filename,ios::out);
    OutFile<<"thetaCM [degree]      CrossSection [mb/sr]        statistical error"<<endl;
    cout<<"--------------Output CrossSection File--------------"<<endl;
    
      // 获取点的数量
    Int_t nPoints = grCrossSection->GetN();

    double theta, Y, ey;//ey误差
    double readX, efficiency;
    for (Int_t i = 0; i < nPoints; i++) {
        
        grCrossSection->GetPoint(i, theta, Y); // 获取第 i 个点的 x 和 y 值
        ey = grCrossSection->GetErrorY(i); // 获取第 i 个点的 y 误差
        OutFile<<theta<<"\t\t\t"<<Y<<"\t\t\t"<<ey<<endl;

    }
    
    OutFile.close();
       
}// ParticleDynamics::WriteCrossSection end here!

TGraphErrors* ParticleDynamics::CutCrossSectionData(TGraphErrors* grCrossSection, double Left, double Right){// Cut CrossSectionData

      // 获取点的数量
    Int_t nPoints = grCrossSection->GetN();


    double theta, Y, ey,newbinContent;//ey误差
    for (Int_t i = 0; i < nPoints; i++) {
        
        grCrossSection->GetPoint(i, theta, Y); // 获取第 i 个点的 x 和 y 值
        ey = grCrossSection->GetErrorY(i); // 获取第 i 个点的 y 误差
 
        if(theta<Left || theta>Right){
            grCrossSection->SetPoint(i,theta,0.0); 
            grCrossSection->SetPointError(i, 0, 0.0);     
        }
        

    }
    
    return grCrossSection;

}// ParticleDynamics::CutCrossSectionData end here!


void ParticleDynamics::GetEfficiency(TH1D *hTheta , TString filename){

    // 获取总 bin 数
    Int_t nBins = hTheta->GetNbinsX();
    cout<<"CrossSection nBins:"<<nBins<<endl;

    const double targetParticle=1000000; // Particle Counts : cm-2
    cout<<"targetParticle:"<<targetParticle<<endl;

    fstream OutFile;
    TString out_filename=filename;
    OutFile.open(out_filename,ios::out);
    OutFile<<"--------------Create Efficiency File--------------"<<endl;
    cout<<"--------------Efficiency File--------------"<<endl;
    
    double binContent;// 获取每个bin的计数
    double cirArea=0.,BinCenterX,dtheta=TMath::Pi()/nBins;
    double efficiency;
    for(Int_t i=0;i<nBins;i++){

        binContent = hTheta->GetBinContent(i+1);// bin start -> 1  提取每个bin的计数

        BinCenterX= hTheta->GetBinCenter(i+1); // 提取每个bin中心的横坐标x。

        cirArea=2*TMath::Pi()*(cos(i*dtheta) - cos((i+1)*dtheta)) ;
            
        efficiency=binContent/(targetParticle*cirArea/(4*TMath::Pi()));// 

        if (std::isnan(efficiency)||efficiency==0){// efficiency is equal to nan
            efficiency=0.;
        }
        else if(efficiency>1){
             cout<<"Error Efficiency "<<"i:"<<BinCenterX<<" "<<efficiency<<endl;
             efficiency=1.;
        }
        
        //cout<<"efficiency "<<"i:"<<i<<" "<<efficiency<<endl;
        
        OutFile<<BinCenterX<<"\t"<<efficiency<<endl;
       
    }//cout<<"allcirArea:"<<all<<endl;

    OutFile.close();
       
}// ParticleDynamics::GetEfficiency end here!

void ParticleDynamics::SetEfficiency(TGraphErrors* grCrossSection, TString filename){

    ifstream inf;//读入数据
    bool ifOpen=this->OpenFile(inf,filename);

    if(!ifOpen){ //文件不存时
         cout<<filename<<"-->DWBA Data File does not exist!"<<endl;
    }
    else{
        cout<<"Successfully Read the Efficiency File !"<<endl;
    }

    string ss;
    getline(inf,ss);//读入文件中开始的描述性文本

    const double targetParticle=1000000; // Particle Counts : cm-2
    cout<<"SetEfficiency targetParticle:"<<targetParticle<<endl;

      // 获取点的数量
    Int_t nPoints = grCrossSection->GetN();
    std::cout << "Total points in CrossSection TGraphErrors: " << nPoints << std::endl;

    TGraph* grEff=new TGraph;// 效率
    double theta, Y, ey,newbinContent;//ey误差
    double readX, efficiency;

    Int_t pointCount = 0; // 读取文件直到结束
    while (inf >> readX >> efficiency) {
        grEff->SetPoint(pointCount, readX, efficiency);// 设置效率图的点
        //if(readX!=theta)cout<<"Error: The number of bins in EfficiencyFile does not equal the number of bins in CrossSection TGraph!"<<endl;
        pointCount++;  
    }

    for (Int_t i = 0; i < nPoints; i++) {

        grCrossSection->GetPoint(i, theta, Y); // 获取第 i 个点的 x 和 y 值
        ey = grCrossSection->GetErrorY(i); // 获取第 i 个点的 y 误差
        efficiency=grEff->Eval(theta);// 读取文件的效率值,并插值
        newbinContent=Y/efficiency;// 计算效率
        
        //cout<<"theta:"<<theta<<" CrossSection:"<<Y<<" Efficiency:"<<efficiency<<" newbinContent:"<<newbinContent<<endl;
        
        if (std::isnan(newbinContent)){//  is equal to nan
            newbinContent=0.;
        }

        grCrossSection->SetPoint(i,theta,newbinContent);     
        //grCrossSection->SetPointError(i, 0, ey);

    }

    inf.close();
       
}// ParticleDynamics::SetEfficiency end here!

TGraph*  ParticleDynamics::ReadDWBA(TGraph* grDWBA,TString filename){

    ifstream inf;//读入数据
    this->OpenFile(inf,filename);

    string ss;
    getline(inf,ss);//读入文件中开始的描述性文本
    getline(inf,ss);//读入文件中开始的描述性文本

    double X[180], data[180];
    for(Int_t i=0;i<180;i++){
        
        inf>>X[i]>>data[i];
        grDWBA->SetPoint(i,X[i],data[i]);
        //cout<<X[i]<<" "<<data[i]<<endl;
        if (!inf||X[i]==180) break;
    }

    grDWBA->SetMarkerColor(2);//红色 圆点
    grDWBA->SetLineColor(kRed);
    grDWBA->SetLineWidth(2);
    grDWBA->SetMarkerStyle(20);
    grDWBA->GetXaxis()->SetTitle("theta [degree]"); 
	grDWBA->GetYaxis()->SetTitle("Cross Section [mb/sr]");
    grDWBA->GetXaxis()->CenterTitle(); 
	grDWBA->GetYaxis()->CenterTitle();
    grDWBA->GetXaxis()->SetTitleSize(0.05);
    grDWBA->GetXaxis()->SetTitleOffset(0.90);
    grDWBA->GetYaxis()->SetTitleSize(0.05);
    grDWBA->GetYaxis()->SetTitleOffset(0.9);

    return grDWBA;
}// ParticleDynamics::ReadDWBA end here!

double ParticleDynamics::Normalize(TGraph* theory, TGraph* experiment,Int_t start,Int_t end){// 拟合 获取归一化刻度系数

  TF1* fitFun= new TF1("fit",Response, 10, 140, 2);// 使用DWBA的数据创建拟合函数TH1
  fitFun->SetParNames("Nor factor", "Deviation");

  // 设置参数限制
  fitFun->SetParLimits(0, 0.0, 10.0); // 限制 par[0]: Pulse Height 的范围
  fitFun->SetParLimits(1, -1, 1); // 限制 par[1]: Deviation 的范围

  experiment->Fit(fitFun,"NML","",start, end);//QNML  N: No fit; 不绘制拟合曲线。
  double* fitPar;
  fitPar=fitFun->GetParameters();


  cout<<"Pulse Height Normalization Factor : "<<fitPar[0]<<"\t Deviation:"<<fitPar[1]<<endl;

  return fitPar[0];

}// ParticleDynamics:::Normalize end here!

double ParticleDynamics::Response(double *x, double *par) {

    ///par[0]: pulse Height --> Normalize
	///par[1]: deviation

    double adjustedX = x[0] - par[1];// 调左右偏移量

    // 在 TGraph 中查找对应的 Y 值
    double val = 0.;
    if (adjustedX > RangeStart && adjustedX < RangeEnd) {
        // 使用 TGraph 的 Eval 方法获取 Y 值
        val = par[0] * grDWBA->Eval(adjustedX);
    }

    return val;

}// ParticleDynamics:::Response end here!


double ParticleDynamics::twoGaus(double *x,double *par){ // 双高斯函数+一次函数
    
    // par[0]: 常数项，表示基线偏移
    // par[1]: 线性项的斜率
    // par[2]: 第一个高斯分布的幅度
    // par[3]: 第一个高斯分布的均值
    // par[4]: 第一个高斯分布的标准差
    // par[5]: 第二个高斯分布的幅度
    // par[6]: 第二个高斯分布的均值
    // par[7]: 第二个高斯分布的标准差

  return par[0]+x[0]*par[1]+par[2]*TMath::Exp(-0.5*((x[0]-par[3])/par[4])*((x[0]-par[3])/par[4]))
          +par[5]*TMath::Exp(-0.5*((x[0]-par[6])/par[7])*((x[0]-par[6])/par[7]));

}// ParticleDynamics:::twoGaus end here!
 

void ParticleDynamics::CalclatePID(Int_t m_A , Int_t m_Z){// 粒子鉴别

    double E=0.;// 总能量 MeV
    double E0=projectile_mass*u_MeV;// 静止能量 MeV
    double Ek=this->proj_energy_start;// 总动能 MeV
    double P=0.;// 动量

    E = Ek + E0;
    P=sqrt(E*E-E0*E0);

    double Bp=P/this->projectile_z;// 磁刚度

    double MZ2=projectile_a*projectile_z*projectile_z;
    double dE=MZ2/Ek;

     cout<<"beam---->"<<projectile_a<<GetNuclearName(projectile_a,projectile_z)<<" Bp:"<<Bp<<" [arb]\t Ek:"
        <<Ek/projectile_mass<<" MeV/u"<<"\tdE:"<<dE<<" [arb]"<<endl;

    cout<<"由IC的幅度比例关系推断目标核 4/4 能损dE:"<<dE<<" [arb]"<<endl;

    TString cwd = gSystem->pwd(); //获取当前系统路径
    TString fInName="/amdc_2020";// 原子核质量评估文件
    ifstream inf;//读入数据
    inf.open(cwd+fInName,ios::in);
    if(!inf.is_open()){ //文件不存时
            std::cout << "Current working directory: " << cwd << std::endl;
            cout<<cwd+fInName<<"-->amdc Data File does not exist!"<<endl;
    }
    else{
           ;// cout<<"Successfully Read the Data File !"<<endl;
    }
    string NuclearName;
    Int_t  NuclearA, NuclearZ;
    double NuclearMass;
    string ss;
    getline(inf,ss);//读入文件中开始的描述性文本
    for(Int_t i=0;i<100;i++){// read AMS data      
        
        inf>>NuclearZ>>NuclearA>>NuclearName>>NuclearMass;
        NuclearMass=NuclearMass-NuclearZ*me; // 减去电子质量

        E0=NuclearMass*u_MeV;// 静止能量 MeV
        P=Bp*NuclearZ;// 使用Bp反解动量
        Ek=sqrt(P*P+E0*E0)-E0;// 总动能 MeV

        MZ2=NuclearA*NuclearZ*NuclearZ;
        dE=MZ2/Ek;

        cout<<NuclearA<<GetNuclearName(NuclearA,NuclearZ)<<" Bp:"<<Bp<<" [arb]\t Ek:"
        <<Ek/NuclearMass<<" MeV/u"<<"\tdE:"<<dE<<" [arb]"<<endl;
                      
    } 
    
}// ParticleDynamics::CalclatePID end here!


double ParticleDynamics::GausIntegral(double A, double mu, double sigma, double x1, double x2, double binWidth) {

    // 计算给定高斯函数参数的区间积分
    if(A==0||sigma==0) return 0.;
    
    // 定义高斯函数，直接在栈上创建 TF1 对象
    TF1 gaussian("gaussian", "[0]*exp(-0.5*((x-[1])/[2])^2)", x1, x2);

    // 设置高斯函数的参数：A (振幅), mu (均值), sigma (标准差)
    gaussian.SetParameters(A, mu, sigma);

    // 计算给定区间 [x1, x2] 内的积分（即总计数）
    double integral = gaussian.Integral(x1, x2);
    
    // 根据 bin 宽度调整积分结果（除以 binWidth 以确保与 bin 的宽度成反比）
    integral /= binWidth;  // 将积分结果除以 bin 宽度
    
    return integral;

}// ParticleDynamics::GausIntegral end here!

void ParticleDynamics::minsGaus(TH1D *hTheta ,double x1, double x2, TString filename){
    
    // 将参数减去给定高斯函数参数的区间积分

    ifstream inf;//读入数据
    this->OpenFile(inf,filename);

    double right[180],A[180],mu[180],sigma[180];
    for(Int_t i=0;i<180;i++){// read data      
        inf>>right[i]>>A[i]>>mu[i]>>sigma[i];
        if (!inf) break; //
    } 

    // 获取总 bin 数
    Int_t nBins = hTheta->GetNbinsX();
    // 获取总 bin 宽度
    TH1D* hex= (TH1D*) gROOT->FindObject("hex");
    double exBinBinWidth =hex->GetBinWidth(1); // bin width;
    cout<<"CrossSection nBins:"<<nBins<<" bin width:"<<exBinBinWidth<<endl;
    
    double binContent, binRightEdge,newBinContent;
    for(Int_t i=0;i<nBins;i++){

        binContent = hTheta->GetBinContent(i+1);// bin start -> 1  提取每个bin的计数
        binRightEdge =hTheta->GetBinLowEdge(i + 1) + hTheta->GetBinWidth(i + 1); //提取每个bin右边界的横坐标x。
        
        // 遍历 right 数组，检查右边界是否与给定的 right 数组相匹配
        for (Int_t j = 0; j < 180; j++) {
            if (TMath::Abs(binRightEdge - right[j]) < 1e-5) { // 判断误差范围内相等
                // 如果 bin 的右边界与 right 数组中的某个值相同，则重新赋值该 bin
                newBinContent=this->GausIntegral(A[j],mu[j],sigma[j],x1,x2,exBinBinWidth);

                hTheta->SetBinContent(i + 1,binContent - newBinContent);
                std::cout << "Bin " << i + 1 << " right edge " << binRightEdge 
                          << " matched with right[" << j << "] = " << right[j] 
                          << ", bin content:"<<binContent<<" new bin content set to " 
                          << binContent - newBinContent << std::endl;
                break;
            }
        }
   
    }

}// ParticleDynamics::minsGaus end here!

void ParticleDynamics::SetGausPar(TF1* fgaus){// 对双高斯函数的参数设定拟合初值

       double par[8] = {
            0.0,  // par[0]: 常数项（基线）
            0.0,   // par[1]: 线性斜率（假设为0）
            50.0, // par[2]: 第一个高斯的幅度
            2.1,   // par[3]: 第一个高斯的均值
            1.0,   // par[4]: 第一个高斯的标准差
            10.0,  // par[5]: 第二个高斯的幅度
            2.71,  // par[6]: 第二个高斯的均值
            1.5    // par[7]: 第二个高斯的标准差
        };

        fgaus->SetParName(0, "A0");     // Baseline
        fgaus->SetParName(1, "Slope");  // 线性斜率
        fgaus->SetParName(2, "A1");     // 第一个高斯幅度
        fgaus->SetParName(3, "mu1");    // 第一个高斯均值
        fgaus->SetParName(4, "sigma1"); // 第一个高斯标准差
        fgaus->SetParName(5, "A2");     // 第二个高斯幅度
        fgaus->SetParName(6, "mu2");    // 第二个高斯均值
        fgaus->SetParName(7, "sigma2"); // 第二个高斯标准差
        fgaus->SetParLimits(0, 0.001, 0.001);  // par限制
        fgaus->SetParLimits(1,  0.001, 0.001); // par限制
        fgaus->SetParLimits(3, 2.0, 2.2);  // par限制
        fgaus->SetParLimits(4, 0.16, 0.18);  // par限制
        fgaus->SetParLimits(6, 2.6, 2.9); // par限制
        fgaus->SetParLimits(7, 0.16, 0.18); // par限制

        fgaus->SetParameters(par);

}// ParticleDynamics::SetGausPar end here!

bool ParticleDynamics::OpenFile(ifstream &inf, TString filename){// 打开文件

    TString cwd = gSystem->pwd(); //获取当前系统路径
    TString fInName=filename;
    char firstChar = fInName[0]; // 获取文件名的第一个字符

    if(firstChar == '/') // 判断第一个字符 TString firstChar
        inf.open(fInName,ios::in);
    else
        inf.open(cwd+"/"+fInName,ios::in);

    if(!inf.is_open()){ //文件不存时
            std::cout << "Current working directory: " << cwd << std::endl;
            std::cout<<cwd+fInName<<"--> Data File does not exist!"<<std::endl;
    }
    else{
           ;// cout<<"Successfully Read the Data File !"<<endl;
    }

    return inf.is_open();

}// ParticleDynamics::OpenFile end here!

void ParticleDynamics::AnaThetaCM(TTree* tree,double Emean, TString filename){// 对质心系角度的分bin拟合

    fstream OutFile;
    TString out_filename=filename;
    OutFile.open(out_filename+".txt",ios::out);
    OutFile<<"--------------Create AnaThetaCM File--------------"<<endl;
    cout<<"--------------Create AnaThetaCM File--------------"<<endl;

    
    // 创建一个新的 TCanvas，并绘制
    TCanvas* canCM = new TCanvas;
    
    TH1D* hex= (TH1D*) gROOT->FindObject("hex");
    TH1D* hBepolar= (TH1D*) gROOT->FindObject("hBepolar");// 获取截面bin的大小
  
    
    TTree* treeCM;
    double minAngle = 0.;
    double maxAngle = 0.;
    TString thetaCUT,excitationCUT;;

    // 获取截面总 bin 数
    Int_t nBins = hBepolar->GetNbinsX();// cout<<"nBins="<<nBins<<endl;

    for(Int_t i=-1;i<nBins;i++){// i=-1 时画不卡角度的激发能谱
        
        // 拷贝tree到新的treeCM
        treeCM= tree->CopyTree("");
        //treeCM= tree->CopyTree("ke<18");
        excitationCUT= Form("%.3f-4*0.1775  < excitation && excitation <%.3f+4*0.1775",Emean,Emean);//默认4*sigma
        treeCM= treeCM->CopyTree(excitationCUT);  // 宽条件CUT
        hex->GetXaxis()->SetRangeUser(-1, 8.);
        
        // 计算质心系角度
        minAngle=i*180./nBins;
        maxAngle=(i+1)*180./nBins;
        if(i==-1){
            minAngle=0; 
            maxAngle=180;
        }
        thetaCUT = Form("%.1f < theta_cm && theta_cm < %.1f", minAngle, maxAngle);// 使用 Form 动态生成条件字符串
        
        // cut质子角度
        treeCM= treeCM->CopyTree(thetaCUT);  cout<<"thetaCUT: "<<thetaCUT<<endl;
        this->loopTree(treeCM);
        
        // 画图：删除旧的 canCM，释放内存
        if (canCM) delete canCM;  
        canCM=new TCanvas;
        hex->Draw();

        
        hex->Fit(gaus,"QML","",-2,10);
        //获取拟合参数（幅度、均值、方差）
        double amplitude = gaus->GetParameter(0);  // 幅度
        double mean = gaus->GetParameter(1);       // 均值
        double sigma = gaus->GetParameter(2);      // 方差
        if(hex->GetEntries()==0||amplitude<0||sigma<0) {
            amplitude=0.;
            mean =0.;
            sigma =0.;
        }

        // static TF1* fitgaus=new TF1("fgaus",twoGaus,-1,10,8);
        // this->SetGausPar(fitgaus);// 对双高斯函数的参数设定拟合初值   
        // hex->Fit(fitgaus,"QMW","",1,6);
        // //获取拟合参数（幅度、均值、方差）
        // double amplitude = fitgaus->GetParameter(2);  // 幅度
        // double mean = fitgaus->GetParameter(3);       // 均值
        // double sigma = fitgaus->GetParameter(4);      // 方差
        // if(hex->GetEntries()==0||amplitude<0||sigma<0) {
        //     amplitude=0.;
        //     mean =0.;
        //     sigma =0.;
        // }
        
        if(i!=-1){
            
            // 将拟合结果输出到文本文件
            OutFile << "Bin: " << minAngle+(maxAngle-minAngle)/2. 
                    << "\t\tA: " << amplitude 
                    << "\t\tmu: " << mean 
                    << "\t\tSigma: " << sigma << endl;

        }
        

        // 使用 TLatex 给图加上标签
        hex->SetTitle(thetaCUT);
        TLatex latex;
        double x = 0.2*(hex->GetXaxis()->GetXmax()-hex->GetXaxis()->GetXmin());
        double y = 0.9*(hex->GetMaximum()-hex->GetMinimum());
        //latex.DrawLatex(x, y, Form("Theta_cm: %.1f - %.1f", minAngle, maxAngle));  // 画出角度范围标签

        if(i==-1)canCM->Print(filename+".pdf(");
        if(i>-1&&i<(nBins-1))canCM->Print(filename+".pdf");
        if(i==(nBins-1))canCM->Print(filename+".pdf)");
    }
    
    // 关闭文件
    OutFile.close();


}// ParticleDynamics::AnaThetaCM end here!

TH1D* ParticleDynamics::AnaThetaCrossSection(TTree* tree,double Emean,TString filename) {

    
    ifstream inf;//读入数据
    this->OpenFile(inf,filename + ".txt");

    string ss;
    getline(inf,ss);//读入文件中开始的描述性文本

    double bin[180],A[180],mu[180],sigma[180];
    for(Int_t i=0;i<180;i++){// read data      
        inf>>ss>>bin[i]>>ss>>A[i]>>ss>>mu[i]>>ss>>sigma[i];
        if (!inf) break; //
    } 

    // 关闭文件
    inf.close();

    TH1D* hex= (TH1D*) gROOT->FindObject("hex");
    TH1D* hBepolar= (TH1D*) gROOT->FindObject("hBepolar");// 获取截面bin的大小
    TH1D* hBepolarAdd= (TH1D*) gROOT->FindObject("hBepolarAdd");// 获取截面bin的大小
    hBepolarAdd->Reset();// 先清空数据


    TTree* treeCM;
    double minAngle = 0.;
    double maxAngle = 0.;
    TString thetaCUT, excitationCUT;

    // 获取截面总 bin 数
    Int_t nBins = hBepolar->GetNbinsX();// cout<<"nBins="<<nBins<<endl;

    for(Int_t i=0;i<nBins;i++){
        
        // 拷贝tree到新的treeCM
        treeCM= tree->CopyTree("");
        //treeCM= tree->CopyTree("ke<18");

        // 激发能CUT
        if(mu[i]<Emean-6*0.1775||mu[i]>Emean+6*0.1775) mu[i]=Emean;
        excitationCUT= Form("%.3f-3*0.1775  < excitation && excitation <%.3f+3*0.1775",mu[i],mu[i]);// cout<<"excitationCUT: "<<excitationCUT<<endl;
        treeCM= treeCM->CopyTree(excitationCUT);  // 细化CUT条件

        // 计算质心系角度
        minAngle=i*180./nBins;
        maxAngle=(i+1)*180./nBins;
        thetaCUT = Form("%.1f < theta_cm && theta_cm < %.1f", minAngle, maxAngle);// 使用 Form 动态生成条件字符串
        
        // cut质子角度
        treeCM= treeCM->CopyTree(thetaCUT); // cout<<"thetaCUT: "<<thetaCUT<<endl;
        this->loopTree(treeCM);
        hBepolarAdd->Add(hBepolar);// 将每个bin的直方图相加
        
       
    }
    
    // 将相加后的直方图赋给hBepolar
    hBepolar->Reset(); 
    hBepolar = (TH1D*)hBepolarAdd->Clone("hBepolar");  // 创建一个深拷贝并赋值给 hBepolar
   
   return hBepolar;

}// ParticleDynamics::AnaThetaCrossSection end here!


TH1D* ParticleDynamics::minsPileUp(TH1D *hTheta ,double Emean, TString Emeanfilename ,TString filename){
    
    // 考虑分角度后的计数还有另一个激发态导致的pileup，计算考虑pileup的微分截面
    // Emeanfilename为分析的激发态文件拟合，故filename应为另一个峰的拟合文件
    // 注意：分析过程中，不要改变bin的大小！！！
    // Emean：参考的激发态能量


    ifstream Einf;//读入数据
    string ss;
    this->OpenFile(Einf,Emeanfilename + ".txt");
    getline(Einf,ss);//读入文件中开始的描述性文本

    double Ebin[180],EA[180],Emu[180],Esigma[180];
    for(Int_t i=0;i<180;i++){// read data
        Einf>>ss>>Ebin[i]>>ss>>EA[i]>>ss>>Emu[i]>>ss>>Esigma[i];
        if(Emu[i]==0)Emu[i]=Emean;// 如果文件的mu为0，则赋值为Emean
        if(!Einf) break; //
    }
    // 关闭文件
    Einf.close();


    ifstream inf;//读入数据-->待减去的峰成分
    this->OpenFile(inf,filename + ".txt");
    getline(inf,ss);//读入文件中开始的描述性文本

    double bin[180],A[180],mu[180],sigma[180];
    for(Int_t i=0;i<180;i++){// read data      
        inf>>ss>>bin[i]>>ss>>A[i]>>ss>>mu[i]>>ss>>sigma[i];
        bin[i]=180.0-bin[i];// 质心系：质子与Be-12角度互为补角
        if (!inf) break; //
    } 
    // 关闭文件
    inf.close();


    // 获取总 bin 数
    Int_t nBins = hTheta->GetNbinsX();
    // 获取总 bin 宽度
    TH1D* hex= (TH1D*) gROOT->FindObject("hex");
    double exBinBinWidth =hex->GetBinWidth(1); // bin width;
    cout<<"CrossSection nBins:"<<nBins<<" bin width:"<<exBinBinWidth<<endl;
    
    double binContent, BinCenterX,newBinContent;
    for(Int_t i=0;i<nBins;i++){

        binContent = hTheta->GetBinContent(i+1);// bin start -> 1  提取每个bin的计数
        BinCenterX =hTheta->GetBinCenter(i+1); // 提取每个bin中心的横坐标x
        
        // 遍历 bin 数组，检查是否与给定的 bin 数组相匹配
        for (Int_t j = 0; j < 180; j++) {
            if (TMath::Abs(BinCenterX - bin[j]) < 1e-5) { // 判断误差范围内相等
                // 如果 bin 的右边界与 bin 数组中的某个值相同，则重新赋值该 bin
                //newBinContent=this->GausIntegral(A[j],mu[j],sigma[j],Emu[j]-3*Esigma[j],Emu[j]+3*Esigma[j],exBinBinWidth);// 3*sigma的积分
                newBinContent=this->GausIntegral(A[j],mu[j],sigma[j],Emu[j]-3*0.1775,Emu[j]+3*0.1775,exBinBinWidth);// 3*sigma的积分
                newBinContent=binContent - newBinContent; // 减去3sigma积分
                if(newBinContent<0) newBinContent=0. ;// 注意！！！！
                hTheta->SetBinContent(i + 1,newBinContent);
                std::cout << "hTheta BinCenterX: " << BinCenterX << " binContent: " << binContent 
                          << "\t matched with theta bin[" << j*180./nBins << "] \t-> " << bin[j] 
                           << " newBinContent: " << newBinContent << std::endl;
                break; // 找到匹配的 bin 后，跳出循环
            }
        }
   
    }// 此处hTheta虽然是Be-12质心系下的计数，但文件中的Be-12激发能谱拟合参数是在相同质心系角度卡出来的，因此可以直接相减。

    return hTheta;

}// ParticleDynamics::minsPileUp end here!

void ParticleDynamics::minsPileUp(TH1D *hexPerBin,double bin_center_cm, double &minValue, Int_t binValue, Int_t peak){// 函数重载
    
    // 考虑分角度后的计数还有另一个激发态导致的pileup，计算考虑pileup的微分截面
    
    //cout<<"minValue:"<<minValue;
    // 获取总 bin 宽度
    double exBinBinWidth =hexPerBin->GetBinWidth(1); // bin width;
    
    double A, mu, sigma; // 获取拟合参数
    A = this->bgCoeffMap.count(bin_center_cm) ? this->bgCoeffMap[bin_center_cm][peak*3-2] : 0.0;// 从 map 中获取容器值
    mu = this->bgCoeffMap.count(bin_center_cm) ? this->bgCoeffMap[bin_center_cm][peak*3-1] : 0.0; // 从 map 中获取容器值
    sigma = this->bgCoeffMap.count(bin_center_cm) ? this->bgCoeffMap[bin_center_cm][peak*3] : 0.0; // 从 map 中获取容器值

    // 计算高斯函数在 bin 边界处的值
    double newBinContent;
    double low_edge = hexPerBin->GetBinLowEdge(binValue); // bin的低边界
    double high_edge = low_edge + exBinBinWidth; // 低边界加上宽度得到高边界

    newBinContent=this->GausIntegral(A,mu,sigma,low_edge,high_edge,exBinBinWidth);// 积分值

    newBinContent= minValue - newBinContent; // 减去积分值
    
    if(newBinContent<0) 
        newBinContent=minValue;// 注意！！！！

    minValue=newBinContent; // 更新返回值
    //cout<<"  nminValue:"<<minValue<<endl;
}// ParticleDynamics::minsPileUp end here!

TGraph* ParticleDynamics::computeDiff(const TGraph* gr) {
   
    if (!gr) return nullptr; // 处理空指针输入
    Int_t n = gr->GetN();
    if (n < 2) return nullptr; // 点数不足，无法计算差分

    double* x = gr->GetX();
    double* y = gr->GetY();

    Int_t newN = 180;
    double* newX = new double[newN];
    double* newY = new double[newN];

    for (Int_t i = 0; i < newN; ++i) {
        newX[i] = (i+ (i + 1)) / 2.0; // 计算相邻点的中点作为新x坐标
        newY[i] =gr->Eval(i) - gr->Eval(i+1);         // 计算y值的差分
        //cout<<"newX["<<i<<"]:"<<newX[i]<<" newY["<<i<<"]:"<<newY[i]<<endl;
    }

    TGraph* diffGr = new TGraph(newN, newX, newY);
    
    delete[] newX;
    delete[] newY;
    
    diffGr->GetXaxis()->SetTitle("#theta "); 
    diffGr->GetYaxis()->SetTitle("Diff   [MeV]");
    diffGr->GetXaxis()->CenterTitle(); 
    diffGr->GetYaxis()->CenterTitle();
    diffGr->GetXaxis()->SetTitleSize(0.05);
    diffGr->GetXaxis()->SetTitleOffset(0.90);
    diffGr->GetYaxis()->SetTitleSize(0.05);
    diffGr->GetYaxis()->SetTitleOffset(0.9);

    return diffGr;

}// ParticleDynamics::computeDiff end here!


/// @brief simulate process of inelastic scattering
/// @param[in] beam_mass mass of beam particles
/// @param[in] target_mass mass of target particle
/// @param[in] beam_energy full energy of beam particle
/// @param[in] excited_energy excited energy of beam particle after reaction
/// @param[in] theta elastic angle in center of mass coordinate ->beam particle
/// @param[out] exit_energy full energy of excited beam particle
/// @param[out] recoil_energy full energy of recoil particle
/// @param[out] exit_angle angle of excited beam particle in lab frame 出射 重核
/// @param[out] recoil_angle angle of recoil particle in lab frame     反冲核 轻
///
#endif // RK4_WITH_LEGACY_ANALYSIS
void ParticleDynamics::InelasticScattering(
	double &exit_energy,
	double &recoil_energy,
	double &exit_angle,
	double &recoil_angle,
    double beam_energy,
	double excited_energy,
    double theta){

    double beam_mass=this->projectile_mass*u_MeV;
    double target_mass=this->target_mass*u_MeV;
             beam_energy=beam_energy+beam_mass; // total energy
	double beam_momentum = sqrt(pow(beam_energy, 2.0) - pow(beam_mass, 2.0));
	double beta_mass_center = beam_momentum / (beam_energy + target_mass);
	double gamma_mass_center =
		1.0 / sqrt(1.0 - beta_mass_center * beta_mass_center);
	double reaction_energy = sqrt(
		(beam_energy + target_mass + beam_momentum)
		* (beam_energy + target_mass - beam_momentum)
	);
	// mass of exit particle
	double exit_mass = residual_mass*u_MeV + excited_energy;
	// mass of recoil particle
	double recoil_mass = ejectile_mass*u_MeV;
	// momentum of exit particle or recoil particle in center of mass frame
	double exit_momentum_center =
		sqrt(
			(reaction_energy - exit_mass - recoil_mass)
			* (reaction_energy - exit_mass + recoil_mass)
			* (reaction_energy + exit_mass - recoil_mass)
			* (reaction_energy + exit_mass + recoil_mass)
		) / (2.0 * reaction_energy);
	// exit momentum parallel and vertical part
	double exit_momentum_center_parallel = exit_momentum_center * cos(theta);
	double exit_momentum_center_vertical = exit_momentum_center * sin(theta);
	// exit energy in c.m.
	double exit_energy_center = sqrt(
		exit_momentum_center * exit_momentum_center
		+ exit_mass * exit_mass
	);
	// exit energy in lab frame
	exit_energy = gamma_mass_center * exit_energy_center
		+ gamma_mass_center * beta_mass_center * exit_momentum_center_parallel;
	// exit momentum parallel part in lab frame
	double exit_momentum_parallel = gamma_mass_center * exit_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * exit_energy_center;
	// exit momentum vertical part in lab frame
	double exit_momentum_vertical = exit_momentum_center_vertical;
	// calculate exit angle in lab frame
	exit_angle = fabs(atan(exit_momentum_vertical / exit_momentum_parallel));
	exit_angle = exit_momentum_parallel > 0 ?
		exit_angle : TMath::Pi() - exit_angle;

	// p of recoil particle in center of mass frame
	double recoil_momentum_center = exit_momentum_center;
	// recoil momentum parallel and vertical part
	double recoil_momentum_center_parallel = -recoil_momentum_center * cos(theta);
	double recoil_momentum_center_vertical = -recoil_momentum_center * sin(theta);
	// recoil energy in c.m.
	double recoil_energy_center = sqrt(
		recoil_momentum_center * recoil_momentum_center
		+ recoil_mass * recoil_mass
	);
	// recoil energy in lab frame
	recoil_energy = gamma_mass_center * recoil_energy_center
		+ gamma_mass_center * beta_mass_center * recoil_momentum_center_parallel;
	// recoil momentum parallel part in lab frame
	double recoil_momentum_parallel =
		gamma_mass_center * recoil_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * recoil_energy_center;
	// recoil momentum vertical part in lab frame
	double recoil_momentum_vertical = recoil_momentum_center_vertical;
	// recoil angle in lab frame
	recoil_angle =
		fabs(atan(recoil_momentum_vertical / recoil_momentum_parallel));
	recoil_angle = recoil_momentum_parallel > 0 ?
		recoil_angle : TMath::Pi() - recoil_angle;

} // ParticleDynamics::InelasticScattering end here!


/// @brief simulate breakup reaction process
/// @param[in] parent_mass mass of parent particle
/// @param[in] fragment1_mass mass of fragment particle 1
/// @param[in] fragment2_mass mass of fragment particle 2
/// @param[in] parent_energy full energy of parent particle
/// @param[in] theta breakup angle in center of mass coordinate  fragment1
/// @param[out] fragment1_energy full energy of fragment particle 1
/// @param[out] fragment2_energy full energy of fragment particle 2
/// @param[out] fragment1_angle angle of fragment particle 1 in lab frame
/// @param[out] fragment2_angle angle of fragment particle 2 in lab frame
///
void ParticleDynamics::BreakupReaction(
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
) {
    parent_mass=parent_mass*u_MeV;
    fragment1_mass=fragment1_mass*u_MeV;
    fragment2_mass=fragment2_mass*u_MeV;
    double parent_energy=parent_Ek+parent_mass; // total energy of parent particle

    // cout<<"parent mass: "<<parent_mass<<endl;
    // cout<<"fragment1 mass: "<<fragment1_mass<<endl; 
    // cout<<"fragment2 mass: "<<fragment2_mass<<endl;     
    // cout<<"parent energy: "<<parent_energy<<endl;

	// momentum of parent particle in lab frame
	double parent_momentum =
		sqrt(pow(parent_energy, 2.0) - pow(parent_mass, 2.0));
        
	// beta of center of mass
	double beta_center = parent_momentum / parent_energy;
	// gamma of center of mass
	double gamma_center = 1.0 / sqrt(1.0 - pow(beta_center, 2.0));
	// total energy in center of mass frame
	double parent_energy_center = parent_mass+parent_Ex;//静止质量+激发能
  
	// momentum of fragments in center of mass frame
	double fragment_momentum =
		sqrt(
			(parent_energy_center - fragment1_mass - fragment2_mass)
			* (parent_energy_center - fragment1_mass + fragment2_mass)
			* (parent_energy_center + fragment1_mass - fragment2_mass)
			* (parent_energy_center + fragment1_mass + fragment2_mass)
		) / (2.0 * parent_mass);

	// fragment1 momentum parallel part in c.m. frame
	double fragment1_momentum_center_parallel = fragment_momentum*cos(theta);
	// fragment1 momentum vertical part in c.m. frame
	double fragment1_momentum_center_vertical = fragment_momentum*sin(theta);
	// fragment1 energy in c.m. frame
	double fragment1_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment1_mass, 2.0)
	);
	// fragment1 energy in lab frame
	fragment1_energy = gamma_center * fragment1_energy_center
		+ gamma_center * beta_center * fragment1_momentum_center_parallel;
	// fragment1 momentum parallel part in lab frame
	double fragment1_momentum_parallel =
		gamma_center * fragment1_momentum_center_parallel
		+ gamma_center * beta_center * fragment1_energy_center;
	// fragment1 momentum vertical part in lab frame
	double fragment1_momentum_vertical = fragment1_momentum_center_vertical;
	// fragment1 angle in lab frame
	fragment1_angle =
		fabs(atan(fragment1_momentum_vertical / fragment1_momentum_parallel));
	fragment1_angle = fragment1_momentum_parallel > 0 ?
		fragment1_angle : TMath::Pi() - fragment1_angle;

	// fragment2 momentum parallel part in c.m. frame
	double fragment2_momentum_center_parallel = -fragment_momentum*cos(theta);
	// fragment2 momentum vertical part in c.m. frame
	double fragment2_momentum_center_vertical = -fragment_momentum*sin(theta);
	// fragment2 energy in c.m. frame
	double fragment2_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment2_mass, 2.0)
	);
	// fragment2 energy in lab frame
	fragment2_energy = gamma_center * fragment2_energy_center
		+ gamma_center * beta_center * fragment2_momentum_center_parallel;
	// fragment2 momentum parallel part in lab frame
	double fragment2_momentum_parallel =
		gamma_center * fragment2_momentum_center_parallel
		+ gamma_center * beta_center * fragment2_energy_center;
	// fragment2 momentum vertical part in lab frame
	double fragment2_momentum_vertical = fragment2_momentum_center_vertical;
	// fragment2 angle in lab frame
	fragment2_angle =
		fabs(atan(fragment2_momentum_vertical / fragment2_momentum_parallel));
	fragment2_angle = fragment2_momentum_parallel > 0 ?
		fragment2_angle : TMath::Pi() - fragment2_angle;

}// ParticleDynamics::BreakupReaction end here!

// 旋转矩阵函数：将粒子坐标系中的向量转换到世界坐标系
TVector3 ParticleDynamics::ParticleToWorld(const TVector3& vec_particle, double theta, double phi) {

    TRotation R;
    R.RotateY(theta);  // 绕Y轴旋转极角θ
    R.RotateZ(phi);    // 绕Z轴旋转方位角φ
    return R * vec_particle;

}// ParticleDynamics::ParticleToWorld end here!

//<---------------------main function end----------------------->//


//----------------------------------------------------------------------
// 探测器模拟用接口
//----------------------------------------------------------------------
void ParticleDynamics::SetParticle(int A, int Z) {
    std::lock_guard<std::mutex> catima_lock(CatimaSetupMutex());
    ejectile_a = A;
    ejectile_z = Z;
    ejectile_mass = MassData::GetNuclearMassNuclear(A, Z);
    proj_eom = catima::Projectile(A, Z);
}

void ParticleDynamics::ConfigureTarget(int A, int Z, double stoich, double density_gcm3,
                                       double length_cm) {
    ConfigureTarget(std::vector<Element>{{A,Z,stoich}},density_gcm3,length_cm);
}

void ParticleDynamics::ConfigureTarget(const std::vector<Element>& compound,
                                       double density_gcm3, double length_cm) {
    std::lock_guard<std::mutex> catima_lock(CatimaSetupMutex());
    const double scale = Utils::stoichiometricScale(compound);
    if (!(density_gcm3>0) || !std::isfinite(density_gcm3))
        throw TargetError("Invalid target mass density");
    target_a = compound.front().a;
    target_z = compound.front().z;
    density = density_gcm3;

    this->material = catima::Material();
    for (const auto& element : compound)
        this->material.add_element(element.a, element.z, std::max(1.0, element.s * scale));
    this->material.density(density_gcm3);
    this->material.thickness(length_cm * density_gcm3); // g/cm^2

    this->proj_eom = catima::Projectile(ejectile_a, ejectile_z);
}

void ParticleDynamics::Propagate(double ke_mev, const double vtx[3], double theta,
                                 double phi, std::vector<TrackPoint>& track,
                                 double dt, int max_steps) {
    double state[6];
    initialState(ke_mev, vtx[0], vtx[1], vtx[2], theta, phi, state);
    Reset();

    track.clear();
    track.reserve(4096);

    rk4_solver::Integrator<6, ParticleDynamics> integrator(
        *this, &ParticleDynamics::equationMotion, dt);

    double t = 0.0;
    for (int i = 0; i < max_steps; ++i) {
        double x_plus[6];
        if (StopEvent(t, state, x_plus)) {
            break;
        }

        const double gv = std::sqrt(state[3] * state[3] + state[4] * state[4] +
                                    state[5] * state[5]);
        double ke_now = 0.0;
        if (gv > 0.0) {
            const double beta = std::sqrt(gv * gv / (1.0 + gv * gv));
            const double gamma = gv / beta;
            ke_now = ejectile_mass * u_MeV * (gamma - 1.0);
        }

        TrackPoint p;
        p.x = state[0];
        p.y = state[1];
        p.z = state[2];
        p.ke = ke_now;
        track.push_back(p);

        double t_next;
        integrator.step(t, state, t_next, state);
        t = t_next;
    }
}
