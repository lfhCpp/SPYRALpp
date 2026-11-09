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
