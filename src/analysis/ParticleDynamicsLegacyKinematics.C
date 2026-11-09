// Optional ROOT-based analysis routines retained from the original Spyral workflow.
// Build with SPYRALPP_WITH_LEGACY_ANALYSIS=ON to include these APIs.
#include "../sim/ParticleDynamics.h"
#include "../core/DataPath.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

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
