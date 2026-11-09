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
