// 运行方式（在 RK4/ 目录下）：
//   root -l -b -q test/DrawTrack.C
// 前置条件：先构建共享库 build/librk4.so（cmake --build build -j）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "sim/ParticleDynamics.h"
#include "../extern/rk4_solver/include/rk4_solver.hpp"



//---------------------------- static define-----------------------//


//--------------------->main function start<-----------------------//
void DrawTrack(){

     
    // 创建动力学对象
    ParticleDynamics D2;

    // 创建RK4积分器和事件对象
    double dt = 1e-10;       // 时间步长（秒）

    rk4_solver::Integrator<6, ParticleDynamics> 
        integrator(D2, &ParticleDynamics::equationMotion, dt);

    rk4_solver::Event<6, ParticleDynamics> 
        event(D2, &ParticleDynamics::StopEvent);

    TCanvas* can = new TCanvas("can", "Particle Trajectory and Projections", 1200, 900);
   
    // 创建水平参考线 
    TLine* hlineUp = new TLine(0, 28,100., 28);
    hlineUp->SetLineColor(kRed);
    hlineUp->SetLineStyle(2); // 虚线
    hlineUp->SetLineWidth(2);
    TLine* hlineDown = new TLine(0, -28,100., -28);
    hlineDown->SetLineColor(kRed);
    hlineDown->SetLineStyle(2); // 虚线
    hlineDown->SetLineWidth(2);


    double theta, thetaCM;// degree
    double  exit_energy,recoil_energy,exit_angle,recoil_angle;// MeV
    int startAngle=1,stopAngle=179;
    for(int thetaIter=startAngle;thetaIter<=stopAngle;thetaIter++){ // 角度抽样循环
        
        // 创建投影图用的TGraph对象
        TGraph* graphXY = new TGraph();      // XY投影
        TGraph* graphXZ = new TGraph();      // XZ投影  
        TGraph* graphYZ = new TGraph();      // YZ投影
        TGraph2D* graph3D = new TGraph2D();  // 3D轨迹
        
        D2.InelasticScattering(exit_energy,
                            recoil_energy,
                            exit_angle,
                            recoil_angle, 
                            D2.GetInitialBeamEnergy(),
                            0.,//激发能
                            TMath::Pi()*thetaIter/180.0);// 束流粒子质心系角度
    
        theta=recoil_angle*180.0/TMath::Pi();  //轻核角度
        thetaCM=180-180*10/180.0;             //计算的是束流粒子的质心系角度  轻核：180-角度

        //cout<<"recoil_energy:"<<recoil_energy-D2.GetEjectileMassMeV()
        //    <<" MeV\t recoil_angle:"<<recoil_angle*180.0/TMath::Pi()<<" degree"<<endl;

        // 初始条件:
        double x0[6];

        D2.initialState(recoil_energy-D2.GetEjectileMassMeV(), 0.0, 0.0, 0.0,recoil_angle, 0.*2*TMath::Pi()/360.0 ,x0);

        double t0 = 0.0; // 初始时间
        
        // 模拟参数
        const int max_steps = 100000; // 最大总步数
        double t_history[max_steps];
        double x_history[max_steps][6];

        // 运行模拟
        D2.SetInitialConditions();// 重置积分器初始条件
        size_t actual_steps = rk4_solver::loop<max_steps>(
            integrator, event, t0, x0, t_history, x_history, true
        );
        //cout << "实际步数: " << actual_steps << endl;
        

        // 重置对象
        graphXY->Set(0);     // 清空XY投影图的数据点
        graphXZ->Set(0);     // 清空XZ投影图的数据点
        graphYZ->Set(0);     // 清空YZ投影图的数据点
        graph3D->Set(0);     // 清空3D轨迹的数据点

        // 填充轨迹数据
        for (int i = 0; i < actual_steps; i++) {
            double x = x_history[i][0]*100; // x坐标
            double y = x_history[i][1]*100; // y坐标
            double z = x_history[i][2]*100; // z坐标
            graph3D->SetPoint(i, x, z, y);
            graphXY->SetPoint(i, x, y);  // XY平面投影
            graphXZ->SetPoint(i, z, x);  // XZ平面投影
            graphYZ->SetPoint(i, z, y);  // YZ平面投影
        }
        cout<<"Z:"<<x_history[actual_steps-1][2]<<" m"<<endl;

        can->Clear();
        can->Divide(2, 2); // 2行2列，共4个区域
        can->cd(1); // 切换到第一个区域
        gPad->SetTheta(32.49337);
        gPad->SetPhi(241.1136);


        graph3D->SetTitle(";X [cm];Z [cm];Y [cm]");
        graph3D->GetXaxis()->CenterTitle(); 
        graph3D->GetYaxis()->CenterTitle();
        graph3D->GetZaxis()->CenterTitle();
        graph3D->GetXaxis()->SetTitleSize(0.05);
        graph3D->GetYaxis()->SetTitleSize(0.05);
        graph3D->GetZaxis()->SetTitleSize(0.05);
        graph3D->GetXaxis()->SetTitleOffset(1.6);
        graph3D->GetYaxis()->SetTitleOffset(1.6);
        graph3D->GetZaxis()->SetTitleOffset(1.05);
        graph3D->SetMarkerColor(kBlue);
        graph3D->SetMarkerStyle(20);
        graph3D->SetMarkerSize(0.5);
        graph3D->SetLineColor(kRed);
        graph3D->SetLineWidth(2);
        graph3D->Draw("PLine");

        TPaveText* legend = new TPaveText(0.05, 0.95, 0.85, 0.95, "NDC");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        //legend->SetTextAlign(22);
        legend->SetTextSize(0.05);
        legend->AddText(Form("Ejectile:%.1f MeV  %.1f#circ  Residual:#theta_{cm}:%.1f#circ", recoil_energy-D2.GetEjectileMassMeV(), 
        recoil_angle*180.0/TMath::Pi(),double(thetaIter)));
        legend->Draw();



        // ==================== 2. XY平面投影 ====================
        can->cd(2); // 切换到第二个区域
        graphXY->SetTitle("XY Projection;X [cm];Y [cm]");
        graphXY->GetXaxis()->CenterTitle(); 
        graphXY->GetYaxis()->CenterTitle();
        graphXY->GetXaxis()->SetTitleSize(0.05);
        graphXY->GetYaxis()->SetTitleSize(0.05);
        graphXY->GetYaxis()->SetTitleOffset(0.9);
        graphXY->SetMarkerColor(kBlue);
        graphXY->SetMarkerStyle(20);
        graphXY->SetMarkerSize(0.5);
        graphXY->SetLineColor(kBlue);
        graphXY->SetLineWidth(1);
        graphXY->GetXaxis()->SetLimits(-30, 30);
        graphXY->GetYaxis()->SetLimits(-30, 30);
        graphXY->GetXaxis()->SetRangeUser(-30, 30);
        graphXY->GetYaxis()->SetRangeUser(-30, 30);
        graphXY->Draw("APL");

        //  绘制参考圆
        double radius = 28; // 圆半径 (cm)  
        TEllipse* ref_circle = new TEllipse(0.0, 0.0, radius, radius);
        ref_circle->SetLineColor(kRed);
        ref_circle->SetLineStyle(2); // 虚线
        ref_circle->SetLineWidth(2);
        ref_circle->SetFillStyle(0); // 空心
        ref_circle->Draw("same");

        //  绘制参考圆
        double radiusC = 1.5; // 圆半径 (cm)
        TEllipse* ref_circleC = new TEllipse(0.0, 0.0, radiusC, radiusC);
        ref_circleC->SetLineColor(kRed);
        ref_circleC->SetLineStyle(2); // 虚线
        ref_circleC->SetLineWidth(2);
        ref_circleC->SetFillStyle(0); // 空心
        ref_circleC->Draw("same");

        // // 添加起始点和结束点标记
        TMarker* startXY = new TMarker(x_history[0][0], x_history[0][1], 29);
        startXY->SetMarkerColor(kGreen);
        startXY->SetMarkerSize(1.2);
        startXY->Draw("same");

        TMarker* endXY = new TMarker(x_history[actual_steps-1][0], 
                                    x_history[actual_steps-1][1], 29);
        endXY->SetMarkerColor(kRed);
        endXY->SetMarkerSize(1.2);
        endXY->Draw("same");

        // ==================== 3. XZ平面投影 ====================
        can->cd(3); // 切换到第三个区域
        
        graphXZ->SetTitle("XZ Projection;Z [cm];X [cm]");
        graphXZ->GetXaxis()->CenterTitle(); 
        graphXZ->GetYaxis()->CenterTitle();
        graphXZ->GetXaxis()->SetTitleSize(0.05);
        graphXZ->GetYaxis()->SetTitleSize(0.05);
        graphXZ->GetYaxis()->SetTitleOffset(0.9);
        graphXZ->SetMarkerColor(kBlue);
        graphXZ->SetMarkerStyle(20);
        graphXZ->SetMarkerSize(0.5);
        graphXZ->SetLineColor(kBlue);
        graphXZ->SetLineWidth(1);
        graphXZ->GetXaxis()->SetLimits(0, 100);
        graphXZ->GetYaxis()->SetLimits(-30, 30);
        graphXZ->GetXaxis()->SetRangeUser(0, 100);
        graphXZ->GetYaxis()->SetRangeUser(-30, 30);
        graphXZ->Draw("APL");
        hlineUp->Draw("same");
        hlineDown->Draw("same");

   
       

        // 添加起始点和结束点标记
        TMarker* startXZ = new TMarker(x_history[0][2], x_history[0][0], 29);
        startXZ->SetMarkerColor(kGreen);
        startXZ->SetMarkerSize(1.2);
        startXZ->Draw("same");

        TMarker* endXZ = new TMarker(x_history[actual_steps-1][2], 
                                    x_history[actual_steps-1][0], 29);
        endXZ->SetMarkerColor(kRed);
        endXZ->SetMarkerSize(1.2);
        endXZ->Draw("same");

        
        // ==================== 4. YZ平面投影 ====================
        can->cd(4); // 切换到第四个区域
        graphYZ->SetTitle("YZ Projection;Z [cm];Y [cm]");
        graphYZ->GetXaxis()->CenterTitle(); 
        graphYZ->GetYaxis()->CenterTitle();
        graphYZ->GetXaxis()->SetTitleSize(0.05);
        graphYZ->GetYaxis()->SetTitleSize(0.05);
        graphYZ->GetYaxis()->SetTitleOffset(0.9);
        graphYZ->SetMarkerColor(kBlue);
        graphYZ->SetMarkerStyle(20);
        graphYZ->SetMarkerSize(0.5);
        graphYZ->SetLineColor(kBlue);
        graphYZ->SetLineWidth(1);
        graphYZ->GetXaxis()->SetLimits(0, 100);
        graphYZ->GetYaxis()->SetLimits(-30, 30);
        graphYZ->GetXaxis()->SetRangeUser(0, 100);
        graphYZ->GetYaxis()->SetRangeUser(-30, 30);
        graphYZ->Draw("APL");
        hlineUp->Draw("same");
        hlineDown->Draw("same");

        // 添加起始点和结束点标记
        TMarker* startYZ = new TMarker(x_history[0][2], x_history[0][1], 29);
        startYZ->SetMarkerColor(kGreen);
        startYZ->SetMarkerSize(1.2);
        startYZ->Draw("same");

        TMarker* endYZ = new TMarker(x_history[actual_steps-1][2], 
                                    x_history[actual_steps-1][1], 29);      
                                
        endYZ->SetMarkerColor(kRed);
        endYZ->SetMarkerSize(1.2);
        endYZ->Draw("same");

     
        //sleep(1); 
        can->Update();

        if(thetaIter==startAngle) can->Print("Trajectory.pdf(");  
        else if(thetaIter==stopAngle) can->Print("Trajectory.pdf)"); 
        else can->Print("Trajectory.pdf"); 

    }// angle loop end  
    
        
}


void SingleTrack(){

     
    // 创建动力学对象
    ParticleDynamics D2;

    // 创建RK4积分器和事件对象
    double dt = 1e-10;       // 时间步长（秒）

    rk4_solver::Integrator<6, ParticleDynamics> 
        integrator(D2, &ParticleDynamics::equationMotion, dt);

    rk4_solver::Event<6, ParticleDynamics> 
        event(D2, &ParticleDynamics::StopEvent);

    double theta, thetaCM;// degree
    double  exit_energy,recoil_energy,exit_angle,recoil_angle;// MeV

    D2.InelasticScattering(exit_energy,
                         recoil_energy,
                         exit_angle,
                         recoil_angle, 
                         D2.GetInitialBeamEnergy(),
                         0.,//激发能
                         TMath::Pi()*70./180.0);// 束流粒子质心系角度
   
    theta=recoil_angle*180.0/TMath::Pi();  //轻核角度
    thetaCM=180-180*10/180.0;             //计算的是束流粒子的质心系角度  轻核：180-角度

    cout<<"recoil_energy:"<<recoil_energy-D2.GetEjectileMassMeV()
        <<" MeV\t recoil_angle:"<<recoil_angle*180.0/TMath::Pi()<<" degree"<<endl;

    // 初始条件:
    double x0[6];

    D2.initialState(recoil_energy-D2.GetEjectileMassMeV(), 0.0, 0.0, 0.0,recoil_angle, 0.*2*TMath::Pi()/360.0 ,x0);

    double t0 = 0.0; // 初始时间
    
    // 模拟参数
    const int max_steps = 100000; // 最大总步数
    double t_history[max_steps];
    double x_history[max_steps][6];

    // 运行模拟
    size_t actual_steps = rk4_solver::loop<max_steps>(
        integrator, event, t0, x0, t_history, x_history, true
    );
    //cout << "实际步数: " << actual_steps << endl;

    // 创建投影图用的TGraph对象
    TGraph* graphXY = new TGraph();      // XY投影
    TGraph* graphXZ = new TGraph();      // XZ投影  
    TGraph* graphYZ = new TGraph();      // YZ投影
    TGraph2D* graph3D = new TGraph2D();  // 3D轨迹
    
    // 填充轨迹数据
    for (int i = 0; i < actual_steps; i++) {
        double x = x_history[i][0]; // x坐标
        double y = x_history[i][1]; // y坐标
        double z = x_history[i][2]; // z坐标
        graph3D->SetPoint(i, x, z, y);
        graphXY->SetPoint(i, x, y);  // XY平面投影
        graphXZ->SetPoint(i, z, x);  // XZ平面投影
        graphYZ->SetPoint(i, z, y);  // YZ平面投影
    }
    cout<<"Z:"<<x_history[actual_steps-1][2]<<" m"<<endl;

    TCanvas* c1 = new TCanvas("c1", "Particle Trajectory and Projections", 1200, 900);
    c1->Divide(2, 2); // 2行2列，共4个区域

    c1->cd(1); // 切换到第一个区域
    gPad->SetTheta(32.49337);
    gPad->SetPhi(241.1136);

    graph3D->SetTitle("3D Particle Trajectory;X [m];Z [m];Y [m]");
    graph3D->GetXaxis()->CenterTitle(); 
    graph3D->GetYaxis()->CenterTitle();
    graph3D->GetZaxis()->CenterTitle();
    graph3D->GetXaxis()->SetTitleSize(0.05);
    graph3D->GetYaxis()->SetTitleSize(0.05);
    graph3D->GetZaxis()->SetTitleSize(0.05);
    graph3D->GetXaxis()->SetTitleOffset(1.6);
    graph3D->GetYaxis()->SetTitleOffset(1.6);
    graph3D->GetZaxis()->SetTitleOffset(1.05);
    graph3D->SetMarkerColor(kBlue);
    graph3D->SetMarkerStyle(20);
    graph3D->SetMarkerSize(0.5);
    graph3D->SetLineColor(kRed);
    graph3D->SetLineWidth(2);
    graph3D->Draw("PLine");



    // ==================== 2. XY平面投影 ====================
    c1->cd(2); // 切换到第二个区域
    graphXY->SetTitle("XY Projection;X [m];Y [m]");
    graphXY->GetXaxis()->CenterTitle(); 
    graphXY->GetYaxis()->CenterTitle();
    graphXY->GetXaxis()->SetTitleSize(0.05);
    graphXY->GetYaxis()->SetTitleSize(0.05);
    graphXY->GetYaxis()->SetTitleOffset(0.9);
    graphXY->SetMarkerColor(kBlue);
    graphXY->SetMarkerStyle(20);
    graphXY->SetMarkerSize(0.5);
    graphXY->SetLineColor(kBlue);
    graphXY->SetLineWidth(1);
    graphXY->GetXaxis()->SetLimits(-0.3, 0.3);
    graphXY->GetYaxis()->SetLimits(-0.3, 0.3);
    graphXY->GetXaxis()->SetRangeUser(-0.3, 0.3);
    graphXY->GetYaxis()->SetRangeUser(-0.3, 0.3);
    graphXY->Draw("APL");

    //  绘制参考圆
    double radius = 0.25; // 圆半径 (m)
    TEllipse* ref_circle = new TEllipse(0.0, 0.0, radius, radius);
    ref_circle->SetLineColor(kRed);
    ref_circle->SetLineStyle(2); // 虚线
    ref_circle->SetLineWidth(2);
    ref_circle->SetFillStyle(0); // 空心
    ref_circle->Draw("same");

    //  绘制参考圆
    double radiusC = 0.015; // 圆半径 (m)
    TEllipse* ref_circleC = new TEllipse(0.0, 0.0, radiusC, radiusC);
    ref_circleC->SetLineColor(kRed);
    ref_circleC->SetLineStyle(2); // 虚线
    ref_circleC->SetLineWidth(2);
    ref_circleC->SetFillStyle(0); // 空心
    ref_circleC->Draw("same");

    // // 添加起始点和结束点标记
    TMarker* startXY = new TMarker(x_history[0][0], x_history[0][1], 29);
    startXY->SetMarkerColor(kGreen);
    startXY->SetMarkerSize(1.2);
    startXY->Draw("same");

    TMarker* endXY = new TMarker(x_history[actual_steps-1][0], 
                                x_history[actual_steps-1][1], 29);
    endXY->SetMarkerColor(kRed);
    endXY->SetMarkerSize(1.2);
    endXY->Draw("same");

    // ==================== 3. XZ平面投影 ====================
    c1->cd(3); // 切换到第三个区域
    graphXZ->SetTitle("XZ Projection;Z [m];X [m]");
    graphXZ->GetXaxis()->CenterTitle(); 
    graphXZ->GetYaxis()->CenterTitle();
    graphXZ->GetXaxis()->SetTitleSize(0.05);
    graphXZ->GetYaxis()->SetTitleSize(0.05);
    graphXZ->GetYaxis()->SetTitleOffset(0.9);
    graphXZ->SetMarkerColor(kBlue);
    graphXZ->SetMarkerStyle(20);
    graphXZ->SetMarkerSize(0.5);
    graphXZ->SetLineColor(kBlue);
    graphXZ->SetLineWidth(1);
    graphXZ->Draw("APL");

    // 创建水平参考线 
    TLine* hlineUp = new TLine(0, 0.25,1., 0.25);
    hlineUp->SetLineColor(kRed);
    hlineUp->SetLineStyle(2); // 虚线
    hlineUp->SetLineWidth(2);
    hlineUp->Draw("same");

    TLine* hlineDown = new TLine(0, -0.25,1., -0.25);
    hlineDown->SetLineColor(kRed);
    hlineDown->SetLineStyle(2); // 虚线
    hlineDown->SetLineWidth(2);
    hlineDown->Draw("same");

    // 添加起始点和结束点标记
    TMarker* startXZ = new TMarker(x_history[0][2], x_history[0][0], 29);
    startXZ->SetMarkerColor(kGreen);
    startXZ->SetMarkerSize(1.2);
    startXZ->Draw("same");

    TMarker* endXZ = new TMarker(x_history[actual_steps-1][2], 
                                x_history[actual_steps-1][0], 29);
    endXZ->SetMarkerColor(kRed);
    endXZ->SetMarkerSize(1.2);
    endXZ->Draw("same");

    
    // ==================== 4. YZ平面投影 ====================
    c1->cd(4); // 切换到第四个区域
    graphYZ->SetTitle("YZ Projection;Z [m];Y [m]");
    graphYZ->GetXaxis()->CenterTitle(); 
    graphYZ->GetYaxis()->CenterTitle();
    graphYZ->GetXaxis()->SetTitleSize(0.05);
    graphYZ->GetYaxis()->SetTitleSize(0.05);
    graphYZ->GetYaxis()->SetTitleOffset(0.9);
    graphYZ->SetMarkerColor(kBlue);
    graphYZ->SetMarkerStyle(20);
    graphYZ->SetMarkerSize(0.5);
    graphYZ->SetLineColor(kBlue);
    graphYZ->SetLineWidth(1);
    graphYZ->Draw("APL");
    hlineUp->Draw("same");
    hlineDown->Draw("same");

    // 添加起始点和结束点标记
    TMarker* startYZ = new TMarker(x_history[0][2], x_history[0][1], 29);
    startYZ->SetMarkerColor(kGreen);
    startYZ->SetMarkerSize(1.2);
    startYZ->Draw("same");

    TMarker* endYZ = new TMarker(x_history[actual_steps-1][2], 
                                x_history[actual_steps-1][1], 29);      
                               
    endYZ->SetMarkerColor(kRed);
    endYZ->SetMarkerSize(1.2);
    endYZ->Draw("same");

    c1->Update();
    
    
}

//<---------------------main function end----------------------->//

