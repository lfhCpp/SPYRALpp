// 运行方式（在 RK4/ 目录下）：
//   root -l -b -q test/testRK4.C
// 前置条件：先构建共享库 build/librk4.so（cmake --build build -j）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "sim/ParticleDynamics.h"
#include "../extern/rk4_solver/include/rk4_solver.hpp"



//---------------------------- static define-----------------------//


//--------------------->main function start<-----------------------//
void testRK4(){

     
    // 创建动力学对象
    ParticleDynamics D2;

    // 创建RK4积分器和事件对象
    double dt = 1e-10;       // 时间步长（秒）

    rk4_solver::Integrator<6, ParticleDynamics> 
        integrator(D2, &ParticleDynamics::equationMotion, dt);

    rk4_solver::Event<6, ParticleDynamics> 
        event(D2, &ParticleDynamics::StopEvent);

    // 初始条件:
    double x0[6];

    //D2.initialState(10, 0.0, 0.0, 0., 0., 0.0, x0);
    D2.initialState(10, 0.0, 0.0, 0.0,60.0 * (TMath::Pi()/ 180.0),0., x0);

    double t0 = 0.0; // 初始时间
    
    // 模拟参数
    const int max_steps = 100000; // 最大总步数
    double t_history[max_steps];
    double x_history[max_steps][6];
    
    // 运行轨迹模拟
    //rk4_solver::loop(integrator, t0, x0, t_history, x_history);

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
        graph3D->SetPoint(i, x, y, z);
        graphXY->SetPoint(i, x, y);  // XY平面投影
        graphXZ->SetPoint(i, z, x);  // XZ平面投影
        graphYZ->SetPoint(i, z, y);  // YZ平面投影
    }
    cout<<"Z:"<<x_history[actual_steps-1][2]<<" m"<<endl;

    TCanvas* c1 = new TCanvas("c1", "Particle Trajectory and Projections", 1200, 800);
    c1->Divide(2, 2); // 2行2列，共4个区域

    c1->cd(1); // 切换到第一个区域
    gPad->SetTheta(30); // 设置视角
    gPad->SetPhi(30);   // 设置视角

    graph3D->SetTitle("3D Particle Trajectory;X (m);Y (m);Z (m)");
    graph3D->SetMarkerColor(kBlue);
    graph3D->SetMarkerStyle(20);
    graph3D->SetMarkerSize(0.5);
    graph3D->SetLineColor(kRed);
    graph3D->SetLineWidth(2);
    graph3D->Draw("PLine");



    // ==================== 2. XY平面投影 ====================
    c1->cd(2); // 切换到第二个区域
    graphXY->SetTitle("XY Projection;X (m);Y (m)");
    graphXY->SetMarkerColor(kBlue);
    graphXY->SetMarkerStyle(20);
    graphXY->SetMarkerSize(0.5);
    graphXY->SetLineColor(kBlue);
    graphXY->SetLineWidth(1);
    graphXY->Draw("APL");

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
    graphXZ->SetTitle("XZ Projection;Z (m);X (m)");
    graphXZ->SetMarkerColor(kBlue);
    graphXZ->SetMarkerStyle(20);
    graphXZ->SetMarkerSize(0.5);
    graphXZ->SetLineColor(kBlue);
    graphXZ->SetLineWidth(1);
    graphXZ->Draw("APL");

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
    graphYZ->SetTitle("YZ Projection;Z (m);Y (m)");
    graphYZ->SetMarkerColor(kBlue);
    graphYZ->SetMarkerStyle(20);
    graphYZ->SetMarkerSize(0.5);
    graphYZ->SetLineColor(kBlue);
    graphYZ->SetLineWidth(1);
    graphYZ->Draw("APL");

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

