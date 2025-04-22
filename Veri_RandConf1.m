clear;
clc;
close all;
%% 根据优化区间绘制末端运动图像
% 与优化构型具有相同的末端位姿，但各关节角不在优化范围内的操作构型，控制末端运动
% 绘制末端位置误差和方位角误差图像
coef_joiwei = [1.40 1.40 1.40 1.05 0.70 0.70 0.35];
% Confi_rand2 = [143.6971574114474; 24.00; -69.3535260543474; -44.5615875940660; 97.9820912575147; -132.9894691309734; 122.5620376277255];
% Confi_rand2 = [138.5051485692268 24.50 -64.6051678816212 -19.6906061872884 -105.0342345610541 133.9433100102033 -54.7221341795219];
Confi_rand2 = [120.0924773068039 167 36.5343861490611 35.3315781098094 97.4317608097554 -83.0799536585081 -53.1704697041120];

conf_rad = Confi_rand2 * pi / 180;
T_sta = xbForKinEnd(conf_rad);      
R_sta = T_sta(1:3,1:3);            % 末端姿态旋转角
rpy_sta = rotm2rpy(R_sta);        % 末端姿态方位角的理论值

%% 五次多项式求解
T_all = 15;      % 运动总时间
t_int = 0.002;  % 采样间隔
xyz_sta = T_sta(1:3,4);     % 运动起点，需要根据具体构型改变
d_dul = 0.3;    % 运动距离，注意区别正反向：正向运动时距离为正值，负向运动时距离为负值
a3 = 10/(T_all^3);
a4 = -15/(T_all^4);
a5 = 6/(T_all^5);

t_sam = 0:t_int:T_all;
s_sca = a3 * t_sam.^3 + a4 * t_sam.^4 + a5 * t_sam.^5;
v_sca = 3*a3 * t_sam.^2 + 4*a4 * t_sam.^3 + 5*a5 * t_sam.^4;
a_sca = 6*a3 * t_sam + 12*a4 * t_sam.^2 + 20*a5 * t_sam.^3;

s_eef = xyz_sta(3) + d_dul * s_sca;
% figure('name','The relationship between displacement and time');
% plot(t_sam,s_xpos);
% xlabel('t/(s)');
% ylabel('Displacement/(m)');
v_eef = d_dul * v_sca;  % 根据位移d_dul决定运动方向
% figure('name','The relationship between velocity and time');
% plot(t_sam,v_pos);
% xlabel('t/(s)');
% ylabel('Velocity/(m/s)');
a_eef = d_dul * a_sca;
% figure('name','The relationship between acceleration and time');
% plot(t_sam,a_xpos);
% xlabel('t/(s)');
% ylabel('Acceleration of EEF/(m/s^2)');
%% 定义机械臂末端速度
mm = size(v_eef,2);
v_end = zeros(6,mm);
v_end(3,:) = v_eef;      % 定义笛卡尔空间的末端速度
v_zposi = v_end;        % z轴正向运动
v_ynega = zeros(6,mm);
v_ynega(2,:) = -v_eef;
%% 求解关节角速度
theta_Orientation = conf_rad ;       % 定义机械臂初始构型    % ini_confi*pi/180      % 单位转化，将°转化为rad

angle_alljoint = zeros(7,mm);      % 记录各关节角的中间变量，用于绘制角度变化曲线
vel_theta = zeros(7,mm);         % 各关节角速度，用于绘制角速度变化曲线
v_allend = zeros(6,mm);          % 根据得到的角速度计算末端速度，用于验证误差
ang_dedai = zeros(7,mm);         % 记录每个采样周期内的角度变化量
angle_alljoint(:,1) = theta_Orientation;          

for ii=1:mm
    Theta_update = angle_alljoint(:,ii);
    matrix_jacobi = JoiAngToJacMatBase(Theta_update);
    matrix_zhuanzhi = matrix_jacobi.';
    matrix_ladder = matrix_jacobi * matrix_zhuanzhi;      % inv(matrix_jacobi*matrix_zhuanzhi);
    jacobi_inv = matrix_zhuanzhi / matrix_ladder;
    v_joint = jacobi_inv * v_ynega(:,ii);      %计算关节速度v_joint
    vel_theta(:,ii) = v_joint;                  % 记录关节角速度，绘制图像
    v_allend(:,ii) = matrix_jacobi * v_joint;    % 根据求出的关节速度计算末端速度进行验证
    
    ang_dedai(:,ii) = vel_theta(:,ii) * t_int;      % 每个采样周期内的角度变化量
    angle_alljoint(:,ii+1) = angle_alljoint(:,ii) + ang_dedai(:,ii);
end
%% 绘制关节角度和角速度变化曲线
angle_alljoint(:,mm+1) = [];         % 去除关节角度中多余的一列
JoiAng_deg = angle_alljoint*180/pi;
%%
figure('position',[100 100 640 540],'name','Joint angle');
plot(t_sam,JoiAng_deg(1,:),'-*','Color','#0072BD','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
hold on;
plot(t_sam,JoiAng_deg(2,:),'-^','Color','r','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiAng_deg(3,:),'-o','Color','#EDB120','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiAng_deg(4,:),'-x','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiAng_deg(5,:),'-v','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiAng_deg(6,:),'-d','Color','#4DBEEE','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiAng_deg(7,:),'-s','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);       % 'Marker','s',
set(gca,'Position',[0.13 0.15 0.80 0.74],'XColor',[1.00,0.00,1.00],'YColor',[1.00,0.00,1.00]);
% figure_FontSize = 8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',16);
hold off;
xlabel('Time/(s)','VerticalAlignment','middle');
ylabel('Joint angle/(\circ)','VerticalAlignment','baseline');      % \circ
lgd = legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','NumColumns',4);
% lgd.Location = 'best';
lgd.Position = [0.140468747600913,0.523278381515723,0.759062504798174,0.06814814959632];
lgd.Orientation = 'horizontal';
lgd.NumColumns = 7;
% legend('boxoff');
%%
vel_limit = zeros(mm,1);
vel_limit(:,1) = 2.2;
JoiVel_deg = vel_theta*180/pi;
Min_JoiVel = max(JoiVel_deg(5,:),[],2)
Max_JoiVel = min(JoiVel_deg(1:4,:),[],2)
figure('position',[500 200 640 540],'name','Joint angular velocities');
% figure('name','Joint angular velocities');
plot(t_sam,JoiVel_deg(1,:),'-*','Color','#0072BD','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
hold on;
plot(t_sam,JoiVel_deg(2,:),'-^','Color','r','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiVel_deg(3,:),'-o','Color','#EDB120','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiVel_deg(4,:),'-x','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiVel_deg(5,:),'-v','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiVel_deg(6,:),'-d','Color','#4DBEEE','MarkerIndices',1:800:length(t_sam),'MarkerSize',12);
plot(t_sam,JoiVel_deg(7,:),'-s','MarkerIndices',1:500:length(t_sam),'MarkerSize',12);
plot(t_sam,vel_limit,'--','Color','k');
plot(t_sam,-vel_limit,'--','Color','k');
% yticks(-4:2:4);
set(gca,'Position',[0.1 0.12 0.80 0.74],'XColor',[1.00,0.00,1.00],'YColor',[1.00,0.00,1.00]);
% figure_FontSize = 8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','baseline');
% set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',16);
hold off;
xlabel('Time/(s)');
ylabel('Joint angular velocities/(\circ/s)');
lgd_2 = legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7');
% lgd_2.Location = 'best';
lgd_2.Position = [0.114218747600914,0.127829431081607,0.759062504798174,0.06814814959632];
lgd_2.Orientation = 'horizontal';
lgd_2.NumColumns = 7;
%% 求解结果正运动学验证
T_end = xbForKinEnd(angle_alljoint(:,mm));       % *180/pi
%% 运动过程误差求解
% 根据求解得到的机械臂构型(angle_alljoint)，求解正运动学
% 比较末端位置与规划位置的误差，和方位角误差

% 规划速度和位置
xyz_plan = zeros(3,mm);
xyz_plan(1,:) = xyz_sta(1);
xyz_plan(2,:) = xyz_sta(2);
xyz_plan(3,:) = s_eef;
rpy_plan = zeros(3,mm);
rpy_plan(1,:) = rpy_sta(1);
rpy_plan(2,:) = rpy_sta(2);
rpy_plan(3,:) = rpy_sta(3);

eef_EulerAng = zeros(3,mm);
eef_posi = zeros(3,mm);
for ii = 1:mm
    T_weizi = xbForKinEnd(angle_alljoint(:,ii));
    eef_posi(:,ii) = T_weizi(1:3,4);
    R_zitai = T_weizi(1:3,1:3);
    eef_EulerAng(:,ii) = rotm2rpy(R_zitai);
end
% % 绘制末端位置变化曲线
% figure('name','the relationship between end-effector position and time');
% plot(t_sam,xyz_plan(1,:),'-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,xyz_plan(2,:),'-+','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,xyz_plan(3,:),'-x','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_posi(1,:),'-o','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_posi(2,:),'-.','Color','k','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_posi(3,:),'--','Color','r');
% hold off;
% xlabel('t/(s)');
% ylabel('End-Effector Position/(m)');
% legend('X-Axis planning','Y-Axis planning','Z-Axis planning','X-Axis Fact','Y-Axis Fact','Z-Axis Fact');
% % 绘制末端方位角
% figure('name','the relationship between Euler angle and time');
% plot(t_sam,rpy_plan(1,:),'-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,rpy_plan(2,:),'-x','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,rpy_plan(3,:),'-+','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_EulerAng(1,:),'-o','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_EulerAng(2,:),'-.','Color','k','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,eef_EulerAng(3,:),'--','Color','r');
% hold off;
% xlabel('t/(s)');
% ylabel('End-Effector Orientation/(rad)');
% legend('r planning','p planning','y planning','r Fact','p Fact','y Fact');

%% 位置和方位角误差计算及图像绘制
% 误差计算：规划值-实际值
err_ori = rpy_plan - eef_EulerAng;
err_posi = (xyz_plan - eef_posi) * 1000;        % 将误差单位从m转化为mm
% 误差图像绘制
% figure('name','End-Effector Position Error');
% plot(t_sam,err_posi(1,:),'-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,err_posi(2,:),'-+','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_posi(3,:),'-x','MarkerIndices',1:300:length(t_sam));
% hold off;
% xlabel('t/(s)');
% ylabel('End-Effector Position Error/(mm)');
% legend('X-Axis error','Y-Axis error','Z-Axis error');
% figure('name','End-Effector Orientation Error');
% plot(t_sam,err_ori(1,:),'-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,err_ori(2,:),'-+','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_ori(3,:),'-x','MarkerIndices',1:300:length(t_sam));
% hold off;
% xlabel('t/(s)');
% ylabel('End-Effector Orientation Error/(rad)');
% legend('r error','p error','y error');

% %% 绘制双y轴图像
% % 末端位置和误差
% figure('name','End-Effector Position and Error');
% yyaxis left;
% set(gca,'YColor','k');          % 定义y轴颜色
% plot(t_sam,xyz_plan(1,:),'r-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,xyz_plan(2,:),'g-x','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,xyz_plan(3,:),'b-+','MarkerIndices',1:300:length(t_sam));
% xlabel('t/(s)');
% ylabel('End-Effector Position/(m)','Color','k');
% yyaxis right;
% set(gca,'YColor','k');
% plot(t_sam,err_posi(1,:),'k-.','Marker','s','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_posi(2,:),'m--','Marker','d','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_posi(3,:),'-o','MarkerIndices',1:300:length(t_sam));
% hold off;
% ylabel('End-Effector Position Error/(mm)','Color','k');
% legend('X-Axis Position','Y-Axis Position','Z-Axis Position','X-Axis Error','Y-Axis Error','Z-Axis Error');
% % 末端方位角和误差
% figure('name','End-Effector Orientation and Error');
% yyaxis left;
% set(gca,'YColor','k');
% plot(t_sam,rpy_plan(1,:),'r-*','MarkerIndices',1:300:length(t_sam));
% hold on;
% plot(t_sam,rpy_plan(2,:),'g-x','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,rpy_plan(3,:),'b-+','MarkerIndices',1:300:length(t_sam));
% xlabel('t/(s)');
% ylabel('End-Effector Orientation/(rad)','Color','k');
% yyaxis right;
% set(gca,'YColor','k');
% plot(t_sam,err_ori(1,:),'k-.','Marker','s','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_ori(2,:),'m--','Marker','d','MarkerIndices',1:300:length(t_sam));
% plot(t_sam,err_ori(3,:),'-o','MarkerIndices',1:300:length(t_sam));
% hold off;
% ylabel('End-Effector Orientation Error/(rad)','Color','k');
% legend('R Orientation','P Orientation','Y Orientation','R Error','P Error','Y Error');
% % legend('boxoff');
%% 绘制模拟运动图像
Poi_Zero = zeros(3,mm);
Poi_1 = zeros(3,mm);
Poi_2 = zeros(3,mm);
Poi_3 = zeros(3,mm);
Poi_34 = zeros(3,mm);
Poi_4 = zeros(3,mm);
Poi_45 = zeros(3,mm);
Poi_5 = zeros(3,mm);
Poi_6 = zeros(3,mm);
Poi_7 = zeros(3,mm);

for ii = 1 : mm
    Medi_JoiAng = angle_alljoint(:,ii);
    [Poi_1(:,ii), Poi_2(:,ii), Poi_3(:,ii), Poi_34(:,ii), Poi_4(:,ii), Poi_45(:,ii), Poi_5(:,ii), Poi_6(:,ii), Poi_7(:,ii)] = ForKinexb(Medi_JoiAng); 
end
kk = 1;

for ii = 1 :500: mm
    Posi_X(kk,:) = [Poi_Zero(1,ii),Poi_1(1,ii),Poi_2(1,ii),Poi_3(1,ii),...
             Poi_34(1,ii),Poi_4(1,ii),Poi_45(1,ii),...
            Poi_5(1,ii),Poi_6(1,ii),Poi_7(1,ii)];
    Posi_Y(kk,:) = [Poi_Zero(2,ii),Poi_1(2,ii),Poi_2(2,ii),Poi_3(2,ii),...
               Poi_34(2,ii),Poi_4(2,ii),Poi_45(2,ii),...
             Poi_5(2,ii),Poi_6(2,ii),Poi_7(2,ii)];
    Posi_Z(kk,:) = [Poi_Zero(3,ii),Poi_1(3,ii),Poi_2(3,ii),Poi_3(3,ii),...
                 Poi_34(3,ii),Poi_4(3,ii),Poi_45(3,ii),...
                 Poi_5(3,ii),Poi_6(3,ii),Poi_7(3,ii)];
    kk = kk + 1;
end

figure;
nn = size(Posi_X,1);
% figure('position',[500 300 260 220],'name','Joint angular velocities');

for ii = 1 : nn
    plot3(Posi_X(ii,1:2), Posi_Y(ii,1:2), Posi_Z(ii,1:2),'Color','r','LineWidth',1.0);
    plot3(Posi_X(ii,2:3), Posi_Y(ii,2:3), Posi_Z(ii,2:3),'Color','g','LineWidth',1.0);
    plot3(Posi_X(ii,3:4), Posi_Y(ii,3:4), Posi_Z(ii,3:4),'Color','c','LineWidth',1.0);
    plot3(Posi_X(ii,4:5), Posi_Y(ii,4:5), Posi_Z(ii,4:5),'Color','m','LineWidth',1.0);
    plot3(Posi_X(ii,5:6), Posi_Y(ii,5:6), Posi_Z(ii,5:6),'Color','m','LineWidth',1.0);

    plot3(Posi_X(ii,6:7), Posi_Y(ii,6:7), Posi_Z(ii,6:7),'Color',[0.72,0.27,1.00],'LineWidth',1.0);
    plot3(Posi_X(ii,7:8), Posi_Y(ii,7:8), Posi_Z(ii,7:8),'Color',[0.72,0.27,1.00],'LineWidth',1.0);     % '#A2142F'

    plot3(Posi_X(ii,8:9), Posi_Y(ii,8:9), Posi_Z(ii,8:9),'Color','#4DBEEE','LineWidth',1.0);

    plot3(Posi_X(ii,9:10), Posi_Y(ii,9:10), Posi_Z(ii,9:10),'Color','#7E2F8E','LineWidth',1.0);
    hold on;
end
plot3(Posi_X(:,10),Posi_Y(:,10),Posi_Z(:,10),'Color','b','LineWidth',1.0);
xlim([-1.11 0.2]);
% ylim([])
zlim([-2 0.4]);
xx_1 = [0.2 0.2 0.2 -1.11];
yy_1 = [-5 -5 0 0 ];
zz_1 = [-2 0.4 0.4 0.4];
plot3(xx_1,yy_1,zz_1,'-','Color',[1.00,0.00,1.00]);
set(gca,'TickDir','in','XColor',[1.00,0.00,1.00],'YColor',[1.00,0.00,1.00],'ZColor',[1.00,0.00,1.00]);
xlabel('X/(m)','Position',[-0.543417401885293,-5.416440107093514,-2.269274350766224]);  % 'HorizontalAlignment','center',"VerticalAlignment","bottom","FontSize",12);
ylabel('Y/(m)',"VerticalAlignment","middle","HorizontalAlignment","right","FontSize",12);
zlabel('Z/(m)',"VerticalAlignment","middle","FontSize",12);
Arrow_1 = annotation('arrow',[0.726071428571424,0.671428571428571],[0.457571428571432,0.459047619047619],'HeadWidth',6,"HeadStyle","plain");
Str_1 = 'Random configuration';
Arrow_2 = annotation('textarrow',[0.431428571428571,0.486071428571429],[0.318095238095238,0.511904761904762],...
    'Interpreter','latex','String',Str_1,'FontSize',14);
Arrow_2.Position = [0.431428571428571,0.318095238095238,0.054642857142858,0.193809523809524];
Arrow_2.LineWidth = 1;
Arrow_2.HeadStyle = 'plain';
Arrow_2.HeadWidth = 6;
grid on;
view(-72,18);
%% 程序所需函数
function [Jaco_base] = JoiAngToJacMatBase(Theta)
%% 计算以基坐标系为参考的雅可比矩阵
%   输入量为包含7个关节角度的数组，单位为弧度rad；
%   输出量为使用旋量法以基坐标系为参考的雅可比矩阵
%% 机械臂关节常量定义



%% 关节角度赋值
% Theta_rad = Theta * pi / 180;
Theta_rad = Theta;

theta1 = Theta_rad(1);
theta2 = Theta_rad(2);
theta3 = Theta_rad(3);
theta4 = Theta_rad(4);
theta5 = Theta_rad(5);
theta6 = Theta_rad(6);
theta7 = Theta_rad(7);
%% 求解公式
JB_S1 = [-1;
          0;
          0;
          0;
          0;
          0];
JB_S2 = [0;
         -sin(theta1);
         -cos(theta1);
         0;
         -a0N*cos(theta1);
         a0N*sin(theta1)];
JB_S3 = [-cos(theta2);
         cos(theta1)*sin(theta2);
         -sin(theta1)*sin(theta2);
         a1N*sin(theta2);
         a1N*cos(theta1)*cos(theta2) - a0N*cos(theta2)*sin(theta1)*sin(theta2) + a0N*sin(theta1)*sin(theta2)*(cos(theta2) - 1);
         a0N*cos(theta1)*sin(theta2)*(cos(theta2) - 1) - a0N*cos(theta1)*cos(theta2)*sin(theta2) - a1N*cos(theta2)*sin(theta1)];
JB_S4 = [-cos(theta2);
         cos(theta1)*sin(theta2);
         -sin(theta1)*sin(theta2);
         cos(theta3)*sin(theta2)*(a1N + a3N) - cos(theta1)*sin(theta2)*(a0N*sin(theta1)*sin(theta2) + a1N*cos(theta1)*(cos(theta3) - 1) - a1N*cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta1)*sin(theta2)*(a1N*sin(theta1)*(cos(theta3) - 1) - a0N*cos(theta1)*sin(theta2) + a1N*cos(theta1)*cos(theta2)*sin(theta3));
         sin(theta1)*sin(theta2)*(a0N*(cos(theta2) - 1) + a1N*sin(theta2)*sin(theta3)) - (a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta2)*(a0N*sin(theta1)*sin(theta2) + a1N*cos(theta1)*(cos(theta3) - 1) - a1N*cos(theta2)*sin(theta1)*sin(theta3));
         cos(theta2)*(a1N*sin(theta1)*(cos(theta3) - 1) - a0N*cos(theta1)*sin(theta2) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) - (a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + cos(theta1)*sin(theta2)*(a0N*(cos(theta2) - 1) + a1N*sin(theta2)*sin(theta3))];
JB_S5 = [-cos(theta2);
         cos(theta1)*sin(theta2);
         -sin(theta1)*sin(theta2);
         - (sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2))*(a1N + a3N + a5N) - cos(theta1)*sin(theta2)*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - a1N*cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta1)*sin(theta2)*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + a1N*cos(theta1)*cos(theta2)*sin(theta3));
         sin(theta1)*sin(theta2)*(a0N*(cos(theta2) - 1) + a1N*sin(theta2)*sin(theta3) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1)) - (cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))*(a1N + a3N + a5N) - cos(theta2)*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - a1N*cos(theta2)*sin(theta1)*sin(theta3));
         cos(theta2)*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) - (cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))*(a1N + a3N + a5N) + cos(theta1)*sin(theta2)*(a0N*(cos(theta2) - 1) + a1N*sin(theta2)*sin(theta3) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1))];
JB_S6 = [sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3));
         sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) - cos(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)));
         sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) - cos(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)));
         (cos(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))*(a1N + a3N + a5N) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))*(a1N + a3N + a5N) - a1N*cos(theta2)*sin(theta1)*sin(theta3)) - (cos(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1))) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))*(a1N + a3N + a5N) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))*(a1N + a3N + a5N) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) + (cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)))*(a0N + a2N + a4N + a6N);
         (cos(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1))) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))*(a0N*(cos(theta2) - 1) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2))*(a1N + a3N + a5N) + (cos(theta5) - 1)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))*(a1N + a3N + a5N) + a1N*sin(theta2)*sin(theta3) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1)) + (cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))*(a0N + a2N + a4N + a6N) - (cos(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))*(a1N + a3N + a5N) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))*(a1N + a3N + a5N) - a1N*cos(theta2)*sin(theta1)*sin(theta3));
         (cos(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))*(a1N + a3N + a5N) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))*(a1N + a3N + a5N) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) - (cos(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))*(a0N*(cos(theta2) - 1) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2))*(a1N + a3N + a5N) + (cos(theta5) - 1)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))*(a1N + a3N + a5N) + a1N*sin(theta2)*sin(theta3) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1)) + (cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1))))*(a0N + a2N + a4N + a6N)];
JB_S7 = [- cos(theta2)*cos(theta6) - sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)));
         cos(theta1)*cos(theta6)*sin(theta2) - sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))));
         - sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))) - cos(theta6)*sin(theta1)*sin(theta2);
         (cos(theta2)*sin(theta6) - cos(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))))*(a1N + a3N + a5N + a7N) - (sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))) + cos(theta6)*sin(theta1)*sin(theta2))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))*(a1N + a3N + a5N) + sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))*(a0N + a2N + a4N + a6N) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))*(a1N + a3N + a5N) - cos(theta1)*sin(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) + (sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))) - cos(theta1)*cos(theta6)*sin(theta2))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))*(a1N + a3N + a5N) + sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1))))*(a0N + a2N + a4N + a6N) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))*(a1N + a3N + a5N) + sin(theta1)*sin(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) - a1N*cos(theta2)*sin(theta1)*sin(theta3));
         (sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))) + cos(theta6)*sin(theta1)*sin(theta2))*(a0N*(cos(theta2) - 1) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2))*(a1N + a3N + a5N) + cos(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) + (cos(theta5) - 1)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))*(a1N + a3N + a5N) + a1N*sin(theta2)*sin(theta3) + sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)))*(a0N + a2N + a4N + a6N) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1)) - (cos(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))) + cos(theta1)*sin(theta2)*sin(theta6))*(a1N + a3N + a5N + a7N) - (cos(theta2)*cos(theta6) + sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) + a0N*sin(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))*(a1N + a3N + a5N) + sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1))))*(a0N + a2N + a4N + a6N) + a1N*cos(theta1)*(cos(theta3) - 1) - sin(theta4)*(a1N + a3N)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))*(a1N + a3N + a5N) + sin(theta1)*sin(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) - a1N*cos(theta2)*sin(theta1)*sin(theta3));
         (cos(theta2)*cos(theta6) + sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))))*((a1N + a3N)*(cos(theta4) - 1)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - a0N*cos(theta1)*sin(theta2) - sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))*(a1N + a3N + a5N) + sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))*(a0N + a2N + a4N + a6N) - sin(theta4)*(a1N + a3N)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + a1N*sin(theta1)*(cos(theta3) - 1) + (cos(theta5) - 1)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))*(a1N + a3N + a5N) - cos(theta1)*sin(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) + a1N*cos(theta1)*cos(theta2)*sin(theta3)) - (cos(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)) - sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))) - sin(theta1)*sin(theta2)*sin(theta6))*(a1N + a3N + a5N + a7N) - (sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))) + sin(theta5)*(cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)) - sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))) - cos(theta1)*cos(theta6)*sin(theta2))*(a0N*(cos(theta2) - 1) - sin(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2))*(a1N + a3N + a5N) + cos(theta2)*(cos(theta6) - 1)*(a0N + a2N + a4N + a6N) + (cos(theta5) - 1)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3))*(a1N + a3N + a5N) + a1N*sin(theta2)*sin(theta3) + sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta3)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta5)*(cos(theta3)*sin(theta2)*sin(theta4) + cos(theta4)*sin(theta2)*sin(theta3)))*(a0N + a2N + a4N + a6N) + cos(theta3)*sin(theta2)*sin(theta4)*(a1N + a3N) + sin(theta2)*sin(theta3)*(a1N + a3N)*(cos(theta4) - 1))];
Jaco_lad = [JB_S1  JB_S2  JB_S3   JB_S4   JB_S5   JB_S6   JB_S7];
Jaco_base = [Jaco_lad(4:6,:);Jaco_lad(1:3,:)];
end

function [Tend_end] = xbForKinEnd(Theta)
%% 建立以末端坐标系为参考的小臂正运动学求解方程；
% Theta为包含七个关节角度的数组，单位是°；
% 求解Tend_1的关节螺旋轴和点的位置均为相对末端坐标系得到；
% 输出结果Tend_end为各关节旋转任意角度后的末端坐标系相对基坐标系的位姿矩阵；
%% 用于单位转化，不能同时运行
Theta_deg = Theta * 180 / pi;
% Theta_deg = Theta;
%% 

theta1 = Theta_deg(1);
theta2 = Theta_deg(2);
theta3 = Theta_deg(3);
theta4 = Theta_deg(4);
theta5 = Theta_deg(5);
theta6 = Theta_deg(6);
theta7 = Theta_deg(7);

M_ori = [-1     0   0   -a0N-a2N-a4N-a6N-a8N;
         0      1   0   0;
         0      0   -1  -a1N-a3N-a5N-a7N;
         0      0   0   1];
Tend_1 = [1     0               0               0;
          0     cosd(theta1)    -sind(theta1)   -(a1N + a3N + a5N + a7N)*sind(theta1);
          0     sind(theta1)    cosd(theta1)    (cosd(theta1) - 1)*(a1N + a3N + a5N + a7N);
          0     0               0               1];
Tend_2 = [cosd(theta2)       -sind(theta2)    0   (cosd(theta2) - 1)*(a2N + a4N + a6N + a8N);
          sind(theta2)       cosd(theta2)     0   sind(theta2)*(a2N + a4N + a6N + a8N);
          0                 0               1   0;
          0                 0               0   1];
Tend_3 = [1     0               0               0;
          0     cosd(theta3)     -sind(theta3)    -sind(theta3)*(a3N + a5N + a7N);
          0     sind(theta3)     cosd(theta3)     (cosd(theta3) - 1)*(a3N + a5N + a7N);
          0     0               0               1];
Tend_4 = [1     0               0               0;
          0     cosd(theta4)     -sind(theta4)    -sind(theta4)*(a5N + a7N);
          0     sind(theta4)     cosd(theta4)     (a5N + a7N)*(cosd(theta4) - 1);
          0     0               0               1];
Tend_5 = [1     0               0               0;
          0     cosd(theta5)     -sind(theta5)    -a7N*sind(theta5);
          0     sind(theta5)     cosd(theta5)     a7N*(cosd(theta5) - 1);
          0     0               0               1];
Tend_6 = [cosd(theta6)       -sind(theta6)        0   a8N*(cosd(theta6) - 1);
          sind(theta6)       cosd(theta6)         0   a8N*sind(theta6);
          0                 0                   1   0;
          0                 0                   0   1];
Tend_7 = [1     0               0               0;
          0     cosd(theta7)     -sind(theta7)    0;
          0     sind(theta7)     cosd(theta7)     0;
          0     0               0               1];

Tend_end = M_ori * Tend_1 * Tend_2 * Tend_3 * Tend_4 * Tend_5 * Tend_6 * Tend_7;
end

function rpy = rotm2rpy( R )
%% 位姿矩阵转化为欧拉角
% 输入R为3*3的位姿矩阵；
% 输出rpy为欧拉角
if abs(R(3 ,1) - 1.0) < 1.0e-12   % singularity
    a = 0.0;
    b = -pi / 2.0;
    c = atan2(-R(1, 2), -R(1, 3));
elseif abs(R(3, 1) + 1.0) < 1.0e-12   % singularity
    a = 0.0;
    b = pi / 2.0;
    c = -atan2(R(1, 2), R(1, 3));
else
    a = atan2(R(3, 2), R(3, 3));
    c = atan2(R(2, 1), R(1, 1));
    %     a = atan2(-R(3, 2), -R(3, 3));  %a另一个解
    %     c = atan2(-R(2, 1), -R(1, 1));  %c另一个解
    cosC = cos(c);
    sinC = sin(c);
    
    if abs(cosC) > abs(sinC)
        b = atan2(-R(3, 1), R(1, 1) / cosC);
    else
        b = atan2(-R(3, 1), R(2, 1) / sinC);
    end
end
rpy = [a; b; c];
end

function [Posi_TDB1, Posi_TDB2, Posi_TDB3, Posi_TDB334, Posi_TDB4, Posi_TDB445, Posi_TDB5, Posi_TDB6, Posi_TDB] = ForKinexb(Theta)
% 求解机械臂各连杆的端点位置坐标

 
 
Theta1=Theta(1);
Theta2=Theta(2);
Theta3=Theta(3);
Theta4=Theta(4);
Theta5=Theta(5);
Theta6=Theta(6);
Theta7=Theta(7);

TDB0=[0 1 0 0;
    -1 0 0 0;
    0 0 1 0;
    0 0 0 1];

TD01=[cos(Theta1) -sin(Theta1) 0 0;
        0           0        -1 -A0N;
      sin(Theta1) cos(Theta1) 0 0;
        0           0          0 1];
  
TD12=[cos(Theta2) -sin(Theta2) 0  0;
      0 0 -1 -A1N;
 sin(Theta2) cos(Theta2) 0 0;
 0 0 0 1];

TD23=[sin(Theta3) cos(Theta3) 0 0;
     0 0 1 A2N;
 cos(Theta3) -sin(Theta3) 0 0;
 0 0 0 1];

TD3_34 = [1     0       0       A3N;
          0     1       0       0;
          0     0       1       0;
          0     0       0       1];

TD34=[cos(Theta4) -sin(Theta4) 0 A3N;
sin(Theta4) cos(Theta4) 0 0;    
0 0 1 A4N;    
0 0 0 1];

TD4_45 = [1     0       0       A5N;
          0     1       0       0;
          0     0       1       0;
          0     0       0       1];

TD45=[-sin(Theta5) -cos(Theta5) 0 A5N;
cos(Theta5) -sin(Theta5) 0 0;    
0 0 1 A6N;    
0 0 0 1];

TD56=[cos(Theta6) -sin(Theta6) 0 0;  
    0 0 -1 -A7N; 
    sin(Theta6) cos(Theta6) 0 0;
      0 0 0 1];
  
TD67=[cos(Theta7) -sin(Theta7) 0 0;  
    0 0 1 A8N; 
    -sin(Theta7) -cos(Theta7) 0 0;  
0 0 0 1];

TD7E=[0 -1 0 0;
    0 0 -1 0;
    1 0 0 0;
    0 0 0 1];

TDB1 = TDB0*TD01;
Posi_TDB1 = TDB1(1:3, 4);       % 坐标系{1}的位置坐标
TDB2 = TDB1*TD12;
Posi_TDB2 = TDB2(1:3, 4);       % 坐标系{2}的位置坐标
TDB3 = TDB2*TD23;
Posi_TDB3 = TDB3(1:3, 4);       % 坐标系{3}的位置坐标
TDB3_34 = TDB3 * TD3_34;
Posi_TDB334 = TDB3_34(1:3, 4);  % 坐标系{3}{4}交点的位置坐标
TDB4 = TDB3*TD34;
Posi_TDB4 = TDB4(1:3, 4);       % 坐标系{4}的位置坐标
TDB4_45 =  TDB4 * TD4_45;
Posi_TDB445 = TDB4_45(1:3, 4);  % 坐标系{4}{5}交点的位置坐标
TDB5 = TDB4*TD45;
Posi_TDB5 = TDB5(1:3, 4);       % 坐标系{5}的位置坐标
TDB6 = TDB5*TD56;
Posi_TDB6 = TDB6(1:3, 4);       % 坐标系{6}的位置坐标
TDB = TDB6*TD67*TD7E;
Posi_TDB = TDB(1:3, 4);         % 末端执行器的{e}的位置坐标，与坐标系{7}的位置坐标相同

end
