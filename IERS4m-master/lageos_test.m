clear;
clc;
mjd = utc2mjd(2019, 12, 29, 0, 0, 0);
% mjd = utc2mjd(2020, 1, 5, 0, 0, 0);
%current is 58846 2020 1 5 is 58853
fMJD_UTC = mjd;

%输入ITRS数据
% PL51  -6081.051297   4312.745884  -9813.822846 999999.999999
% VL51  52005.689664   8852.857209 -28503.945088 999999.999999
%data of 2019/12/29
r_itrf = [-6081.051297; 4312.745884; -9813.822846];
v_itrf = [5.2005689664; 0.8852857209; -2.8503945088];
% PL51  -8459.396554  -1494.796210   8701.915129 999999.999999
% VL51 -41126.340246  27177.508498 -35021.859535 999999.999999
%data of 2020/01/05
% r_itrf = [-8459.396554; -1494.796210; 8701.915129];
% v_itrf = [-4.1126340246; 2.7177508498; -3.5021859535];
%使用服务器下载最新EOP参数
eopobj = USNO(); 
eopobj = eopobj.initWIthEOPdata(); 
[xp,yp,du,dt] = eopobj.getEOP(fMJD_UTC);
%% 计算正确的角速度（基于UT1）
theta1 = IERS.ERA(fMJD_UTC, du);
theta2 = IERS.ERA(fMJD_UTC + 1/86400, du);
omega_ut1 = theta2 - theta1;  % 弧度/秒

%% 计算状态转换矩阵
% 位置转换矩阵（GCRS -> ITRS）
GC2IT = IERS.GCRS2ITRS(fMJD_UTC, dt, du, xp, yp);

% 构建正向状态转换矩阵（GCRS -> ITRS）
w_fwd = [0, 0, omega_ut1];
GC2IT_full = IERS.wrot(w_fwd, GC2IT);
GC2IT_state = GC2IT_full(1:6, 1:6);

% 反向状态转换矩阵（ITRS -> GCRS）
IT2GC_state = inv(GC2IT_state);

%% 执行反向转换：ITRS -> GCRS
%state_gcrs_back = IT2GC_state * state_itrf;
state_gcrs_back = IT2GC_state * [r_itrf; v_itrf];
format long g
r = state_gcrs_back(1:3)
v = state_gcrs_back(4:6)
[a,e,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(r,v);
%SMA a 半长轴 ECC c 偏心率 INC i 轨道倾角，这里是弧度制 RAAN 升交点赤经 O AOP 近地点幅角 o 以上均为弧度制
% mu 真近点角，需转换为平近点角度 M
E = 2 * atan2(sqrt(1 - e) * sin(nu / 2), sqrt(1 + e) * cos(nu / 2));

% 步骤2：通过开普勒方程计算平近点角 M (Mean Anomaly, 弧度)
M_rad = E - e * sin(E); % 经典开普勒方程

% 步骤3：将M从弧度转换为度
M_deg = rad2deg(M_rad); % 这就是轨道传播器需要的 M 参数

% ============== 2. 其他参数的单位转换 (弧度 -> 度) ==============
INC_deg = rad2deg(i);   % 轨道倾角
RAAN_deg = rad2deg(O);  % 升交点赤经
AOP_deg = rad2deg(o);   % 近地点幅角

% ============== 3. 打印转换结果，准备输入轨道传播器 ==============
fprintf('========== 轨道传播器输入参数 ==========\n');
fprintf('SMA  (半长轴 a): %.10f km\n', a);
fprintf('ECC  (偏心率 e): %.10f\n', e);
fprintf('INC  (轨道倾角): %.10f deg\n', INC_deg);
fprintf('RAAN (升交点赤经): %.10f deg\n', RAAN_deg);
fprintf('AOP  (近地点幅角): %.10f deg\n', AOP_deg);
fprintf('M    (平近点角): %.10f deg\n\n', M_deg);




