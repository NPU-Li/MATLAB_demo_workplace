%% GCRS ↔ ITRS 双向坐标转换验证
% 使用CIO方法和IAU 2000/2006模型
% 基于矩阵求逆实现高精度双向转换
clc;
clear;
% 输入时间：2004年4月6日 07:51:28.386 UTC
fMJD_UTC = 53101.3274118751;

% GCRS位置和速度（km, km/s）- Vallado论文数据
r_gcrs = [5102.5089530; 6123.0113955; 6378.1369371];
v_gcrs = [-4.743220161; 0.790536492; 5.533755724];

% ITRS预期结果（用于验证）- Vallado论文数据
r_itrf_expected = [-1033.4793830; 7901.2952754; 6380.3565958];
v_itrf_expected = [-3.225636520; -2.872451450; 5.531924446];

% 使用Vallado论文中的精确EOP参数
% xp = -0.140682;   % 极移X（角秒）
% yp = 0.333309;    % 极移Y（角秒）
% du = -0.439962;   % UT1-UTC（秒）
% dt = 32;          % TAI-UTC（秒）

%使用服务器下载最新EOP参数
eopobj = USNO(); 
%eopobj = eopobj.initWithFinalsHttp(); 

eopobj = eopobj.initWIthEOPdata(); %采用之前下载好的文件
[xp,yp,du,dt] = eopobj.getEOP(fMJD_UTC);
% xp 
% yp
% du
% dt
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

%% 执行正向转换：GCRS -> ITRS
state_itrf = GC2IT_state * [r_gcrs; v_gcrs];
r_itrf = state_itrf(1:3);
v_itrf = state_itrf(4:6);

%% 执行反向转换：ITRS -> GCRS
%state_gcrs_back = IT2GC_state * state_itrf;
state_gcrs_back = IT2GC_state * [r_itrf_expected; v_itrf_expected];
r_gcrs_back = state_gcrs_back(1:3);
v_gcrs_back = state_gcrs_back(4:6);

%% 计算误差
% 正向转换误差
pos_error_fwd_mm = norm(r_itrf - r_itrf_expected) * 1e6;
vel_error_fwd_mmps = norm(v_itrf - v_itrf_expected) * 1e6;

% 反向转换误差
pos_error_back_mm = norm(r_gcrs_back - r_gcrs) * 1e6;
vel_error_back_mmps = norm(v_gcrs_back - v_gcrs) * 1e6;

%% 输出结果
fprintf('GCRS ↔ ITRS 双向转换验证结果\n');
fprintf('时间: 2004-04-06 07:51:28.386 UTC\n');
fprintf('角速度: %.12f rad/s\n\n', omega_ut1);

fprintf('正向转换 (GCRS -> ITRS):\n');
fprintf('  位置误差: %.3f 毫米\n', pos_error_fwd_mm);
fprintf('  速度误差: %.3f 毫米/秒\n\n', vel_error_fwd_mmps);

fprintf('反向转换 (ITRS -> GCRS):\n');
fprintf('  位置误差: %.3f 毫米\n', pos_error_back_mm);
fprintf('  速度误差: %.3f 毫米/秒\n\n', vel_error_back_mmps);

%% 双向一致性验证
% GCRS -> ITRS -> GCRS
state_gcrs_roundtrip = IT2GC_state * (GC2IT_state * [r_gcrs; v_gcrs]);
roundtrip_pos_error = norm(state_gcrs_roundtrip(1:3) - r_gcrs) * 1e6;
roundtrip_vel_error = norm(state_gcrs_roundtrip(4:6) - v_gcrs) * 1e6;

fprintf('双向一致性验证:\n');
fprintf('  GCRS->ITRS->GCRS:\n');
fprintf('    位置误差: %.3f 毫米\n', roundtrip_pos_error);
fprintf('    速度误差: %.3f 毫米/秒\n', roundtrip_vel_error);
