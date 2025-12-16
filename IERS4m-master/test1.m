%% GCRS ↔ ITRS 双向坐标转换验证（最终简化版）
% 使用CIO方法和IAU 2000/2006模型
% 基于矩阵求逆实现高精度双向转换

% 输入时间：2004年4月6日 07:51:28.386 UTC
fMJD_UTC = 58846;
eopobj = USNO(); 
eopobj = eopobj.initWithFinalsHttp(); 
[xp,yp,du,dt] = eopobj.getEOP(fMJD_UTC);
%% 计算状态转换矩阵
% 位置转换矩阵（GCRS -> ITRS）
GC2IT = IERS.GCRS2ITRS(fMJD_UTC, dt, du, xp, yp);

% 构建正向状态转换矩阵（GCRS -> ITRS）
w_fwd = [0, 0, omega_ut1];
GC2IT_full = IERS.wrot(w_fwd, GC2IT);
GC2IT_state = GC2IT_full(1:6, 1:6);

% 反向状态转换矩阵（ITRS -> GCRS）
IT2GC_state = inv(GC2IT_state);
%% 计算正确的角速度（基于UT1）
theta1 = IERS.ERA(fMJD_UTC, du);
theta2 = IERS.ERA(fMJD_UTC + 1/86400, du);
omega_ut1 = theta2 - theta1;  % 弧度/秒
%%

% ITRS预期结果（用于验证）- Vallado论文数据
r_itrf = [-6081.051297; 4312.745884; -9813.822846];
v_itrf = [5200.5689664; 885.2857209; -2850.3945088];

%% 执行反向转换：ITRS -> GCRS
%state_gcrs_back = IT2GC_state * state_itrf;
state_gcrs_back = IT2GC_state * [r_itrf; v_itrf];
r_gcrs_back = state_gcrs_back(1:3)
v_gcrs_back = state_gcrs_back(4:6)
