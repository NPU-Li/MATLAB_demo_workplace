%% GCRS到ITRS坐标转换（含速度）- 完整示例
% 使用CIO方法和IAU 2000/2006模型
% 输入时间：2004年4月6日 07:51:28.386 UTC
fMJD_UTC = 53101.3274118751;
% GCRS位置和速度（km, km/s）
r_gcrs = [5102.5089530; 6123.0113955; 6378.1369371];
v_gcrs = [-4.743220161; 0.790536492; 5.533755724];
% ITRF预期结果（用于验证）
r_itrf_expected = [-1033.4793830; 7901.2952754; 6380.3565958];
v_itrf_expected = [-3.225636520; -2.872451450; 5.531924446];
% 使用Vallado论文中的精确EOP参数
eopobj = USNO(); 
%eopobj = eopobj.initWithFinalsHttp(); 
% [xp,yp,du,dt] = eopobj.getEOP(fMJD_UTC);
xp = -0.140682;   % 极移X（角秒）
yp = 0.333309;    % 极移Y（角秒）
du = -0.439962;   % UT1-UTC（秒）
dt = 32;          % TAI-UTC（秒）
%% 计算GCRS到ITRS的转换矩阵
% 1. 计算3x3位置转换矩阵（GCRS -> ITRS）
GC2IT = IERS.GCRS2ITRS(fMJD_UTC, dt, du, xp, yp);

%% 扩展为6x6状态转换矩阵（包含速度）
% 对于GCRS到ITRS转换，角速度方向为正
%w = [0, 0, IERS.DS2R];  % GCRS -> ITRS方向
theta1 = IERS.ERA(fMJD_UTC, du);
theta2 = IERS.ERA(fMJD_UTC + 1/86400, du);  % 1秒后
omega_correct = theta2 - theta1;  % 弧度/秒
w = [0, 0, omega_correct];
% 扩展为9x9矩阵以包含加速度项
GC2IT_full = IERS.wrot(w, GC2IT);

% 缩减为6x6用于位置和速度
GC2IT_state = GC2IT_full(1:6, 1:6);

%% 执行转换
state_gcrs = [r_gcrs; v_gcrs];
state_itrf = GC2IT_state * state_gcrs;

r_itrf = state_itrf(1:3);
v_itrf = state_itrf(4:6);

%% 显示结果
fprintf('===== GCRS到ITRS坐标转换 =====\n\n');
fprintf('输入时间 (UTC): %.10f MJD\n', fMJD_UTC);
fprintf('对应: 2004-04-06 07:51:28.386 UTC\n\n');

fprintf('EOP参数:\n');
fprintf('  xp = %.6f 角秒\n', xp);
fprintf('  yp = %.6f 角秒\n', yp);
fprintf('  UT1-UTC = %.6f 秒\n', du);
fprintf('  TAI-UTC = %d 秒\n\n', dt);

fprintf('GCRS输入:\n');
fprintf('  位置: [%12.7f, %12.7f, %12.7f] km\n', r_gcrs);
fprintf('  速度: [%12.7f, %12.7f, %12.7f] km/s\n\n', v_gcrs);

fprintf('ITRS输出:\n');
fprintf('  位置: [%12.7f, %12.7f, %12.7f] km\n', r_itrf);
fprintf('  速度: [%12.7f, %12.7f, %12.7f] km/s\n\n', v_itrf);

%% 验证
fprintf('验证（与Vallado论文结果比较）:\n');
fprintf('  预期ITRS位置: [%12.7f, %12.7f, %12.7f] km\n', r_itrf_expected);
fprintf('  预期ITRS速度: [%12.7f, %12.7f, %12.7f] km/s\n\n', v_itrf_expected);

% 计算误差
pos_error_m = norm(r_itrf - r_itrf_expected) * 1000;  % 转换为米
vel_error_mmps = norm(v_itrf - v_itrf_expected) * 1e6;  % 转换为mm/s

fprintf('转换误差:\n');
fprintf('  位置误差: %.3f 米\n', pos_error_m);
fprintf('  速度误差: %.3f 毫米/秒\n\n', vel_error_mmps);

%% 性能指标
pos_rel_error = norm(r_itrf - r_itrf_expected) / norm(r_itrf_expected);
vel_rel_error = norm(v_itrf - v_itrf_expected) / norm(v_itrf_expected);

fprintf('相对误差:\n');
fprintf('  位置相对误差: %.2e\n', pos_rel_error);
fprintf('  速度相对误差: %.2e\n\n', vel_rel_error);

%% 转换矩阵信息
fprintf('转换矩阵信息:\n');
fprintf('  3x3位置转换矩阵行列式: %.12f\n', det(GC2IT));
fprintf('  3x3位置转换矩阵正交性误差: %.2e\n', norm(GC2IT*GC2IT' - eye(3), 'fro'));
fprintf('  地球自转角速度: %.10f rad/s\n', w);