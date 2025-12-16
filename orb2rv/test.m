%% 添加 THALASSA MEX 接口路径（改成你的实际安装目录）
addpath('/home/user/thalassa-1.4/thalassa-1.4/lib');

%% 地球引力常数 [km^3/s^2]
muEarth = 398600.4418;

%% 六颗卫星在 MJD 58834 的 Kepler 元素
% [a(km), e, i(deg), Omega(deg), omega(deg), M(deg)]
elements = [...
    7577.6, 0.001,  87.9,   3.0,   45.2, 126.8;  % Target
    9209.8, 0.218, 106.3, 347.2,   56.7, 255.5;  % Field1
    8486.6, 0.189,  34.7,  48.7,  224.3, 188.5;  % Field2
   10434.8, 0.313, 171.4, 327.5,  256.3, 291.2;  % Field3
    7988.8, 0.143,  57.1,  38.4,   59.9, 142.7;  % Field4
    7958.9, 0.047, 132.2, 313.4,  316.1, 188.3]; % Field5

%% 将轨道根数转换成笛卡尔初始状态 ECI
states0 = zeros(6,6);
for k = 1:6
    a      = elements(k,1);
    e      = elements(k,2);
    i_rad  = deg2rad(elements(k,3));
    Om_rad = deg2rad(elements(k,4));
    w_rad  = deg2rad(elements(k,5));
    M_rad  = deg2rad(elements(k,6));

    % 半通径 p
    p = a * (1 - e^2);

    % M -> E (Newton-Raphson)
    E_rad = M_rad;
    for iter = 1:50
        E_rad = E_rad - (E_rad - e*sin(E_rad) - M_rad) / (1 - e*cos(E_rad));
    end

    % E -> ν
    nu_rad = 2 * atan2( sqrt(1+e) * sin(E_rad/2), ...
                        sqrt(1-e) * cos(E_rad/2) );

    % 调用 orb2rv 得到 ECI 笛卡尔状态
    [r, v] = orb2rv(p, e, i_rad, Om_rad, w_rad, nu_rad, [], [], [], muEarth);
    states0(k,:) = [r', v'];
end

parameters.paths.phys_path   = '/home/user/thalassa_dir/data/physical_constants.txt';
parameters.paths.earth_path  = '/home/user/thalassa_dir/data/earth_potential/GRIM5-S1.txt';
parameters.paths.kernel_path = '/home/user/thalassa_dir/data/kernels_to_load.furnsh';
%% 构造 THALASSA parameters（注意路径需改成你实际文件位置）
parameters.model.insgrav = 1;   % 非球形重力
parameters.model.isun    = 3;   % 太阳摄动
parameters.model.imoon   = 3;   % 月球摄动
parameters.model.idrag   = 4;   % 大气阻力
parameters.model.iF107   = 1;   % 恒定 F10.7
parameters.model.iSRP    = 2;   % 太阳光压
parameters.model.iephem  = 1;   % DE431 星历
parameters.model.gdeg    = 7;
parameters.model.gord    = 7;



parameters.settings.tol    = 1e-11;
parameters.settings.imcoll = 0;
parameters.settings.eqs    = 2;

parameters.spacecraft.mass      = 150;
parameters.spacecraft.area_drag = 1.84;
parameters.spacecraft.area_srp  = 1.84;
parameters.spacecraft.cd        = 1.28;
parameters.spacecraft.cr        = 1.0;

%% 时间设置
tStart_MJD = 58834.0;
tEnd_MJD   = 58850.1;
dt_days    = 0.00001;  % 步长(天) 可调

times = tStart_MJD:dt_days:tEnd_MJD;
Nt    = length(times);

%% 传播 Target
targetOut = mthalassa(times, states0(1,:), parameters);
% 自动适配输出格式
if size(targetOut,1) ~= Nt
    targetOut = targetOut.'; % 转置
end

posTAll = targetOut(:,1:3);

%% 初始化结果
minDist = inf(1,5);
minTime = nan(1,5);
names   = {'Field1','Field2','Field3','Field4','Field5'};

%% 传播并比较其它卫星
for k = 2:6
    fieldOut = mthalassa(times, states0(k,:), parameters);
    if size(fieldOut,1) ~= Nt
        fieldOut = fieldOut.'; % 转置
    end
    posFAll = fieldOut(:,1:3);

    % 计算每步距离
    for idx = 1:Nt
        dist_km = norm(posFAll(idx,:) - posTAll(idx,:));
        if dist_km < minDist(k-1)
            minDist(k-1) = dist_km;
            minTime(k-1) = times(idx);
        end
    end
end

%% 输出结果
fprintf('\n===== 最近距离结果 =====\n');
for k = 1:5
    fprintf('%s 最近距离: %.3f km, 发生在 MJD %.5f \n', ...
            names{k}, minDist(k), minTime(k));
end
