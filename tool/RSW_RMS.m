% 理论/参考状态 (作为“真实值”)
r_true = [-3096.9430505926; 10180.0305303511; 6019.4889249234]; % km, J2000
v_true = [2.90200089988995; -1.84357274438334; 4.57259024549662]; % km/s, J2000

% 计算/预报状态 (多组数据，这里是5组)
% 格式： [MJD, X, Y, Z, Vx, Vy, Vz]
computed_data = [
    5.893000000000000E+04, -3.060120669094199E+03,  1.013567673668771E+04,  6.111831984807768E+03,  2.921017775085653E+00, -1.881392830956273E+00,  4.545146217218048E+00;
    % 在此添加更多行数据...
];

% 提取位置和速度
r_computed = computed_data(:, 2:4)';  % 每列是一个位置矢量
v_computed = computed_data(:, 5:7)';  % 每列是一个速度矢量
num_points = size(r_computed, 2);     % 数据点数量

% 初始化误差数组
errors_r = zeros(1, num_points);
errors_s = zeros(1, num_points);
errors_w = zeros(1, num_points);

% 1. 使用理论状态构建RSW坐标系的单位基向量
R_hat = r_true / norm(r_true);                     % 径向
W_hat = cross(r_true, v_true); W_hat = W_hat / norm(W_hat); % 法向
S_hat = cross(W_hat, R_hat);                      % 迹向

% 2. 循环计算每个数据点在RSW坐标系下的误差
for i = 1:num_points
    % 计算当前位置偏差
    delta_r = r_computed(:, i) - r_true;
    
    % 将偏差投影到RSW各轴上，得到误差分量
    errors_r(i) = dot(delta_r, R_hat); % 径向误差
    errors_s(i) = dot(delta_r, S_hat); % 迹向误差
    errors_w(i) = dot(delta_r, W_hat); % 法向误差
end

% 3. 计算均方根(RMS)误差
rms_r = sqrt(mean(errors_r .^ 2));
rms_s = sqrt(mean(errors_s .^ 2));
rms_w = sqrt(mean(errors_w .^ 2));

% 4. 输出结果，格式参考你提供的表格
fprintf('RMS误差分析结果：\n');
fprintf('%-30s %-12s %-12s %-12s\n', '条件', 'RMS r (km)', 'RMS s (km)', 'RMS w (km)');
fprintf('%-30s %-12.4f %-12.4f %-12.4f\n', '你的计算结果', rms_r, rms_s, rms_w);

% 5. (可选) 与文献中的参考值对比
fprintf('\n与文献参考值对比：\n');
fprintf('%-30s %-12s %-12s %-12s\n', '条件', 'RMS r (km)', 'RMS s (km)', 'RMS w (km)');
fprintf('%-30s %-12.4f %-12.4f %-12.4f\n', '未建模岁差-章动', 0.0739, 28.1057, 10.8734);
fprintf('%-30s %-12.4f %-12.4f %-12.4f\n', '已建模岁差-章动', 0.0015, 0.06931, 0.0194);
fprintf('%-30s %-12.4f %-12.4f %-12.4f\n', '你的计算结果', rms_r, rms_s, rms_w);

% 6. 计算并显示总位置误差 (3D RMS)
rms_total = sqrt(rms_r^2 + rms_s^2 + rms_w^2);
fprintf('\n总位置误差 (3D RMS): %.4f km\n', rms_total);