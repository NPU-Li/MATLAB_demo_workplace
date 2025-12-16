function mjd = utc2mjd(y, m, d, h, mi, s)
% UTC2MJD 将UTC时间转换为简化儒略日(MJD)
% 输入：年, 月, 日, 时, 分, 秒（可以是标量或等长向量）
% 输出：简化儒略日 Modified Julian Date

% 对于1月和2月，当作前一年的13月和14月
mask = (m <= 2);
y(mask) = y(mask) - 1;
m(mask) = m(mask) + 12;

% 计算儒略日
A = floor(y / 100);
B = 2 - A + floor(A / 4);

jd = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + B - 1524.5;

% 加上日的小数
day_frac = (h + mi/60 + s/3600) / 24;
jd = jd + day_frac;

% 转换为简化儒略日
mjd = jd - 2400000.5;

end