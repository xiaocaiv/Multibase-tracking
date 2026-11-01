function [x_init, P_init] = triangulate_init(bearing1, bearing2, plat1, plat2, P0)
%TRIANGULATE_INIT 双平台三角测量法初始化目标状态
%   [x_init, P_init] = TRIANGULATE_INIT(bearing1, bearing2, plat1, plat2, P0)
%   bearing1: 平台1 方位角 (rad)
%   bearing2: 平台2 方位角 (rad)
%   plat1:    平台1 位置 [px; py]
%   plat2:    平台2 位置 [px; py]
%   P0:       初始协方差矩阵 [4x4]
%
%   求解两条方位线的交点:
%     tan(theta1) = (y - py1) / (x - px1)
%     tan(theta2) = (y - py2) / (x - px2)
%
%   返回:
%     x_init: 初始状态 [x; y; 0; 0]
%     P_init: 初始协方差 (速度分量设较大值)

% 方向向量
d1 = [cos(bearing1); sin(bearing1)];
d2 = [cos(bearing2); sin(bearing2)];

% 线性方程组: plat1 + t1*d1 = plat2 + t2*d2
% [d1, -d2] * [t1; t2] = plat2 - plat1
A = [d1, -d2];
b = plat2 - plat1;

% 检查是否平行 (行列式接近零)
if abs(det(A)) < 1e-10 || abs(bearing1 - bearing2) < 1e-6 ...
                      || abs(abs(bearing1 - bearing2) - pi) < 1e-6
    % 平行情况: 无法三角测量, 返回平台中点
    x_init = [(plat1(1) + plat2(1)) / 2;
              (plat1(2) + plat2(2)) / 2;
              0;
              0];
    P_init = P0;
    return;
end

t = A \ b;
intersection = plat1 + t(1) * d1;

x_init = [intersection(1); intersection(2); 0; 0];
P_init = P0;
end
