function z = measurement_bearing(state, plat_pos)
%MEASUREMENT_BEARING 计算平台到目标的方位角 (非线性量测函数)
%   z = MEASUREMENT_BEARING(state, plat_pos)
%   state:    状态矩阵 [n x N], 每列一个状态向量, 行1/2 为位置 x/y
%   plat_pos: 平台位置 [px; py]
%
%   z: 方位角 [1 x N], 范围 [-pi, pi]

dx = state(1, :) - plat_pos(1);
dy = state(2, :) - plat_pos(2);
z = atan2(dy, dx);
end
