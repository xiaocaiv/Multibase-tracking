function [z_pred, S, P_xz] = ukf_measurement_predict(x, P, R, plat_pos, alpha, beta, kappa)
%UKF_MEASUREMENT_PREDICT UKF 量测预测步
%   [z_pred, S, P_xz] = UKF_MEASUREMENT_PREDICT(x, P, R, plat_pos, alpha, beta, kappa)
%   将 sigma 点通过非线性方位角量测函数传播, 计算:
%     z_pred: 预测量测 (标量)
%     S:      新息协方差 (标量)
%     P_xz:   状态-量测互协方差 [n x 1]

n = length(x);
[X, Wm, Wc] = generate_sigma_points(x, P, alpha, beta, kappa);
n_sigma = 2 * n + 1;

% sigma 点通过量测函数
Z = measurement_bearing(X, plat_pos);  % [1 x n_sigma]

% 预测量测 (加权均值)
z_pred = 0;
for i = 1:n_sigma
    z_pred = z_pred + Wm(i) * Z(i);
end

% 角度修正: 确保 z_pred 在 [-pi, pi]
z_pred = mod(z_pred + pi, 2 * pi) - pi;

% 新息协方差 S (标量)
S = 0;
for i = 1:n_sigma
    dz = Z(i) - z_pred;
    dz = mod(dz + pi, 2 * pi) - pi;  % 角度差归一化
    S = S + Wc(i) * dz^2;
end
S = S + R;

% 互协方差 P_xz [n x 1]
P_xz = zeros(n, 1);
for i = 1:n_sigma
    dx = X(:, i) - x;
    dz = Z(i) - z_pred;
    dz = mod(dz + pi, 2 * pi) - pi;
    P_xz = P_xz + Wc(i) * dx * dz;
end
end
