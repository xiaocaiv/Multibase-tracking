function [x_pred, P_pred] = ukf_predict(x, P, F, Q, alpha, beta, kappa)
%UKF_PREDICT UKF 预测步 (时间更新)
%   [x_pred, P_pred] = UKF_PREDICT(x, P, F, Q, alpha, beta, kappa)
%   将 sigma 点通过线性状态转移矩阵 F 传播, 计算预测均值和协方差

[X, Wm, Wc] = generate_sigma_points(x, P, alpha, beta, kappa);
n = length(x);
n_sigma = 2 * n + 1;

% 通过过程模型传播 sigma 点
X_pred = F * X;

% 预测均值 (加权求和)
x_pred = zeros(n, 1);
for i = 1:n_sigma
    x_pred = x_pred + Wm(i) * X_pred(:, i);
end

% 预测协方差
P_pred = zeros(n, n);
for i = 1:n_sigma
    diff = X_pred(:, i) - x_pred;
    P_pred = P_pred + Wc(i) * (diff * diff');
end
P_pred = P_pred + Q;
end
