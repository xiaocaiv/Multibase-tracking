function [X, Wm, Wc] = generate_sigma_points(x, P, alpha, beta, kappa)
%GENERATE_SIGMA_POINTS 生成无迹变换 Sigma 点 (2n+1 个)
%   [X, Wm, Wc] = GENERATE_SIGMA_POINTS(x, P, alpha, beta, kappa)
%   x:     状态向量 [n x 1]
%   P:     协方差矩阵 [n x n]
%   alpha: UT 参数, 控制 sigma 点扩散
%   beta:   UT 参数, 包含先验分布信息 (高斯分布取 2)
%   kappa:  UT 参数, 通常取 0
%
%   X:  sigma 点矩阵 [n x (2n+1)]
%   Wm: 均值权重   [1 x (2n+1)]
%   Wc: 协方差权重 [1 x (2n+1)]

n = length(x);
lambda = alpha^2 * (n + kappa) - n;

% Cholesky 分解求 sqrt(P)
L = chol(P, 'lower');

% 生成 2n+1 个 sigma 点
X = zeros(n, 2 * n + 1);
X(:, 1) = x;
sqrt_n_lambda = sqrt(n + lambda);
for i = 1:n
    X(:, i + 1)     = x + sqrt_n_lambda * L(:, i);
    X(:, i + n + 1) = x - sqrt_n_lambda * L(:, i);
end

% 权重计算
Wm = zeros(1, 2 * n + 1);
Wc = zeros(1, 2 * n + 1);
Wm(1) = lambda / (n + lambda);
Wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
Wm(2:end) = 1 / (2 * (n + lambda));
Wc(2:end) = 1 / (2 * (n + lambda));
end
