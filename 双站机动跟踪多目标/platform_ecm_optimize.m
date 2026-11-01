function [best_heading, best_pos_next] = platform_ecm_optimize(...
    tracks, current_pos, current_heading, plat_speed, params)
%PLATFORM_ECM_OPTIMIZE 一步贪婪 ECM 行列式最小化, 优化平台航向
%   [best_heading, best_pos_next] = PLATFORM_ECM_OPTIMIZE(...
%       tracks, current_pos, current_heading, plat_speed, params)
%
%   候选航向: current_heading + [-30, -20, -10, 0, 10, 20, 30]°
%   指标: 最小化 det(P_post) ≡ 最大化 Σ log(1 + H·P_pred·Hᵀ / R)
%   无确认航迹时保持当前航向直行

deltas = params.ecm_candidate_deltas;
n_cand = length(deltas);
candidates = current_heading + deltas;

% 筛选确认航迹
conf_idx = [];
for t = 1:length(tracks)
    if tracks{t}.confirmed
        conf_idx(end + 1) = t;
    end
end

% 无确认航迹 → 直行
if isempty(conf_idx)
    best_heading = current_heading;
    best_pos_next = current_pos + plat_speed * params.dt * [cos(best_heading); sin(best_heading)];
    return;
end

% 预测所有确认航迹
n_conf = length(conf_idx);
x_preds = cell(n_conf, 1);
P_preds = cell(n_conf, 1);
for i = 1:n_conf
    t = conf_idx(i);
    [x_preds{i}, P_preds{i}] = ukf_predict(...
        tracks{t}.x, tracks{t}.P, params.F, params.Q, ...
        params.alpha, params.beta, params.kappa);
end

% 评估每个候选航向
scores = zeros(1, n_cand);
for c = 1:n_cand
    theta_c = candidates(c);
    pos_c = current_pos + plat_speed * params.dt * [cos(theta_c); sin(theta_c)];

    info_gain = 0;
    for i = 1:n_conf
        dx = x_preds{i}(1) - pos_c(1);
        dy = x_preds{i}(2) - pos_c(2);
        r2 = max(dx^2 + dy^2, 1e-6);

        P = P_preds{i};
        % s = H * P_pred * Hᵀ = (dy²·P₁₁ − 2dx·dy·P₁₂ + dx²·P₂₂) / r⁴
        s = (dy^2 * P(1,1) - 2 * dx * dy * P(1,2) + dx^2 * P(2,2)) / (r2 * r2);
        s = max(s, 0);

        % 最小距离惩罚: 阻止平台过度逼近目标 (归一化到 ~对数项量级)
        r = sqrt(r2);
        penalty = params.ecm_lambda_penalty * max(0, (params.ecm_r_min - r) / 100)^2;

        info_gain = info_gain + log(1 + s / params.R) - penalty;
    end

    scores(c) = info_gain;
end

% 选最优
[~, best_idx] = max(scores);
best_heading = candidates(best_idx);
best_pos_next = current_pos + plat_speed * params.dt * [cos(best_heading); sin(best_heading)];
end
