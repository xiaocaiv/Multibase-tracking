function [beta, v_combined, spread, beta0, unassociated] = jpda_associate(...
    z_preds, Ss, P_xzs, measurements, tracks, plat_pos, params)
%JPDA_ASSOCIATE 联合概率数据关联 (JPDA), 支持频谱特征辅助
%   [beta, v_combined, spread, beta0, unassociated] = ...
%       JPDA_ASSOCIATE(z_preds, Ss, P_xzs, measurements, tracks, plat_pos, params)
%
%   输入:
%     z_preds:      各航迹预测量测 (元胞数组, 每元素为标量)
%     Ss:           各航迹新息协方差 (元胞数组, 每元素为标量)
%     P_xzs:        各航迹互协方差 (元胞数组, 每元素 [4x1])
%     measurements: 当前平台量测 (结构体数组, 字段 .bearing, .spectrum)
%     tracks:       航迹元胞数组
%     plat_pos:     当前平台位置
%     params:       系统参数结构体

num_tracks = length(z_preds);
num_meas = length(measurements);

% 预分配输出
beta = cell(num_tracks, 1);
v_combined = cell(num_tracks, 1);
spread = cell(num_tracks, 1);
beta0 = ones(num_tracks, 1);
unassociated = true(1, max(num_meas, 0));

% 无航迹或无量测: 所有量测视为杂波
if num_tracks == 0 || num_meas == 0
    for t = 1:num_tracks
        beta{t} = [1; zeros(num_meas, 1)];
        v_combined{t} = 0;
        spread{t} = 0;
    end
    return;
end

% ---------- 第1步: 计算所有航迹-量测对的新息 ----------
v = zeros(num_tracks, num_meas);
for t = 1:num_tracks
    for j = 1:num_meas
        v(t, j) = mod(measurements(j).bearing - z_preds{t} + pi, 2 * pi) - pi;
    end
end

% ---------- 第2步: 构建确认矩阵 (波门筛选) ----------
validation_matrix = false(num_tracks, num_meas);
for t = 1:num_tracks
    for j = 1:num_meas
        if v(t, j)^2 / Ss{t} < params.gate_gamma
            validation_matrix(t, j) = true;
        end
    end
end

% 全部未通过波门
if ~any(validation_matrix(:))
    for t = 1:num_tracks
        beta{t} = [1; zeros(num_meas, 1)];
        v_combined{t} = 0;
        spread{t} = 0;
    end
    return;
end

% ---------- 第3步: 计算关联似然 ----------
likelihood = zeros(num_tracks, num_meas);
for t = 1:num_tracks
    for j = 1:num_meas
        if validation_matrix(t, j)
            likelihood(t, j) = exp(-0.5 * v(t, j)^2 / Ss{t}) / sqrt(2 * pi * Ss{t});

            if params.use_features
                likelihood(t, j) = likelihood(t, j) * ...
                    compute_spectrum_likelihood(tracks{t}, measurements(j), params);
            end
        end
    end
end

% ---------- 第4步: 枚举联合事件并计算概率 ----------
events = feasible_events(validation_matrix);
num_events = length(events);
event_weight = zeros(num_events, 1);

% 杂波空间密度 (方位角空间, 体积 = 2*pi)
lambda_spatial = params.lambda_c / (2 * pi);

for e = 1:num_events
    ev = events{e};
    w = 1.0;
    num_clutter = 0;

    for t = 1:num_tracks
        j = ev(t);
        if j > 0
            w = w * params.Pd * likelihood(t, j);
        else
            w = w * (1 - params.Pd);
        end
    end

    % 统计杂波量测数
    used = false(1, num_meas);
    for t = 1:num_tracks
        if ev(t) > 0
            used(ev(t)) = true;
        end
    end
    num_clutter = sum(~used);

    % 杂波因子
    w = w * lambda_spatial^num_clutter;

    event_weight(e) = w;
end

% 归一化
total_weight = sum(event_weight);
if total_weight > 0
    event_prob = event_weight / total_weight;
else
    event_prob = ones(num_events, 1) / num_events;
end

% ---------- 第5步: 计算边缘关联概率 ----------
for t = 1:num_tracks
    beta{t} = zeros(num_meas + 1, 1);
    for e = 1:num_events
        j = events{e}(t);
        beta{t}(j + 1) = beta{t}(j + 1) + event_prob(e);
    end
    beta0(t) = beta{t}(1);

    % 组合新息
    v_combined{t} = 0;
    for j = 1:num_meas
        v_combined{t} = v_combined{t} + beta{t}(j + 1) * v(t, j);
    end

    % 散布项 (用于协方差更新)
    spread{t} = 0;
    for j = 1:num_meas
        spread{t} = spread{t} + beta{t}(j + 1) * v(t, j)^2;
    end
    spread{t} = spread{t} - v_combined{t}^2;
end

% ---------- 第6步: 识别未关联量测 ----------
for j = 1:num_meas
    total_beta = 0;
    for t = 1:num_tracks
        total_beta = total_beta + beta{t}(j + 1);
    end
    unassociated(j) = (total_beta < 0.5);
end
end

% ========================================================================
function spec_lik = compute_spectrum_likelihood(track, meas, params)
% 计算频谱特征似然: LOFAR 与 DEMON 分别计算 ρ, 相乘得到联合似然
%   避免 LOFAR 的高相关性淹没 DEMON 的区分能力

if ~track.feature_init
    spec_lik = 1.0;
    return;
end

n_lofar = params.N_lofar;
lofar_idx = 1:n_lofar;
demon_idx = n_lofar + 1 : params.N_features;

% LOFAR 段 Pearson
xl = track.spectrum(lofar_idx); xl = xl - mean(xl);
yl = meas.spectrum(lofar_idx);  yl = yl - mean(yl);
rho_lofar = (xl' * yl) / (norm(xl) * norm(yl) + 1e-10);

% DEMON 段 Pearson
xd = track.spectrum(demon_idx); xd = xd - mean(xd);
yd = meas.spectrum(demon_idx);  yd = yd - mean(yd);
rho_demon = (xd' * yd) / (norm(xd) * norm(yd) + 1e-10);

spec_lik = exp(params.feature_gamma * (max(0, rho_lofar)^2 + max(0, rho_demon)^2));
end
