function [tracks, unassociated] = jpda_process_platform(tracks, measurements, plat_pos, params)
%JPDA_PROCESS_PLATFORM 对单个平台执行完整 JPDA + UKF 更新流程
%   [tracks, unassociated] = JPDA_PROCESS_PLATFORM(tracks, measurements, plat_pos, params)
%
%   unassociated: 未关联量测索引 (逻辑数组)
%
%   对每个航迹:
%     1. 计算量测预测 (UKF measurement predict)
%     2. 运行 JPDA 获取边缘关联概率
%     3. 用组合新息执行 UKF 状态与协方差更新
%     4. 保存 beta0_last 供调用方做融合后的 hit/miss 合并

num_tracks = length(tracks);
num_meas = length(measurements);
unassociated = true(1, max(num_meas, 0));

if num_tracks == 0
    return;
end

% ---------- 第1步: 各航迹量测预测 ----------
z_preds = cell(num_tracks, 1);
Ss = cell(num_tracks, 1);
P_xzs = cell(num_tracks, 1);

for t = 1:num_tracks
    [z_preds{t}, Ss{t}, P_xzs{t}] = ukf_measurement_predict(...
        tracks{t}.x, tracks{t}.P, params.R, plat_pos, ...
        params.alpha, params.beta, params.kappa);
end

% ---------- 第2步: JPDA 关联 ----------
[beta, v_combined, spread, beta0, unassociated] = jpda_associate(...
    z_preds, Ss, P_xzs, measurements, tracks, plat_pos, params);

% ---------- 第3步: UKF 更新各航迹 ----------
for t = 1:num_tracks
    if num_meas == 0
        tracks{t}.beta0_last = 1.0;
        continue;
    end

    % 卡尔曼增益
    K = P_xzs{t} / Ss{t};

    % 状态更新
    tracks{t}.x = tracks{t}.x + K * v_combined{t};

    % JPDA 协方差更新
    tracks{t}.P = tracks{t}.P ...
        - (1 - beta0(t)) * (K * Ss{t} * K') ...
        + K * spread{t} * K';

    % 强制对称
    tracks{t}.P = (tracks{t}.P + tracks{t}.P') / 2;

    % 保存检测状态 (hit/miss 由调用方融合后统一处理)
    tracks{t}.beta0_last = beta0(t);

    % 更新频谱模板 (仅在特征辅助时生效)
    if params.use_features
        tracks{t} = update_track_spectrum(tracks{t}, measurements, beta{t}, params);
    end
end
end

% ========================================================================
function track = update_track_spectrum(track, measurements, beta, params)
% 利用 JPDA 权重更新航迹频谱模板 (加权平均 + 指数平滑)

num_meas = length(measurements);
if num_meas == 0, return; end

spec_weighted = zeros(params.N_features, 1);
total_w = 0;
for j = 1:num_meas
    w = beta(j + 1);  % beta(1)=beta0, beta(2:end) 对应各量测
    if w > 0.01
        spec_weighted = spec_weighted + w * measurements(j).spectrum;
        total_w = total_w + w;
    end
end
if total_w < 0.01, return; end

spec_new = spec_weighted / total_w;
spec_new = spec_new / norm(spec_new);

if ~track.feature_init
    track.spectrum = spec_new;
    track.feature_init = true;
else
    alpha = 0.3;
    track.spectrum = (1 - alpha) * track.spectrum + alpha * spec_new;
    track.spectrum = track.spectrum / norm(track.spectrum);
end
end
