function [z1, z2] = generate_measurements(truth, params, k)
%GENERATE_MEASUREMENTS 生成当前时刻的双平台方位角量测 (含杂波 + 频谱特征)
%   [z1, z2] = GENERATE_MEASUREMENTS(truth, params, k)
%
%   输出 (结构体数组, 可能为空 0x0):
%     z1(j).bearing  — 平台1 第j个量测的方位角 (rad)
%     z1(j).spectrum — 平台1 第j个量测的 96x1 含噪频谱
%     z2: 同上, 平台2

num_targets = truth.num_targets;

% 初始化空结构体数组
z1 = struct('bearing', {}, 'spectrum', {});
z2 = struct('bearing', {}, 'spectrum', {});

%% 目标量测
for tgt = 1:num_targets
    tx = truth.X{tgt}(1, k);
    ty = truth.X{tgt}(2, k);

    % 真实方位角
    bearing1 = atan2(ty - params.plat1_pos(2), tx - params.plat1_pos(1));
    bearing2 = atan2(ty - params.plat2_pos(2), tx - params.plat2_pos(1));

    % 从缓存取频谱模板 (generate_truth 中已生成)
    F_tgt = truth.spectrum{tgt};

    % ---- 平台1 ----
    if rand() < params.Pd
        m.bearing = mod(bearing1 + params.sigma_bearing * randn() + pi, 2*pi) - pi;

        if params.use_features
            m.spectrum = add_spectrum_noise(F_tgt, params);
        else
            m.spectrum = [];
        end
        z1(end + 1) = m;  %#ok<AGROW>
    end

    % ---- 平台2 ----
    if rand() < params.Pd
        m.bearing = mod(bearing2 + params.sigma_bearing * randn() + pi, 2*pi) - pi;

        if params.use_features
            m.spectrum = add_spectrum_noise(F_tgt, params);
        else
            m.spectrum = [];
        end
        z2(end + 1) = m;  %#ok<AGROW>
    end
end

%% 杂波量测
n_clutter1 = poissrnd(params.lambda_c);
n_clutter2 = poissrnd(params.lambda_c);

for c = 1:n_clutter1
    m.bearing = params.clutter_range(1) + ...
        (params.clutter_range(2) - params.clutter_range(1)) * rand();
    if params.use_features
        m.spectrum = generate_clutter_spectrum(params);
    else
        m.spectrum = [];
    end
    z1(end + 1) = m;  %#ok<AGROW>
end

for c = 1:n_clutter2
    m.bearing = params.clutter_range(1) + ...
        (params.clutter_range(2) - params.clutter_range(1)) * rand();
    if params.use_features
        m.spectrum = generate_clutter_spectrum(params);
    else
        m.spectrum = [];
    end
    z2(end + 1) = m;  %#ok<AGROW>
end
end

% ========================================================================
function F_noisy = add_spectrum_noise(F_tgt, params)
% 对频谱模板加噪声 (dB 域), 模拟谱估计误差

% 转到 dB 域
F_dB = 20 * log10(max(F_tgt, 1e-10));

% 添加高斯噪声
F_dB_noisy = F_dB + params.sigma_spectral * randn(params.N_features, 1);

% 微小频率偏移: 对相邻 bin 做轻度平滑 (模拟 FFT 频率分辨率误差)
if rand() > 0.5
    alpha = 0.05 * rand();
    F_dB_noisy = F_dB_noisy + alpha * [0; F_dB_noisy(1:end-1) - F_dB_noisy(1:end-1)];
end

% 转回线性并归一化
F_noisy = 10.^(F_dB_noisy / 20);
F_noisy = F_noisy / (norm(F_noisy) + 1e-10);
end

% ========================================================================
function F_clutter = generate_clutter_spectrum(params)
% 生成杂波频谱: 各频点独立均匀采样, 归一化后近似平坦谱

F_dB = params.clutter_spec_min + ...
    (params.clutter_spec_max - params.clutter_spec_min) * rand(params.N_features, 1);
F_clutter = 10.^(F_dB / 20);
F_clutter = F_clutter / (norm(F_clutter) + 1e-10);
end
