function [spectrum, info] = generate_ship_spectrum(params, target_specs, wav_file)
%GENERATE_SHIP_SPECTRUM 生成舰船辐射噪声的 LOFAR+DEMON 联合频谱
%   [spectrum, info] = GENERATE_SHIP_SPECTRUM(params, target_specs, wav_file)
%
%   支持两种模式:
%     - 仿真模式: 根据目标参数 (桨叶数、轴频、线谱数等) 合成频谱
%     - 文件导入: 从 .wav 音频文件提取 LOFAR + DEMON 特征
%
%   输入:
%     params:       系统参数结构体 (来自 init_parameters.m)
%                   需包含: N_lofar, N_demon, lofar_freqs, demon_freqs,
%                           sigma_spectral, clutter_spec_min/max
%     target_specs: 目标特征参数结构体
%       .blade_count    - 桨叶数 (仿真模式, 默认 4)
%       .shaft_freq     - 轴频 Hz (仿真模式, 默认 3)
%       .n_lines        - LOFAR 线谱数 (仿真模式, 默认 8)
%       .seed           - 随机种子 (保证可复现)
%     wav_file:     .wav 文件路径 (字符串); 空 '' 或不存在则使用仿真模式
%
%   输出:
%     spectrum: N_features×1 归一化联合频谱 [LOFAR(64); DEMON(32)]
%     info:     包含元数据的结构体
%       .mode          - 'simulated' 或 'from_wav'
%       .lofar         - 64×1 LOFAR 频段 (线性功率, 归一化前)
%       .demon         - 32×1 DEMON 频段 (线性功率, 归一化前)

%% 默认参数
if nargin < 3 || isempty(wav_file) || ~exist(wav_file, 'file')
    use_wav = false;
else
    use_wav = true;
end

if nargin < 2
    target_specs = struct();
end

% 填充 target_specs 默认值
if ~isfield(target_specs, 'blade_count'), target_specs.blade_count = 4; end
if ~isfield(target_specs, 'shaft_freq'),  target_specs.shaft_freq  = 3; end
if ~isfield(target_specs, 'n_lines'),     target_specs.n_lines     = 8; end
if ~isfield(target_specs, 'seed'),        target_specs.seed        = 1; end

%% 模式选择
if use_wav
    [spectrum_raw, info] = extract_from_wav(wav_file, target_specs, params);
else
    [spectrum_raw, info] = simulate_spectrum(target_specs, params);
end

%% 归一化
total_norm = norm(spectrum_raw);
if total_norm > 0
    spectrum = spectrum_raw / total_norm;
else
    spectrum = ones(params.N_features, 1) / sqrt(params.N_features);
end

info.lofar = spectrum_raw(1:params.N_lofar);
info.demon = spectrum_raw(params.N_lofar + 1:end);

end

% ========================================================================
function [spec, info] = simulate_spectrum(specs, params)
% 仿真生成 LOFAR + DEMON 频谱

% 保存并设置随机种子
rng_state = rng;
rng(specs.seed);

N_lofar = params.N_lofar;
N_demon = params.N_demon;
f_lofar = params.lofar_freqs(:);   % N_lofar × 1
f_demon = params.demon_freqs(:);   % N_demon × 1

% ---- LOFAR 仿真: 离散线谱 + 连续背景 ----
lofar_fine = linspace(f_lofar(1), f_lofar(end), 2048)';
lofar_spec = -40 * ones(size(lofar_fine));  % 背景噪声基底 (dB)

% 生成线谱: 对数均匀分布的随机频率
r = rand(specs.n_lines, 1);
line_freqs = 10.^(log10(f_lofar(1)) + r * log10(f_lofar(end) / f_lofar(1)));
line_amps  = -15 + 15 * rand(specs.n_lines, 1);  % [-15, 0] dB

% 将线谱加到频谱上 (使用高斯展宽模拟有限频率分辨率)
for i = 1:specs.n_lines
    [~, idx] = min(abs(lofar_fine - line_freqs(i)));
    sigma_f = 1.5 + 3 * rand();  % 1.5–4.5 Hz 展宽
    gaussian = exp(-(lofar_fine - lofar_fine(idx)).^2 / (2 * sigma_f^2));
    lofar_spec = lofar_spec + line_amps(i) * gaussian;
end

% 转换到线性功率，插值到标准频点
lofar_linear = 10.^(lofar_spec / 20);
lofar_binned = interp1(lofar_fine, lofar_linear, f_lofar, 'linear', 0);

% ---- DEMON 仿真: 轴频 + 叶频谐波序列 ----
f_shaft = specs.shaft_freq;
B = specs.blade_count;
f_blade = B * f_shaft; 

demon_fine = linspace(0, f_demon(end), 1024)';
demon_spec = zeros(size(demon_fine));
noise_floor = 0.02;  % 低噪声基底

% 轴频谐波 (1×, 2×, 3× shaft rate)
for n = 1:3
    fn = n * f_shaft;
    amp = 1.0 * (0.6)^(n - 1);
    demon_spec = demon_spec + amp * exp(-(demon_fine - fn).^2 / (2 * 0.3^2));
end

% 叶频及其谐波 (主要调制能量)
for n = 1:4
    fn = n * f_blade;
    amp = 1.5 * (0.5)^(n - 1);
    demon_spec = demon_spec + amp * exp(-(demon_fine - fn).^2 / (2 * 0.5^2));
end

% 添加基底噪声
demon_spec = demon_spec + noise_floor;

% 插值到标准频点
demon_binned = interp1(demon_fine, demon_spec, f_demon, 'linear', 0);

% ---- 拼接 ----
spec = [lofar_binned; demon_binned];

% 保证非负
spec = max(spec, 1e-6);

% 恢复随机状态
rng(rng_state);

% 元数据
info = struct();
info.mode = 'simulated';
info.params = specs;
info.f_lofar = f_lofar;
info.f_demon = f_demon;

end

% ========================================================================
function [spec, info] = extract_from_wav(wav_file, target_specs, params)
% 从 .wav 文件提取 LOFAR + DEMON 频谱

N_lofar = params.N_lofar;
N_demon = params.N_demon;
f_lofar_out = params.lofar_freqs(:);
f_demon_out = params.demon_freqs(:);

% ---- 读取音频 ----
try
    [y, fs] = audioread(wav_file);
catch ME
    warning('无法读取音频文件 %s: %s。回退到仿真模式。', wav_file, ME.message);
    [spec, info] = simulate_spectrum(target_specs, params);
    return;
end

% 转为单声道
if size(y, 2) > 1
    y = mean(y, 2);
end

% 限制处理时长 (前 30 秒)
max_dur = 30;
n_max = min(length(y), round(max_dur * fs));
y = y(1:n_max);

% ---- LOFAR: STFT 平均功率谱 ----
win_len = min(2048, round(fs * 0.1));
overlap = round(win_len * 0.75);
[S, f_stft] = spectrogram(y, hamming(win_len), overlap, win_len, fs);
S_mag = mean(abs(S), 2);  % 时间平均幅度谱

% 提取 [f_lofar_out(1), f_lofar_out(end)] 频段
idx_stft = f_stft >= f_lofar_out(1) & f_stft <= f_lofar_out(end);
if sum(idx_stft) > 1
    lofar_binned = interp1(f_stft(idx_stft), S_mag(idx_stft), ...
        f_lofar_out, 'linear', 0);
else
    lofar_binned = zeros(N_lofar, 1);
end

% ---- DEMON: 包络解调 ----
% 带通滤波: 提取宽带空化噪声频段
bp_low  = max(500, f_lofar_out(end) + 50);
bp_high = min(fs / 2 - 100, 5000);

if bp_high > bp_low && fs > 2 * bp_high
    % 设计带通滤波器
    [b_bp, a_bp] = butter(4, [bp_low, bp_high] / (fs / 2), 'bandpass');
    y_bp = filtfilt(b_bp, a_bp, y);

    % 包络检波 (Hilbert)
    envelope = abs(hilbert(y_bp));

    % 低通滤波 (50 Hz 截止 → 保留调制信息)
    [b_lp, a_lp] = butter(4, f_demon_out(end) / (fs / 2), 'low');
    envelope_lp = filtfilt(b_lp, a_lp, envelope);

    % 去均值 → 只保留调制分量
    envelope_ac = envelope_lp - mean(envelope_lp);

    % FFT 得到 DEMON 谱
    N_fft = max(2048, 2^nextpow2(length(envelope_ac)));
    Y_demon = abs(fft(envelope_ac, N_fft));
    f_demon = (0:N_fft - 1)' * fs / N_fft;

    % 提取 [1, f_demon_out(end)] Hz
    idx_demon = f_demon >= f_demon_out(1) & f_demon <= f_demon_out(end);
    if sum(idx_demon) > 1
        demon_binned = interp1(f_demon(idx_demon), Y_demon(idx_demon), ...
            f_demon_out, 'linear', 0);
    else
        demon_binned = zeros(N_demon, 1);
    end
else
    % 采样率不足，回退: 对全频带信号做简易 DEMON
    envelope = abs(hilbert(y));
    [b_lp, a_lp] = butter(4, f_demon_out(end) / (fs / 2), 'low');
    envelope_lp = filtfilt(b_lp, a_lp, envelope);
    envelope_ac = envelope_lp - mean(envelope_lp);

    N_fft = max(2048, 2^nextpow2(length(envelope_ac)));
    Y_demon = abs(fft(envelope_ac, N_fft));
    f_demon = (0:N_fft - 1)' * fs / N_fft;

    idx_demon = f_demon >= f_demon_out(1) & f_demon <= f_demon_out(end);
    if sum(idx_demon) > 1
        demon_binned = interp1(f_demon(idx_demon), Y_demon(idx_demon), ...
            f_demon_out, 'linear', 0);
    else
        demon_binned = zeros(N_demon, 1);
    end
end

% ---- 拼接 ----
spec = [lofar_binned(:); demon_binned(:)];
spec = max(spec, 1e-6);

% 元数据
info = struct();
info.mode = 'from_wav';
info.wav_file = wav_file;
info.fs = fs;
info.f_lofar = f_lofar_out;
info.f_demon = f_demon_out;

end
