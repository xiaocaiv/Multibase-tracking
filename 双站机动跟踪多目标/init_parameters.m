function params = init_parameters(overrides)
%INIT_PARAMETERS 初始化双平台协同跟踪系统全部参数
%   params = INIT_PARAMETERS() 使用默认参数
%   params = INIT_PARAMETERS(overrides) 用 overrides 结构体的字段覆盖默认值

if nargin < 1, overrides = struct(); end

%% 仿真参数
params.T_sim = 100;             % 仿真总时长 (s)
params.dt = 1;                  % 时间步长 (s)
params.K = params.T_sim / params.dt;  % 总步数

%% 观测平台参数
params.plat1_pos = [0; 0];      % 平台1 位置 (m)
params.plat2_pos = [1000; 0];   % 平台2 位置 (m)

%% 移动平台参数
params.mobile_platforms = true;          % true = 移动平台 + ECM 优化
params.plat1_speed = 6;                  % 平台1 航速 (m/s)
params.plat2_speed = 6;                  % 平台2 航速 (m/s)
params.plat1_init_heading = pi/4;        % 平台1 初始航向 (指向监视区中心)
params.plat2_init_heading = 3*pi/4;      % 平台2 初始航向
params.max_turn_rate = 30 * pi / 180;    % 最大转弯率 (rad/s)
params.ecm_candidate_deltas = (-30:10:30) * pi / 180;  % ECM 候选航向偏移
params.ecm_r_min = 300;                  % 最小逼近距离 (m), 防止绕圈
params.ecm_lambda_penalty = 0.5;         % 距离惩罚强度

%% 目标参数
params.num_targets = 2;

% 目标1: 初始状态 [x; y; vx; vy]
params.target1_init = [200; 300; 1; 5];
% 目标2: 初始状态 —— 与目标1交叉
params.target2_init = [800; 300; 0; 8];

%% 运动模型参数 (CV - 匀速运动)
% 状态转移矩阵 F
dt = params.dt;
params.F = [1 0 dt 0;
            0 1 0 dt;
            0 0 1  0;
            0 0 0  1];

% 过程噪声协方差 (离散白噪声加速度模型)
sigma_a = 0.5;  % 加速度噪声标准差 (m/s^2)
q11 = dt^3 / 3 * sigma_a^2;
q12 = dt^2 / 2 * sigma_a^2;
q22 = dt      * sigma_a^2;
params.Q = [q11, 0,   q12, 0;
            0,   q11, 0,   q12;
            q12, 0,   q22, 0;
            0,   q12, 0,   q22];

%% 量测参数
params.sigma_bearing = 1 * pi / 180;  % 方位角标准差 1°
params.R = params.sigma_bearing^2;    % 量测噪声协方差 (标量)

%% 检测与杂波
params.Pd = 0.9;                % 检测概率
params.lambda_c = 3;            % 每帧每平台杂波数期望值
params.clutter_range = [-pi, pi];  % 杂波方位角范围
params.surveil_region = [-500, 1500, -500, 1500];  % [xmin, xmax, ymin, ymax]

%% 频谱特征参数 (LOFAR + DEMON 联合)
params.use_features = true;           % 主开关: false 时完全兼容原系统

% 频谱结构定义
params.N_lofar = 64;
params.N_demon = 32;
params.N_features = params.N_lofar + params.N_demon;
params.lofar_freqs = logspace(log10(10), log10(500), params.N_lofar)';
params.demon_freqs = linspace(1, 50, params.N_demon)';

% 量测噪声
params.sigma_spectral = 3;            % 频谱高斯噪声标准差 (dB)

% 相似度参数
params.feature_gamma = 8;             % 锐化因子 (分段 LOFAR/DEMON)

% 杂波频谱参数
params.clutter_spec_min = -60;        % 杂波谱功率下界 (dB)
params.clutter_spec_max = 0;          % 杂波谱功率上界 (dB)

% 目标特征参数
% 目标特征参数
% 仿真模式 (空字符串)
% params.target_specs(1) = struct('blade_count', 4, 'shaft_freq', 3, ...
%     'n_lines', 8, 'seed', 100, 'wav_file', '');
% params.target_specs(2) = struct('blade_count', 5, 'shaft_freq', 4.5, ...
%     'n_lines', 10, 'seed', 200, 'wav_file', '');

% WAV 文件导入模式
params.target_specs(1) = struct('blade_count', 4, 'shaft_freq', 3, ...
    'n_lines', 8, 'seed', 100, 'wav_file', 'Cargo_40s.wav');
params.target_specs(2) = struct('blade_count', 5, 'shaft_freq', 4.5, ...
    'n_lines', 10, 'seed', 200, 'wav_file', 'Motor_40s.wav');

%% UKF 参数
params.alpha = 1e-3;
params.beta = 2;
params.kappa = 0;
params.state_dim = 4;

%% JPDA 参数
params.Pg = 0.99;               % 波门概率
params.gate_gamma = chi2inv(params.Pg, 1);  % 1-DOF 波门阈值

%% 航迹管理参数
params.confirm_hits = 3;        % 确认航迹所需连续关联次数
params.max_misses = 3;          % 航迹删除所需连续丢失次数
params.init_cov = diag([1000, 1000, 100, 100]);  % 新航迹初始协方差
params.max_new_tracks = 3;      % 每步最多新建航迹数

%% GOSPA 评价参数
% GOSPA 同时惩罚"定位误差""漏检目标""虚假航迹", 用来弥补 RMSE 只评价
% 已正确关联航迹的缺陷 (RMSE 会忽略漏检与虚警, 因而偏乐观)
params.gospa_c = 100;           % 截断距离 (m): 建议取典型位置误差的 3~5 倍;
                                % 超过 c 的偏差不如直接判为"漏检 + 虚警"
params.gospa_p = 2;             % 距离阶数 (p=2 时 GOSPA 满足度量公理)
params.gospa_alpha = 2;         % 基数分配系数: 漏检/虚警各承担 c^p 的一半
params.gospa_use_confirmed = true;  % true = 只用已确认航迹参与评估
params.gospa_conv_step = 20;    % "收敛后均值"的起始步 (之前航迹尚未建立)

%% 应用覆盖参数
fns = fieldnames(overrides);
for i = 1:length(fns)
    if isfield(params, fns{i})
        params.(fns{i}) = overrides.(fns{i});
    end
end
end
