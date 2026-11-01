function truth = generate_truth(params)
%GENERATE_TRUTH 生成目标真实轨迹 (纯匀速直线运动, 无随机扰动)
%   truth = GENERATE_TRUTH(params)
%
%   truth 结构体:
%     truth.X:      元胞数组 {num_targets}, 每个 [4 x K] 完整状态历史
%     truth.pos:     元胞数组 {num_targets}, 每个 [2 x K] 位置历史
%     truth.num_targets: 目标数

K = params.K;
num_targets = params.num_targets;

truth = struct();
truth.num_targets = num_targets;
truth.X = cell(num_targets, 1);
truth.pos = cell(num_targets, 1);

% 目标初始状态
targets_init = {params.target1_init, params.target2_init};

for tgt = 1:num_targets
    X_hist = zeros(4, K);
    X_hist(:, 1) = targets_init{tgt};

    for k = 2:K
        X_hist(:, k) = params.F * X_hist(:, k - 1);  % 纯匀速直线
    end

    truth.X{tgt} = X_hist;
    truth.pos{tgt} = X_hist(1:2, :);
end

%% 频谱特征真值 (调用 generate_ship_spectrum 生成一次, 后续帧从缓存读取)
truth.spectrum = cell(num_targets, 1);
truth.spec_info = cell(num_targets, 1);
for tgt = 1:num_targets
    [truth.spectrum{tgt}, truth.spec_info{tgt}] = ...
        generate_ship_spectrum(params, params.target_specs(tgt), ...
        params.target_specs(tgt).wav_file);
end
end
