% sweep_clutter.m — 蒙特卡洛: 不同杂波密度下的跟踪性能对比
clear; close all;

%% 参数设置
USE_FEATURES = true;
N_MC = 50;                               % 蒙特卡洛次数
lambda_c_values = [0, 1, 2, 3, 5];       % 杂波密度扫描点
colors = lines(2);                        % 目标1蓝, 目标2红

%% 蒙特卡洛仿真
fprintf('蒙特卡洛杂波密度扫描 (N=%d)...\n', N_MC);
n_lc = length(lambda_c_values);

pos_mean = NaN(2, n_lc);  pos_std = NaN(2, n_lc);
vel_mean = NaN(2, n_lc);  vel_std = NaN(2, n_lc);
gos_mean = NaN(1, n_lc);  gos_std = NaN(1, n_lc);   % GOSPA (收敛段均值)

for idx = 1:n_lc
    lc = lambda_c_values(idx);
    fprintf('  lambda_c = %d [', lc);

    pos_all = NaN(2, N_MC);
    vel_all = NaN(2, N_MC);
    gos_all = NaN(1, N_MC);

    for mc = 1:N_MC
        rng(mc);

        overrides = struct('lambda_c', lc, 'use_features', USE_FEATURES);
        params = init_parameters(overrides);
        try
            [pos_all(:, mc), vel_all(:, mc), g] = run_simulation(params);
            gos_all(mc) = g.mean_total_conv;
        catch
            pos_all(:, mc) = NaN(2, 1);
            vel_all(:, mc) = NaN(2, 1);
            gos_all(mc) = NaN;
        end

        if mod(mc, 10) == 0, fprintf('.'); end
    end

    % 统计 (排除 NaN = 未跟踪到/崩溃)
    for tgt = 1:2
        valid_pos = pos_all(tgt, ~isnan(pos_all(tgt, :)));
        valid_vel = vel_all(tgt, ~isnan(vel_all(tgt, :)));
        pos_mean(tgt, idx) = mean(valid_pos);
        pos_std(tgt, idx)  = std(valid_pos);
        vel_mean(tgt, idx) = mean(valid_vel);
        vel_std(tgt, idx)  = std(valid_vel);
    end
    n_ok_t1 = sum(~isnan(pos_all(1, :)));
    n_ok_t2 = sum(~isnan(pos_all(2, :)));

    valid_gos = gos_all(~isnan(gos_all));
    if ~isempty(valid_gos)
        gos_mean(idx) = mean(valid_gos);
        gos_std(idx)  = std(valid_gos);
    end

    fprintf('] OK=%d/%d 目标1 %.1f±%.1fm | 目标2 %.1f±%.1fm | GOSPA %.1f±%.1fm\n', ...
        n_ok_t1, N_MC, pos_mean(1, idx), pos_std(1, idx), ...
        pos_mean(2, idx), pos_std(2, idx), gos_mean(idx), gos_std(idx));
end

%% 图1: 位置RMSE vs 杂波密度 (带误差棒)
figure('Name', '杂波密度-位置RMSE', 'Position', [100, 100, 750, 450]);
for tgt = 1:2
    errorbar(lambda_c_values, pos_mean(tgt, :), pos_std(tgt, :), 'o-', ...
        'Color', colors(tgt, :), 'LineWidth', 1.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors(tgt, :), ...
        'DisplayName', sprintf('目标%d (均值±1σ)', tgt));
    hold on;
end
xlabel('杂波密度 \lambda_c (每帧每平台)');
ylabel('位置 RMSE (m)');
title(sprintf('不同杂波密度下的位置估计 RMSE (N=%d Monte Carlo)', N_MC));
legend('Location', 'best'); grid on;

%% 图2: 速度RMSE vs 杂波密度 (带误差棒)
figure('Name', '杂波密度-速度RMSE', 'Position', [100, 100, 750, 450]);
for tgt = 1:2
    errorbar(lambda_c_values, vel_mean(tgt, :), vel_std(tgt, :), 'o-', ...
        'Color', colors(tgt, :), 'LineWidth', 1.5, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors(tgt, :), ...
        'DisplayName', sprintf('目标%d (均值±1σ)', tgt));
    hold on;
end
xlabel('杂波密度 \lambda_c (每帧每平台)');
ylabel('速度 RMSE (m/s)');
title(sprintf('不同杂波密度下的速度估计 RMSE (N=%d Monte Carlo)', N_MC));
legend('Location', 'best'); grid on;

%% 图3: GOSPA vs 杂波密度 (带误差棒)
% GOSPA 同时计入漏检与虚警, 杂波增大时通常比 RMSE 上升更明显,
% 因此更能反映关联算法在密集杂波下的真实退化
figure('Name', '杂波密度-GOSPA', 'Position', [100, 100, 750, 450]);
errorbar(lambda_c_values, gos_mean, gos_std, 'o-', ...
    'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.85, 0.33, 0.10], ...
    'DisplayName', 'GOSPA (收敛段均值±1σ)');
xlabel('杂波密度 \lambda_c (每帧每平台)');
ylabel('GOSPA (m)');
title(sprintf('不同杂波密度下的 GOSPA 指标 (N=%d Monte Carlo)', N_MC));
legend('Location', 'best'); grid on;

fprintf('完成。\n');
