function plot_results(truth, tracks_history, measurements_history, params, plat_history, gos)
%PLOT_RESULTS 生成全部可视化图表
%   PLOT_RESULTS(truth, tracks_history, measurements_history, params, plat_history)
%   PLOT_RESULTS(..., gos)  gos 为 EVALUATE_GOSPA 的输出, 省略则不画 GOSPA 图
%
%   图表:
%     Fig 1: 真实轨迹 vs 估计轨迹 (二维俯视图)
%     Fig 2: 方位角量测时间序列
%     Fig 3: 位置 RMSE 随时间变化
%     Fig 4: 速度 RMSE 随时间变化 (Vx / Vy 子图)
%     Fig 5: 频谱模板对比 (真值 vs 估计, 仅 use_features=true 时)
%     Fig 6: GOSPA 及其定位/基数分量 + 漏检与虚警数 (仅传入 gos 时)

K = params.K;

if nargin < 6
    gos = [];
end

%% 提取估计航迹 (航迹到目标的关联基于最近邻)
% 收集被确认为对应真实目标的航迹
confirmed_tracks = cell(K, 1);
for k = 1:K
    trk = tracks_history{k};
    confirmed_tracks{k} = {};
    for t = 1:length(trk)
        if trk{t}.confirmed
            confirmed_tracks{k}{end + 1} = trk{t};
        end
    end
end

% 将航迹与目标关联: 在末时刻选择距离最近的两个确认航迹
final_tracks = confirmed_tracks{K};
est_x = zeros(params.num_targets, K);
est_y = zeros(params.num_targets, K);
est_vx = zeros(params.num_targets, K);
est_vy = zeros(params.num_targets, K);
track_assigned = false(1, length(final_tracks));

for tgt = 1:params.num_targets
    best_dist = inf;
    best_trk = 0;
    for t = 1:length(final_tracks)
        if track_assigned(t)
            continue;
        end
        dist = norm(final_tracks{t}.x(1:2) - truth.pos{tgt}(:, K));
        if dist < best_dist
            best_dist = dist;
            best_trk = t;
        end
    end
    if best_trk > 0
        track_assigned(best_trk) = true;
        % 回溯该航迹在所有时刻的状态
        trk_id = final_tracks{best_trk}.id;
        for k = 1:K
            found = false;
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == trk_id
                    est_x(tgt, k)  = trk_list{t}.x(1);
                    est_y(tgt, k)  = trk_list{t}.x(2);
                    est_vx(tgt, k) = trk_list{t}.x(3);
                    est_vy(tgt, k) = trk_list{t}.x(4);
                    found = true;
                    break;
                end
            end
            if ~found
                est_x(tgt, k)  = NaN;
                est_y(tgt, k)  = NaN;
                est_vx(tgt, k) = NaN;
                est_vy(tgt, k) = NaN;
            end
        end
    end
end

%% 图1: 真实轨迹 vs 估计轨迹
figure('Name', '轨迹对比', 'Position', [100, 100, 900, 700]);
hold on;

colors = lines(params.num_targets);
h_legend = gobjects(0);  % 收集图例句柄
leg_labels = {};

% 平台 (移动轨迹 或 静态标记)
if nargin >= 5 && isfield(params, 'mobile_platforms') && params.mobile_platforms
    % 移动平台: 绘制完整轨迹
    h_legend(end+1) = plot(plat_history.pos1(1, :), plat_history.pos1(2, :), 'k-', 'LineWidth', 1.5);
    leg_labels{end+1} = '平台1';
    plot(plat_history.pos1(1, 1), plat_history.pos1(2, 1), 'ks', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(plat_history.pos1(1, end), plat_history.pos1(2, end), 'ko', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'g');

    h_legend(end+1) = plot(plat_history.pos2(1, :), plat_history.pos2(2, :), 'k--', 'LineWidth', 1.5);
    leg_labels{end+1} = '平台2';
    plot(plat_history.pos2(1, 1), plat_history.pos2(2, 1), 'kd', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(plat_history.pos2(1, end), plat_history.pos2(2, end), 'ko', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'g');
else
    % 静态平台: 固定标记
    h_legend(end+1) = plot(params.plat1_pos(1), params.plat1_pos(2), 'ks', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'k');
    leg_labels{end+1} = '平台1';
    h_legend(end+1) = plot(params.plat2_pos(1), params.plat2_pos(2), 'kd', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'k');
    leg_labels{end+1} = '平台2';
end

% 起点/终点标记 (仅第1个目标各画一次入图例, 其余目标只画不入图例)
h_start = gobjects(1, params.num_targets);
h_end   = gobjects(1, params.num_targets);

for tgt = 1:params.num_targets
    % 真实轨迹
    h_legend(end+1) = plot(truth.pos{tgt}(1, :), truth.pos{tgt}(2, :), '-', ...
        'Color', colors(tgt, :), 'LineWidth', 2);
    leg_labels{end+1} = sprintf('目标%d 真实', tgt);

    % 起点 (圆圈), 终点 (方块)
    h_start(tgt) = plot(truth.pos{tgt}(1, 1), truth.pos{tgt}(2, 1), 'o', ...
        'Color', colors(tgt, :), 'MarkerSize', 8, 'MarkerFaceColor', colors(tgt, :));
    h_end(tgt) = plot(truth.pos{tgt}(1, end), truth.pos{tgt}(2, end), 's', ...
        'Color', colors(tgt, :), 'MarkerSize', 10, 'MarkerFaceColor', colors(tgt, :));

    % 估计轨迹
    valid = ~isnan(est_x(tgt, :));
    if any(valid)
        h_legend(end+1) = plot(est_x(tgt, valid), est_y(tgt, valid), '--', ...
            'Color', colors(tgt, :), 'LineWidth', 1.5);
        leg_labels{end+1} = sprintf('目标%d 估计', tgt);
    end
end

% 图例加入起点/终点标记
h_legend(end+1) = h_start(1);
leg_labels{end+1} = '起点';
h_legend(end+1) = h_end(1);
leg_labels{end+1} = '终点';

legend(h_legend, leg_labels, 'Location', 'best'); grid on; axis equal;
xlim(params.surveil_region(1:2));
ylim(params.surveil_region(3:4));
xlabel('X轴 (m)'); ylabel('Y轴 (m)');

%% 图2: 方位角量测时间序列
figure('Name', '方位角量测', 'Position', [100, 100, 900, 400]);

time = params.dt * (1:K);

subplot(2, 1, 1);
hold on;
% 真实方位角 (无噪声)
for tgt = 1:params.num_targets
    true_bearing = zeros(1, K);
    for k = 1:K
        tx = truth.X{tgt}(1, k);
        ty = truth.X{tgt}(2, k);
        true_bearing(k) = atan2(ty - params.plat1_pos(2), tx - params.plat1_pos(1));
    end
    plot(time, true_bearing * 180 / pi, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
end
% 量测
for k = 1:K
    z1 = measurements_history{k, 1};
    if ~isempty(z1)
        plot(time(k) * ones(1, length(z1)), [z1.bearing] * 180 / pi, 'k.', 'MarkerSize', 4);
    end
end
xlabel('时间 (s)'); ylabel('方位角 (deg)');
title('平台1 方位角量测'); grid on;
legend_str1 = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str1{tgt} = sprintf('目标%d 真实', tgt);
end
legend([legend_str1, {'量测'}], 'Location', 'best');

subplot(2, 1, 2);
hold on;
for tgt = 1:params.num_targets
    true_bearing = zeros(1, K);
    for k = 1:K
        tx = truth.X{tgt}(1, k);
        ty = truth.X{tgt}(2, k);
        true_bearing(k) = atan2(ty - params.plat2_pos(2), tx - params.plat2_pos(1));
    end
    plot(time, true_bearing * 180 / pi, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
end
for k = 1:K
    z2 = measurements_history{k, 2};
    if ~isempty(z2)
        plot(time(k) * ones(1, length(z2)), [z2.bearing] * 180 / pi, 'k.', 'MarkerSize', 4);
    end
end
xlabel('时间 (s)'); ylabel('方位角 (deg)');
title('平台2 方位角量测'); grid on;

%% 图3: 位置 RMSE
figure('Name', '位置 RMSE', 'Position', [100, 100, 600, 400]);

for tgt = 1:params.num_targets
    rmse = sqrt((est_x(tgt, :) - truth.pos{tgt}(1, :)).^2 + ...
                (est_y(tgt, :) - truth.pos{tgt}(2, :)).^2);
    plot(time, rmse, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;
end
xlabel('时间 (s)'); ylabel('位置 RMSE (m)');
title('位置估计 RMSE');
legend_str = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str{tgt} = sprintf('目标%d', tgt);
end
legend(legend_str, 'Location', 'best');
grid on;

%% 图4: 速度 RMSE
figure('Name', '速度 RMSE', 'Position', [100, 100, 700, 500]);

for tgt = 1:params.num_targets
    rmse_vx = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    rmse_vy = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));

    subplot(2, 1, 1);
    plot(time, rmse_vx, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;

    subplot(2, 1, 2);
    plot(time, rmse_vy, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;
end

subplot(2, 1, 1);
xlabel('时间 (s)'); ylabel('V_x RMSE (m/s)');
title('速度估计 RMSE — X 方向');
legend_str = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str{tgt} = sprintf('目标%d', tgt);
end
legend(legend_str, 'Location', 'best'); grid on;

subplot(2, 1, 2);
xlabel('时间 (s)'); ylabel('V_y RMSE (m/s)');
title('速度估计 RMSE — Y 方向');
legend(legend_str, 'Location', 'best'); grid on;

%% 图5: 频谱模板对比 (真值 vs 估计)
if isfield(params, 'use_features') && params.use_features && isfield(truth, 'spectrum')

    figure('Name', '频谱对比', 'Position', [100, 100, 1000, 250 * params.num_targets]);

    for tgt = 1:params.num_targets
        F_true = truth.spectrum{tgt};

        % 找该目标对应的估计频谱
        F_est = NaN(params.N_features, 1);
        final_trk = tracks_history{K};
        for t = 1:length(final_trk)
            if final_trk{t}.confirmed && isfield(final_trk{t}, 'feature_init') ...
                    && final_trk{t}.feature_init
                dist = norm(final_trk{t}.x(1:2) - truth.pos{tgt}(:, K));
                if dist < 300
                    F_est = final_trk{t}.spectrum;
                    break;
                end
            end
        end

        % LOFAR 子图
        subplot(params.num_targets, 2, 2*tgt - 1);
        fl = params.lofar_freqs;
        if any(~isnan(F_est))
            plot(fl, F_est(1:params.N_lofar), 'b-', 'LineWidth', 1.5); hold on;
        end
        plot(fl, F_true(1:params.N_lofar), 'r--', 'LineWidth', 1.5);
        set(gca, 'XScale', 'log');
        xlabel('频率 (Hz)'); ylabel('归一化幅度');
        title(sprintf('目标%d LOFAR (蓝=估计, 红=真值)', tgt));
        grid on;

        % DEMON 子图
        subplot(params.num_targets, 2, 2*tgt);
        fd = params.demon_freqs;
        demon_idx = params.N_lofar + (1:params.N_demon);
        if any(~isnan(F_est))
            plot(fd, F_est(demon_idx), 'b-', 'LineWidth', 1.5); hold on;
        end
        plot(fd, F_true(demon_idx), 'r--', 'LineWidth', 1.5);
        xlabel('频率 (Hz)'); ylabel('归一化幅度');
        title(sprintf('目标%d DEMON (蓝=估计, 红=真值)', tgt));
        grid on;
    end
end

%% 图6: GOSPA 指标
if ~isempty(gos)
    figure('Name', 'GOSPA 指标', 'Position', [120, 120, 900, 620]);

    % 上: 总 GOSPA 及其可加分解
    subplot(2, 1, 1);
    plot(gos.time, gos.total, 'k-',  'LineWidth', 2.0); hold on;
    plot(gos.time, gos.loc,   'b--', 'LineWidth', 1.3);
    plot(gos.time, gos.card,  'r-.', 'LineWidth', 1.3);
    xlabel('时间 (s)'); ylabel('GOSPA (m)');
    title(sprintf('GOSPA 及其分量 (p=%d, \\alpha=%d, c=%.0f m)', ...
        gos.p, gos.alpha, gos.c));
    legend({'总 GOSPA', '定位分量', '基数分量'}, 'Location', 'best');
    grid on;

    % 下: 基数误差来源
    subplot(2, 1, 2);
    stairs(gos.time, gos.n_missed, 'r-', 'LineWidth', 1.5); hold on;
    stairs(gos.time, gos.n_false,  'm-', 'LineWidth', 1.5);
    xlabel('时间 (s)'); ylabel('个数');
    title('漏检目标数与虚假航迹数');
    legend({'漏检目标 (missed)', '虚假航迹 (false)'}, 'Location', 'best');
    ylim([0, max([1, gos.n_missed, gos.n_false]) + 0.5]);
    grid on;
end

%% 打印平均 RMSE 汇总
fprintf('\n========== 平均 RMSE 汇总 ==========\n');
fprintf('%-8s %-16s %-16s %-16s\n', '目标', '位置RMSE(m)', '速度Vx RMSE(m/s)', '速度Vy RMSE(m/s)');
fprintf('-----------------------------------------------\n');

for tgt = 1:params.num_targets
    pos_rmse = sqrt((est_x(tgt, :) - truth.X{tgt}(1, :)).^2 + ...
                    (est_y(tgt, :) - truth.X{tgt}(2, :)).^2);
    vx_rmse  = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    vy_rmse  = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));

    % 仅用有效估计值 (航迹存在的时间段)
    valid_pos = ~isnan(pos_rmse);
    valid_vx  = ~isnan(vx_rmse);
    valid_vy  = ~isnan(vy_rmse);

    avg_pos = mean(pos_rmse(valid_pos));
    avg_vx  = mean(vx_rmse(valid_vx));
    avg_vy  = mean(vy_rmse(valid_vy));

    fprintf('目标%-3d  %-16.2f %-16.2f %-16.2f\n', tgt, avg_pos, avg_vx, avg_vy);
end

% 总体平均值
all_pos = []; all_vx = []; all_vy = [];
for tgt = 1:params.num_targets
    pos_rmse = sqrt((est_x(tgt, :) - truth.X{tgt}(1, :)).^2 + ...
                    (est_y(tgt, :) - truth.X{tgt}(2, :)).^2);
    vx_rmse  = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    vy_rmse  = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));
    all_pos = [all_pos, pos_rmse(~isnan(pos_rmse))];  %#ok<AGROW>
    all_vx  = [all_vx,  vx_rmse(~isnan(vx_rmse))];    %#ok<AGROW>
    all_vy  = [all_vy,  vy_rmse(~isnan(vy_rmse))];    %#ok<AGROW>
end
fprintf('-----------------------------------------------\n');
fprintf('%-8s %-16.2f %-16.2f %-16.2f\n', '总体平均', mean(all_pos), mean(all_vx), mean(all_vy));
fprintf('===============================================\n');

%% 打印 GOSPA 汇总
if ~isempty(gos)
    if gos.use_confirmed
        scope = '仅确认航迹';
    else
        scope = '全部航迹';
    end
    fprintf('\n========== GOSPA 汇总 ==========\n');
    fprintf('参数: p=%d, alpha=%d, 截断 c=%.0f m, 评估对象=%s\n', ...
        gos.p, gos.alpha, gos.c, scope);
    fprintf('%-24s %-14s %-14s %-14s\n', '统计区间', '总GOSPA(m)', '定位分量(m)', '基数分量(m)');
    fprintf('------------------------------------------------------------\n');
    fprintf('%-24s %-14.2f %-14.2f %-14.2f\n', ...
        sprintf('全时段 (1~%d)', K), ...
        gos.mean_total, gos.mean_loc, gos.mean_card);
    fprintf('%-24s %-14.2f %-14.2f %-14.2f\n', ...
        sprintf('收敛段 (%d~%d)', gos.conv_step, K), ...
        gos.mean_total_conv, gos.mean_loc_conv, gos.mean_card_conv);
    fprintf('%-24s %-14.2f %-14.2f %-14.2f\n', ...
        '末帧', gos.final_total, gos.loc(end), gos.card(end));
    fprintf('------------------------------------------------------------\n');
    fprintf('最大总 GOSPA: %.2f m\n', gos.max_total);
    fprintf('平均漏检目标数: %.2f 个/帧   平均虚假航迹数: %.2f 条/帧\n', ...
        gos.mean_missed, gos.mean_false);
    fprintf('说明: RMSE 只评价已正确关联的航迹, 忽略漏检与虚警; GOSPA 把三者\n');
    fprintf('      统一折算为距离代价, 因此数值更大且更严格, 两者应配合使用。\n');
    fprintf('===============================================\n');
end

fprintf('可视化完成。\n');
end
