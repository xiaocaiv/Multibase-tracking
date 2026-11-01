function [pos_rmse, vel_rmse, gos] = run_simulation(params)
%RUN_SIMULATION 运行一次仿真, 返回两目标的位置/速度 RMSE 以及 GOSPA 统计
%   [pos_rmse, vel_rmse] = RUN_SIMULATION(params)
%   [pos_rmse, vel_rmse, gos] = RUN_SIMULATION(params)
%   pos_rmse: [2x1] 每个目标的位置 RMSE (m)
%   vel_rmse: [2x1] 每个目标的速度 RMSE (m/s)
%   gos:      EVALUATE_GOSPA 的输出结构体 (仅在请求第 3 个输出时计算)
%
%   注意: RMSE 只评价"已正确关联"的航迹, 忽略漏检目标与虚假航迹; GOSPA 把
%   两者都折算进距离代价, 更严格。蒙特卡洛统计建议同时报告两者。

gos = [];

truth = generate_truth(params);

%% 移动平台初始化
if params.mobile_platforms
    params.plat1_heading_hist = zeros(1, params.K);
    params.plat2_heading_hist = zeros(1, params.K);
    params.plat1_pos_hist = zeros(2, params.K);
    params.plat2_pos_hist = zeros(2, params.K);
    params.plat1_heading_hist(1) = params.plat1_init_heading;
    params.plat2_heading_hist(1) = params.plat2_init_heading;
    params.plat1_pos_hist(:, 1) = params.plat1_pos;
    params.plat2_pos_hist(:, 1) = params.plat2_pos;
end

tracks = {};
next_track_id = 1;
tracks_history = cell(params.K, 1);

for k = 1:params.K
    % 0. 移动平台: 覆写当前步位置
    if params.mobile_platforms
        params.plat1_pos = params.plat1_pos_hist(:, k);
        params.plat2_pos = params.plat2_pos_hist(:, k);
    end

    % 1. 生成量测
    [z1, z2] = generate_measurements(truth, params, k);

    % 2. UKF 预测 + 保存预测状态
    num_trk = length(tracks);
    x_preds = cell(num_trk, 1);
    P_preds = cell(num_trk, 1);
    for t = 1:num_trk
        [tracks{t}.x, tracks{t}.P] = ukf_predict(...
            tracks{t}.x, tracks{t}.P, params.F, params.Q, ...
            params.alpha, params.beta, params.kappa);
        x_preds{t} = tracks{t}.x;
        P_preds{t} = tracks{t}.P;
    end

    % 3. 深拷贝航迹
    tracks2 = cell(num_trk, 1);
    for t = 1:num_trk
        tracks2{t} = tracks{t};
    end

    % 4. 双平台并行处理
    [tracks,  unassoc1] = jpda_process_platform(tracks,  z1, params.plat1_pos, params);
    [tracks2, unassoc2] = jpda_process_platform(tracks2, z2, params.plat2_pos, params);

    % 5. 协方差交叉融合
    for t = 1:num_trk
        P_pred_inv = inv(P_preds{t});
        P1_inv = inv(tracks{t}.P);
        P2_inv = inv(tracks2{t}.P);

        P_new_inv = P1_inv + P2_inv - P_pred_inv;
        tracks{t}.P = inv(P_new_inv);
        tracks{t}.P = (tracks{t}.P + tracks{t}.P') / 2;
        tracks{t}.x = tracks{t}.P * (P1_inv * tracks{t}.x ...
            + P2_inv * tracks2{t}.x - P_pred_inv * x_preds{t});

        det1 = tracks{t}.beta0_last < 0.5;
        det2 = tracks2{t}.beta0_last < 0.5;
        if det1 || det2
            tracks{t}.hit_count = tracks{t}.hit_count + 1;
            tracks{t}.miss_count = 0;
        else
            tracks{t}.miss_count = tracks{t}.miss_count + 1;
        end

        if params.use_features && tracks{t}.feature_init && tracks2{t}.feature_init
            spec_avg = (tracks{t}.spectrum + tracks2{t}.spectrum) / 2;
            tracks{t}.spectrum = spec_avg / norm(spec_avg);
        end

        tracks{t} = rmfield(tracks{t}, 'beta0_last');
    end

    % 6. 航迹管理
    [tracks, next_track_id] = manage_tracks(tracks, z1, z2, ...
        unassoc1, unassoc2, params, next_track_id);

    tracks_history{k} = tracks;

    % ECM 优化: 决定下一步平台位置
    if params.mobile_platforms && k < params.K
        h1 = params.plat1_heading_hist(k);
        [h1_next, p1_next] = platform_ecm_optimize(tracks, ...
            params.plat1_pos, h1, params.plat1_speed, params);
        params.plat1_heading_hist(k+1) = h1_next;
        params.plat1_pos_hist(:, k+1) = p1_next;

        h2 = params.plat2_heading_hist(k);
        [h2_next, p2_next] = platform_ecm_optimize(tracks, ...
            params.plat2_pos, h2, params.plat2_speed, params);
        params.plat2_heading_hist(k+1) = h2_next;
        params.plat2_pos_hist(:, k+1) = p2_next;
    end
end

%% 计算 RMSE: 末帧最近邻匹配航迹到目标
pos_rmse = NaN(params.num_targets, 1);
vel_rmse = NaN(params.num_targets, 1);

final_trk = tracks_history{params.K};
confirmed_trk = {};
for t = 1:length(final_trk)
    if final_trk{t}.confirmed
        confirmed_trk{end + 1} = final_trk{t};
    end
end

assigned = false(1, length(confirmed_trk));
for tgt = 1:params.num_targets
    best_dist = inf;
    best_idx = 0;
    for t = 1:length(confirmed_trk)
        if assigned(t), continue; end
        d = norm(confirmed_trk{t}.x(1:2) - truth.pos{tgt}(:, params.K));
        if d < best_dist
            best_dist = d;
            best_idx = t;
        end
    end
    if best_idx > 0
        assigned(best_idx) = true;
        trk_id = confirmed_trk{best_idx}.id;

        est_x = NaN(4, params.K);
        for k = 1:params.K
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == trk_id
                    est_x(:, k) = trk_list{t}.x;
                    break;
                end
            end
        end

        valid = ~any(isnan(est_x), 1);
        if any(valid)
            pos_err = sqrt((est_x(1, valid) - truth.X{tgt}(1, valid)).^2 ...
                         + (est_x(2, valid) - truth.X{tgt}(2, valid)).^2);
            vel_err = sqrt((est_x(3, valid) - truth.X{tgt}(3, valid)).^2 ...
                         + (est_x(4, valid) - truth.X{tgt}(4, valid)).^2);
            pos_rmse(tgt) = mean(pos_err);
            vel_rmse(tgt) = mean(vel_err);
        end
    end
end

%% GOSPA 评价 (仅在调用方请求第 3 个输出时计算)
if nargout >= 3
    gos = evaluate_gospa(tracks_history, truth, params);
end
end
