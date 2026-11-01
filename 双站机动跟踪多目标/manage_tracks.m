function [tracks, next_id] = manage_tracks(tracks, z1, z2, unassoc1, unassoc2, ...
    params, next_id)
%MANAGE_TRACKS 航迹生命周期管理 (删除丢失航迹 + 三角测量初始化新航迹)
%   [tracks, next_id] = MANAGE_TRACKS(tracks, z1, z2, unassoc1, unassoc2, ...
%       params, next_id)
%
%   删除条件:
%     - 连续丢失次数 >= max_misses
%     - 航迹位置超出监视区域
%
%   初始化条件:
%     - 利用两平台未关联量测对进行三角测量
%     - 位置须在监视区域内
%     - 不与现有航迹重复

% ---------- 删除丢失航迹 ----------
keep_mask = true(length(tracks), 1);
for t = 1:length(tracks)
    % 连续丢失过多
    if tracks{t}.miss_count >= params.max_misses
        keep_mask(t) = false;
        continue;
    end
    % 超出监视区域
    x_region = params.surveil_region;
    if tracks{t}.x(1) < x_region(1) || tracks{t}.x(1) > x_region(2) ...
        || tracks{t}.x(2) < x_region(3) || tracks{t}.x(2) > x_region(4)
        keep_mask(t) = false;
    end
end
tracks = tracks(keep_mask);

% ---------- 确认现有航迹 ----------
for t = 1:length(tracks)
    if tracks{t}.hit_count >= params.confirm_hits
        tracks{t}.confirmed = true;
    end
end

% ---------- 三角测量初始化新航迹 ----------
unassoc_z1 = z1(unassoc1);
unassoc_z2 = z2(unassoc2);

if isempty(unassoc_z1) || isempty(unassoc_z2)
    return;
end

new_count = 0;
candidates = [];  % 暂存候选新航迹位置

for i = 1:length(unassoc_z1)
    for j = 1:length(unassoc_z2)
        if new_count >= params.max_new_tracks
            break;
        end

        [x_init, P_init] = triangulate_init(unassoc_z1(i).bearing, ...
            unassoc_z2(j).bearing, ...
            params.plat1_pos, params.plat2_pos, params.init_cov);

        % 检查位置是否在监视区域内
        xr = params.surveil_region;
        if x_init(1) < xr(1) || x_init(1) > xr(2) ...
            || x_init(2) < xr(3) || x_init(2) > xr(4)
            continue;
        end

        % 检查与现有航迹是否重复 (距离 > 200m)
        is_duplicate = false;
        for t = 1:length(tracks)
            dist = norm(x_init(1:2) - tracks{t}.x(1:2));
            if dist < 200
                is_duplicate = true;
                break;
            end
        end
        for c = 1:size(candidates, 2)
            dist = norm(x_init(1:2) - candidates(1:2, c));
            if dist < 200
                is_duplicate = true;
                break;
            end
        end
        if is_duplicate
            continue;
        end

        % 频谱一致性校验: LOFAR/DEMON 分段 (仅在特征辅助时生效)
        if params.use_features
            n_lofar = params.N_lofar;

            xl = unassoc_z1(i).spectrum(1:n_lofar); xl = xl - mean(xl);
            yl = unassoc_z2(j).spectrum(1:n_lofar); yl = yl - mean(yl);
            rho_lofar = (xl' * yl) / (norm(xl) * norm(yl) + 1e-10);

            xd = unassoc_z1(i).spectrum(n_lofar+1:end); xd = xd - mean(xd);
            yd = unassoc_z2(j).spectrum(n_lofar+1:end); yd = yd - mean(yd);
            rho_demon = (xd' * yd) / (norm(xd) * norm(yd) + 1e-10);

            % DEMON+LOFAR 均需通过 (AND): 交叉配对 DEMON 低、LOFAR 高 ← 会被滤掉
            if rho_demon < 0.3 || rho_lofar < 0.6
                continue;
            end
        end

        % 创建新航迹
        new_track.x = x_init;
        new_track.P = P_init;
        new_track.id = next_id;
        next_id = next_id + 1;
        new_track.hit_count = 0;
        new_track.miss_count = 0;
        new_track.confirmed = false;

        % 频谱初始化
        if params.use_features
            spec_avg = (unassoc_z1(i).spectrum + unassoc_z2(j).spectrum) / 2;
            new_track.spectrum = spec_avg / norm(spec_avg);
            new_track.feature_init = true;
        else
            new_track.spectrum = NaN;
            new_track.feature_init = false;
        end

        tracks{end + 1} = new_track;
        candidates(:, end + 1) = x_init(1:2);
        new_count = new_count + 1;
    end
    if new_count >= params.max_new_tracks
        break;
    end
end
end
