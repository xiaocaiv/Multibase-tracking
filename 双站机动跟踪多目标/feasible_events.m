function events = feasible_events(validation_matrix)
%FEASIBLE_EVENTS 枚举所有可行的联合关联事件
%   events = FEASIBLE_EVENTS(validation_matrix)
%   validation_matrix: 确认矩阵 [num_tracks x num_meas], 1 表示落入波门
%
%   events: 元胞数组, 每个事件是长度为 num_tracks 的向量
%           events{i}(t) = j 表示第 t 个目标关联第 j 个量测 (j=0 表示未检测)
%
%   约束:
%     - 每个目标至多关联 1 个量测
%     - 每个量测至多关联 1 个目标

[num_tracks, num_meas] = size(validation_matrix);
events = {};
enumerate([], 1, false(1, num_meas));

    function enumerate(assignment, track_idx, used_meas)
        if track_idx > num_tracks
            events{end + 1} = assignment;  %#ok<AGROW>
            return;
        end

        % 选项1: 目标未被检测 (j = 0)
        enumerate([assignment, 0], track_idx + 1, used_meas);

        % 选项2: 目标关联到一个落入波门的量测
        for j = 1:num_meas
            if validation_matrix(track_idx, j) && ~used_meas(j)
                new_used = used_meas;
                new_used(j) = true;
                enumerate([assignment, j], track_idx + 1, new_used);
            end
        end
    end
end
