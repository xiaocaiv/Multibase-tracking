function gos = evaluate_gospa(tracks_history, truth, params)
%EVALUATE_GOSPA 逐帧计算 GOSPA 指标及其分量, 并汇总统计量
%
%   gos = EVALUATE_GOSPA(tracks_history, truth, params)
%
%   逐帧把"确认航迹的位置"作为估计点集 Y_k, 把"目标真实位置"作为真值点集
%   X_k, 调用 GOSPA_METRIC 得到该时刻的 GOSPA 及其定位/基数分量。
%
%   输入:
%     tracks_history  {K} 元胞数组, 每帧为该时刻的航迹元胞 (来自主循环)
%     truth           由 GENERATE_TRUTH 生成的真值结构体 (使用 truth.pos)
%     params          参数结构体, 需要以下字段:
%                       gospa_c             截断距离 (m)
%                       gospa_p             阶数
%                       gospa_alpha         基数分配系数
%                       gospa_use_confirmed 是否只用已确认航迹 (默认 true)
%                       gospa_conv_step     计算"收敛后均值"的起始步 (默认 20)
%
%   输出结构体 gos:
%     ---- 逐帧序列 (1×K) ----
%     gos.time        时间轴 (s)
%     gos.total       总 GOSPA (m)
%     gos.loc         定位分量 (m)
%     gos.card        基数分量 (m)
%     gos.n_missed    每帧漏检目标数
%     gos.n_false     每帧虚假航迹数
%     gos.n_tracks    每帧参与评估的航迹数
%     ---- 汇总 ----
%     gos.mean_total / gos.mean_loc / gos.mean_card   全时段平均
%     gos.mean_total_conv / gos.mean_loc_conv / gos.mean_card_conv
%                    第 conv_step 步之后 (收敛段) 的平均
%     gos.final_total 末帧总 GOSPA
%     gos.max_total   最大总 GOSPA
%     gos.mean_missed / gos.mean_false  平均漏检目标数 / 虚假航迹数
%     gos.c / gos.p / gos.alpha / gos.use_confirmed / gos.conv_step  参数副本
%
%   注意: 第 1 帧尚未建立航迹, GOSPA 必然等于 c (全部漏检), 这是正常现象。
%   因此在论文中报告指标时, 建议用 mean_total_conv (收敛段均值), 或同时给出
%   全时段曲线以体现收敛过程。

if nargin < 3
    error('evaluate_gospa: 需要 3 个输入参数。');
end

K = numel(tracks_history);

% ---- 参数读取 (带缺省值, 兼容未加 GOSPA 字段的旧 params) ----
c     = get_field(params, 'gospa_c', 100);
p     = get_field(params, 'gospa_p', 2);
alpha = get_field(params, 'gospa_alpha', 2);
use_confirmed = get_field(params, 'gospa_use_confirmed', true);
conv_step     = get_field(params, 'gospa_conv_step', 20);
dt    = get_field(params, 'dt', 1);

conv_step = min(max(round(conv_step), 1), K);

% ---- 逐帧评估 ----
total     = zeros(1, K);
loc       = zeros(1, K);
card      = zeros(1, K);
n_missed  = zeros(1, K);
n_false   = zeros(1, K);
n_tracks  = zeros(1, K);

for k = 1:K
    % 真值点集 X_k
    X = zeros(2, truth.num_targets);
    for tgt = 1:truth.num_targets
        X(:, tgt) = truth.pos{tgt}(:, k);
    end

    % 估计点集 Y_k (预分配, 避免循环内扩维)
    trk = tracks_history{k};
    Y = zeros(2, numel(trk));
    n_est = 0;
    for t = 1:numel(trk)
        if use_confirmed && ~trk{t}.confirmed
            continue;
        end
        n_est = n_est + 1;
        Y(:, n_est) = trk{t}.x(1:2);
    end
    Y = Y(:, 1:n_est);

    g = gospa_metric(X, Y, c, p, alpha);

    total(k)    = g.total;
    loc(k)      = g.loc;
    card(k)     = g.card;
    n_missed(k) = g.n_missed;
    n_false(k)  = g.n_false;
    n_tracks(k) = g.n_est;
end

% ---- 汇总 ----
idx_conv = conv_step:K;

gos = struct();
gos.time          = (1:K) * dt;
gos.total         = total;
gos.loc           = loc;
gos.card          = card;
gos.n_missed      = n_missed;
gos.n_false       = n_false;
gos.n_tracks      = n_tracks;

gos.mean_total    = mean(total);
gos.mean_loc      = mean(loc);
gos.mean_card     = mean(card);

gos.mean_total_conv = mean(total(idx_conv));
gos.mean_loc_conv   = mean(loc(idx_conv));
gos.mean_card_conv  = mean(card(idx_conv));

gos.final_total   = total(end);
gos.max_total     = max(total);
gos.mean_missed   = mean(n_missed);
gos.mean_false    = mean(n_false);

gos.c             = c;
gos.p             = p;
gos.alpha         = alpha;
gos.use_confirmed = use_confirmed;
gos.conv_step     = conv_step;
end

% ======================================================================
function val = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
