function g = gospa_metric(X, Y, c, p, alpha, solver)
%GOSPA_METRIC 计算单个时刻的 GOSPA (Generalised Optimal Sub-Pattern Assignment) 距离
%
%   g = GOSPA_METRIC(X, Y, c, p, alpha)
%   g = GOSPA_METRIC(X, Y, c, p, alpha, solver)
%
%   GOSPA 是 Rahmani & Vihola (2019) 提出的多目标跟踪统一评价指标, 同时惩罚
%   三类误差: 定位误差、漏检目标 (missed target)、虚假航迹 (false track)。
%   与 RMSE 的本质区别在于: RMSE 必须"事先知道哪条航迹对应哪个目标", 只能
%   评价已正确关联的那些量; 而 GOSPA 在所有可能的指派上取最优, 因此虚假航迹
%   和漏检目标会被直接计入代价, 不会像 RMSE 那样被忽略。
%
%   定义 (阶数 p, 截断距离 c, 基数分配系数 alpha):
%
%     GOSPA = [ min_{gamma} ( sum_{(i,j) in gamma} d(x_i, y_j)^p
%                             + (c^p / alpha) * (n + m - 2*|gamma|) ) ]^(1/p)
%
%   其中 X = {x_1..x_n} 为真值点集, Y = {y_1..y_m} 为估计点集, gamma 为 X 到 Y
%   的任意"部分单射"(即每个点最多匹配一次), |gamma| 为匹配对数。没有被匹配的
%   点——漏检目标或虚假航迹——每个贡献 c^p/alpha。
%
%   输入:
%     X      d×n 真值点集, 每列一个点 (位置评价时 d = 2)
%     Y      d×m 估计点集, 每列一个点
%     c      截断距离 (>0), 单位与坐标一致。取值的物理含义: 超过 c 的误差
%            不如直接判为"漏检+虚警", 因此 c 应取典型定位误差的 3~5 倍
%     p      距离阶数, 通常取 2 (此时 GOSPA 满足度量公理)
%     alpha  基数分配系数, 通常取 2 (漏检与虚警各承担 c^p 的一半)
%     solver (可选) 指派求解器: 'auto'(默认) | 'matchpairs' | 'hungarian'
%            默认优先使用 MATLAB 内置 matchpairs, 缺失时回退到内置 Hungarian;
%            两种实现结果应完全一致, 可用 TEST_GOSPA 交叉验证
%
%   输出结构体 g:
%     g.total      总 GOSPA 距离
%     g.loc        定位分量 = (匹配对 d^p 之和)^(1/p)
%     g.card       基数分量 = ((c^p/alpha)*(n+m-2|gamma|))^(1/p)
%     g.n_matched  最优指派下的匹配对数 |gamma|
%     g.n_missed   漏检目标数 = n - |gamma|
%     g.n_false    虚假航迹数 = m - |gamma|
%     g.n_true     真值点数 n
%     g.n_est      估计点数 m
%     g.assign     n×1, 第 i 个真值点匹配到的估计点列号 (0 表示未匹配)
%     g.dist       n×1, 对应匹配距离 (未匹配处置为 c)
%     g.loc_sum    匹配距离 p 次幂之和 (内部量, 用于再组合)
%
%   注意: total^p = loc^p + card^p, 即总误差可加性分解为"定位"与"基数"两部分,
%   这正是 GOSPA 相对 RMSE 的核心优势——能把"跟踪不准"和"跟踪漏/多"区分开。
%
%   参考文献:
%     A. S. Rahmani, M. Vihola, "A Two-Level Implementation of GOSPA Metric
%     for Multi-Target Tracking Performance Evaluation", IEEE FUSION, 2019.
%
%   另见 EVALUATE_GOSPA (逐帧批量评估), TEST_GOSPA (正确性自检)。

if nargin < 3 || isempty(c),      c = 100; end
if nargin < 4 || isempty(p),      p = 2;   end
if nargin < 5 || isempty(alpha),  alpha = 2; end
if nargin < 6 || isempty(solver), solver = 'auto'; end

if c <= 0
    error('gospa_metric: 截断距离 c 必须为正数。');
end
if p <= 0
    error('gospa_metric: 阶数 p 必须为正数。');
end
if alpha <= 0
    error('gospa_metric: 基数系数 alpha 必须为正数。');
end
if size(X, 1) ~= size(Y, 1) && ~isempty(X) && ~isempty(Y)
    error('gospa_metric: 真值点集与估计点集的维数必须一致。');
end

n = size(X, 2);
m = size(Y, 2);

g = struct();
g.n_true    = n;
g.n_est     = m;
g.n_matched = 0;
g.n_missed  = n;
g.n_false   = m;
g.loc_sum   = 0;
g.loc       = 0;
g.card      = (c^p / alpha * (n + m))^(1 / p);
g.total     = g.card;
g.assign    = zeros(n, 1);
g.dist      = c * ones(n, 1);

% 两侧都为空: 无目标也无航迹, 距离为 0 (注意不能走上面那行公式)
if n == 0 && m == 0
    g.card  = 0;
    g.total = 0;
    return;
end

w = c^p / alpha;          % 单个未匹配点 (漏检或虚警) 的代价

% ---------- 1. 距离矩阵 ----------
D = pairwise_dist(X, Y);

% ---------- 2. 最优部分指派 ----------
if n > 0 && m > 0
    [row_idx, col_idx] = solve_assignment(D.^p, w, solver);
else
    row_idx = zeros(0, 1);
    col_idx = zeros(0, 1);
end

% ---------- 3. 分量与总量 ----------
n_matched = numel(row_idx);
loc_sum = 0;
for t = 1:n_matched
    i = row_idx(t);
    j = col_idx(t);
    d = D(i, j);
    loc_sum = loc_sum + d^p;
    g.assign(i) = j;
    g.dist(i)   = d;
end

n_unmatched = (n - n_matched) + (m - n_matched);

g.loc_sum   = loc_sum;
g.n_matched = n_matched;
g.n_missed  = n - n_matched;
g.n_false   = m - n_matched;

g.loc   = loc_sum^(1 / p);
g.card  = (w * n_unmatched)^(1 / p);
g.total = (loc_sum + w * n_unmatched)^(1 / p);
end

% ======================================================================
function D = pairwise_dist(X, Y)
%PAIRWISE_DIST n×m 欧氏距离矩阵 (不依赖 Statistics Toolbox 的 pdist2)
n = size(X, 2);
m = size(Y, 2);
D = zeros(n, m);
for i = 1:n
    D(i, :) = sqrt(sum((Y - X(:, i)).^2, 1));
end
end

% ======================================================================
function [rows, cols] = solve_assignment(C, w, solver)
%SOLVE_ASSIGNMENT 最小代价部分指派: 允许行/列不被匹配, 每个未匹配代价 w
%   返回匹配对的行、列下标向量 (长度 |gamma|)
%   solver: 'auto' (默认, 有 matchpairs 就用它) | 'matchpairs' | 'hungarian'
if strcmpi(solver, 'auto')
    if exist('matchpairs', 'file') ~= 0
        solver = 'matchpairs';
    else
        solver = 'hungarian';
    end
end

switch lower(solver)
    case 'matchpairs'
        if exist('matchpairs', 'file') == 0
            error('gospa_metric: 当前 MATLAB 版本没有 matchpairs, 请改用 hungarian 求解器。');
        end
        % MATLAB R2019b+ 内置, 无需工具箱
        [M, ~, ~] = matchpairs(C, w);
        rows = M(:, 1);
        cols = M(:, 2);
    case 'hungarian'
        [rows, cols] = hungarian_partial(C, w);
    otherwise
        error('gospa_metric: 未知求解器 ''%s'' (可选 auto/matchpairs/hungarian)。', solver);
end
end

% ======================================================================
function [rows, cols] = hungarian_partial(C, w)
%HUNGARIAN_PARTIAL 用"虚拟行/列"把部分指派化为方阵完美指派
%   构造 (n+m)×(n+m) 方阵:
%     [  C         w*I_n  ]   行 i 与自己的虚拟列配对 => 行 i 未匹配, 代价 w
%     [ w*I_m       0     ]   列 j 与自己的虚拟行配对 => 列 j 未匹配, 代价 w
%   非对角块必须置为大数 M, 否则任意行都能"免费"配对到别人的虚拟列,
%   最优解会退化成"全部不匹配"(这是本函数曾经的 bug, 已在 TEST_GOSPA 中覆盖)。
n = size(C, 1);
m = size(C, 2);
N = n + m;

% 大数 M: 远大于任何合法代价, 但保持有限以避免数值问题
max_valid = max([C(:); w]);
if ~isfinite(max_valid) || max_valid <= 0
    max_valid = 1;
end
M = 1e6 * max_valid * N;

A = M * ones(N, N);

% 左上: 真实匹配代价
A(1:n, 1:m) = C;

% 右上 (n×n): 行 i -> 自己的虚拟列, 代价 w; 其余 M
for i = 1:n
    A(i, m + i) = w;
end

% 左下 (m×m): 列 j -> 自己的虚拟行, 代价 w; 其余 M
for j = 1:m
    A(n + j, j) = w;
end

% 右下 (m×n): 虚拟-虚拟配对, 代价 0
A(n + 1:N, m + 1:N) = 0;

asg = hungarian_square(A);   % asg(i) = 行 i 匹配的列号

rows = zeros(0, 1);
cols = zeros(0, 1);
for i = 1:n
    j = asg(i);
    if j >= 1 && j <= m
        rows(end + 1, 1) = i;   %#ok<AGROW>
        cols(end + 1, 1) = j;   %#ok<AGROW>
    end
end
end

% ======================================================================
function asg = hungarian_square(A)
%HUNGARIAN_SQUARE Kuhn-Munkres 指派 (O(n^3))
%   asg(i) = 第 i 行分配到的列号
n = size(A, 1);

u = zeros(n + 1, 1);    % u(i+1) = 行 i 的对偶势 (u(1) 供虚拟行使用)
v = zeros(n + 1, 1);    % v(j+1) = 列 j 的对偶势 (v(1) 供虚拟列使用)
p = ones(n + 1, 1);     % p(k) = 已分配给列 k-1 的行号+1, 值为 1 表示未分配
way = ones(n + 1, 1);

for i = 1:n
    p(1) = i + 1;
    j0 = 1;
    minv = inf(n + 1, 1);
    used = false(n + 1, 1);

    while true
        used(j0) = true;
        i0 = p(j0);
        delta = inf;
        j1 = 0;
        for j = 1:n
            k = j + 1;
            if ~used(k)
                cur = A(i0 - 1, j) - u(i0) - v(k);
                if cur < minv(k)
                    minv(k) = cur;
                    way(k) = j0;
                end
                if minv(k) < delta
                    delta = minv(k);
                    j1 = k;
                end
            end
        end
        if j1 == 0
            break;
        end
        for k = 1:(n + 1)
            if used(k)
                u(p(k)) = u(p(k)) + delta;
                v(k) = v(k) - delta;
            else
                minv(k) = minv(k) - delta;
            end
        end
        j0 = j1;
        if p(j0) == 1
            break;
        end
    end

    while j0 ~= 1
        j1 = way(j0);
        p(j0) = p(j1);
        j0 = j1;
    end
end

asg = zeros(n, 1);
for k = 2:(n + 1)
    if p(k) ~= 1
        asg(p(k) - 1) = k - 1;
    end
end
end
