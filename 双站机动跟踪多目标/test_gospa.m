% test_gospa.m — GOSPA 实现正确性自检
%
% 运行方式:
%   >> test_gospa
%   或
%   matlab -batch "test_gospa"
%
% 用"独立途径"验证 GOSPA 实现, 共 6 组测试:
%   T1 已知解析值 (空集 / 完全重合 / 单点 / 全漏检 / 全虚警 / 超截断 / alpha 影响)
%   T2 可加分解一致性   total^p == loc^p + card^p
%   T3 与穷举法对比     枚举所有"部分指派", 小规模随机实例
%   T4 两种求解器一致   matchpairs vs 内置 Hungarian
%   T5 度量性质         非负 / 对称 / 三角不等式
%   T6 EVALUATE_GOSPA   漏检、虚警、确认航迹的统计是否正确
%
% 全部通过打印 ALL TESTS PASSED; 任一失败则抛错 (便于批处理判断)。

clear; clc;
fprintf('===== GOSPA 实现正确性自检 =====\n');

T = struct('n', 0, 'fail', 0);

c = 100; p = 2; a = 2;

%% ---------- T1 已知解析值 ----------
s = T;
fprintf('\n[T1] 已知解析值\n');

% 1) 两侧都空 -> 0
g = gospa_metric(zeros(2, 0), zeros(2, 0), c, p, a);
T = t_check(T, abs(g.total) < 1e-12, '空集 vs 空集 = 0', g.total, 0);

% 2) 完全重合 -> 0
X = [0 100; 0 50];
g = gospa_metric(X, X, c, p, a);
T = t_check(T, abs(g.total) < 1e-12, '完全重合 = 0', g.total, 0);

% 3) 单点距离 30 (< c) -> 30
g = gospa_metric([0; 0], [30; 0], c, p, a);
T = t_check(T, abs(g.total - 30) < 1e-9, '单点 30m = 30', g.total, 30);

% 4) 2 个目标全漏检 -> sqrt(c^2/a * 2) = c = 100
g = gospa_metric([0 100; 0 0], zeros(2, 0), c, p, a);
T = t_check(T, abs(g.total - 100) < 1e-9, '2 目标全漏检 = 100', g.total, 100);

% 5) 单点距离 150 (> c) -> 判为漏检+虚警, 退化成 c
g = gospa_metric([0; 0], [150; 0], c, p, a);
T = t_check(T, abs(g.total - 100) < 1e-9, '单点 150m 超截断 = 100', g.total, 100);
T = t_check(T, g.n_missed == 1 && g.n_false == 1, ...
    '同上应计为 1 漏检 + 1 虚警', sprintf('%d/%d', g.n_missed, g.n_false), '1/1');

% 6) 0 目标 1 虚警 -> sqrt(c^2/a) = c/sqrt(2) = 70.71
g = gospa_metric(zeros(2, 0), [10; 10], c, p, a);
T = t_check(T, abs(g.total - c / sqrt(a)) < 1e-9, ...
    '0 目标 1 虚警 = c/sqrt(alpha)', g.total, c / sqrt(a));

% 7) 1 目标 0 航迹 -> 同样 sqrt(c^2/a)
g = gospa_metric([0; 0], zeros(2, 0), c, p, a);
T = t_check(T, abs(g.total - c / sqrt(a)) < 1e-9, ...
    '1 目标 0 航迹 = c/sqrt(alpha)', g.total, c / sqrt(a));

% 8) alpha = 1 时全漏检代价翻倍: sqrt(c^2*n/1)
g = gospa_metric([0 100; 0 0], zeros(2, 0), c, p, 1);
T = t_check(T, abs(g.total - sqrt(c^2 * 2)) < 1e-9, ...
    'alpha=1 时 2 目标全漏检 = sqrt(2)*c', g.total, sqrt(2) * c);

% 9) 交叉指派: 必须选全局最优而非贪婪最近邻
Xc = [0 100; 0 0];
Yc = [95 5; 5 5];
g = gospa_metric(Xc, Yc, c, p, a);
expect = sqrt((norm(Xc(:, 1) - Yc(:, 2))^2 + norm(Xc(:, 2) - Yc(:, 1))^2));
T = t_check(T, abs(g.total - expect) < 1e-9, ...
    '交叉指派取全局最优', g.total, expect);

T = t_report(T, s, 'T1');

%% ---------- T2 可加分解一致性 ----------
s = T;
fprintf('\n[T2] 可加分解 total^p = loc^p + card^p (200 组随机)\n');
rng(7);
bad = 0;
for trial = 1:200
    n = randi(5) - 1;
    m = randi(5) - 1;
    X = rand(2, n) * 500;
    Y = rand(2, m) * 500;
    g = gospa_metric(X, Y, c, p, a);
    lhs = g.total^p;
    rhs = g.loc^p + g.card^p;
    ok = abs(lhs - rhs) < 1e-9 * max(1, lhs);
    if ~ok
        bad = bad + 1;
        fprintf('    FAIL  trial %d (n=%d,m=%d): %.6g vs %.6g\n', trial, n, m, lhs, rhs);
    end
    T = t_count(T, ok);
end
T = t_report(T, s, 'T2');

%% ---------- T3 与穷举法对比 ----------
s = T;
fprintf('\n[T3] 与穷举法对比 (300 组随机, n,m <= 5)\n');
rng(11);
max_err = 0;
for trial = 1:300
    n = randi(5);
    m = randi(5);
    X = rand(2, n) * 400;
    Y = rand(2, m) * 400;
    g = gospa_metric(X, Y, c, p, a);
    best = brute_force_gospa(X, Y, c, p, a);
    err = abs(g.total^p - best);
    max_err = max(max_err, err);
    ok = err < 1e-7 * max(1, best);
    if ~ok
        fprintf('    FAIL  trial %d (n=%d,m=%d): %.6g vs %.6g\n', trial, n, m, g.total^p, best);
    end
    T = t_count(T, ok);
end
fprintf('    最大偏差 (GOSPA^p): %.3e\n', max_err);
T = t_report(T, s, 'T3');

%% ---------- T4 两种求解器一致性 ----------
s = T;
fprintf('\n[T4] matchpairs 与内置 Hungarian 一致性 (300 组随机)\n');
if exist('matchpairs', 'file') == 0
    fprintf('    跳过: 当前 MATLAB 无 matchpairs。\n');
else
    rng(13);
    mismatch = 0;
    for trial = 1:300
        n = randi(6);
        m = randi(6);
        X = rand(2, n) * 300;
        Y = rand(2, m) * 300;
        g1 = gospa_metric(X, Y, c, p, a, 'matchpairs');
        g2 = gospa_metric(X, Y, c, p, a, 'hungarian');
        ok = abs(g1.total - g2.total) < 1e-9 * max(1, g1.total);
        if ~ok
            mismatch = mismatch + 1;
            fprintf('    FAIL  trial %d (n=%d,m=%d): %.9g vs %.9g\n', ...
                trial, n, m, g1.total, g2.total);
        end
        T = t_count(T, ok);
    end
    fprintf('    不一致实例: %d / 300\n', mismatch);
end
T = t_report(T, s, 'T4');

%% ---------- T5 度量性质 ----------
s = T;
fprintf('\n[T5] 度量性质 (200 组随机)\n');
rng(17);
sym_err = 0;
tri_viol = 0;
nonneg_viol = 0;
for trial = 1:200
    X = rand(2, randi(5)) * 300;
    Y = rand(2, randi(5)) * 300;
    Z = rand(2, randi(5)) * 300;

    gxy = gospa_metric(X, Y, c, p, a);
    gyx = gospa_metric(Y, X, c, p, a);
    gxz = gospa_metric(X, Z, c, p, a);
    gyz = gospa_metric(Y, Z, c, p, a);

    sym_err = max(sym_err, abs(gxy.total - gyx.total));
    if gxy.total < -1e-12
        nonneg_viol = nonneg_viol + 1;
    end
    if gxz.total > gxy.total + gyz.total + 1e-9
        tri_viol = tri_viol + 1;
        fprintf('    FAIL  三角不等式: G(x,z)=%.4f > G(x,y)+G(y,z)=%.4f\n', ...
            gxz.total, gxy.total + gyz.total);
    end
end
T = t_count(T, sym_err < 1e-9);
T = t_count(T, nonneg_viol == 0);
T = t_count(T, tri_viol == 0);
fprintf('    最大对称性偏差: %.3e,  三角不等式违例: %d\n', sym_err, tri_viol);
T = t_report(T, s, 'T5');

%% ---------- T6 EVALUATE_GOSPA 计数 ----------
s = T;
fprintf('\n[T6] EVALUATE_GOSPA 逐帧统计\n');
K = 5;
truth = struct();
truth.num_targets = 2;
truth.pos = {[[0; 0], [10; 0], [20; 0], [30; 0], [40; 0]], ...
             [[0; 300], [0; 300], [0; 300], [0; 300], [0; 300]]};

params = struct('gospa_c', c, 'gospa_p', p, 'gospa_alpha', a, ...
    'gospa_use_confirmed', true, 'gospa_conv_step', 1, 'dt', 1);

% 每帧: 目标1 有确认航迹 (偏差 5m), 目标2 无航迹, 另有 1 条未确认航迹
th = cell(K, 1);
for k = 1:K
    t1 = struct('x', [truth.pos{1}(:, k) + [3; 4]; 0; 0], 'confirmed', true);
    t2 = struct('x', [900; 900; 0; 0], 'confirmed', false);
    th{k} = {t1, t2};
end

gos = evaluate_gospa(th, truth, params);
T = t_check(T, all(gos.n_tracks == 1), '只统计确认航迹 (应为 1)', gos.n_tracks(1), 1);
T = t_check(T, all(gos.n_missed == 1), '目标2 全部漏检 (应为 1)', gos.n_missed(1), 1);
T = t_check(T, all(gos.n_false == 0), '无虚假航迹 (应为 0)', gos.n_false(1), 0);

expect = (5^p + c^p / a)^(1 / p);
T = t_check(T, abs(gos.mean_total - expect) < 1e-9, '总 GOSPA 解析值', ...
    gos.mean_total, expect);

params2 = params;
params2.gospa_use_confirmed = false;
gos2 = evaluate_gospa(th, truth, params2);
T = t_check(T, all(gos2.n_tracks == 2), 'use_confirmed=false 统计全部航迹 (应为 2)', ...
    gos2.n_tracks(1), 2);
T = t_check(T, all(gos2.n_false == 1), 'use_confirmed=false 时多出 1 虚警', ...
    gos2.n_false(1), 1);
T = t_report(T, s, 'T6');

%% ---------- 汇总 ----------
fprintf('\n=================================\n');
fprintf('测试项: %d,  失败: %d\n', T.n, T.fail);
if T.fail == 0
    fprintf('ALL TESTS PASSED\n');
else
    error('test_gospa: %d 项测试失败。', T.fail);
end

%% ===================== 本地函数 =====================

function T = t_count(T, ok)
T.n = T.n + 1;
if ~ok
    T.fail = T.fail + 1;
end
end

function T = t_check(T, cond, name, got, expect)
T.n = T.n + 1;
if ~cond
    T.fail = T.fail + 1;
    if nargin >= 5
        fprintf('    FAIL  %s  (得到 %.6g, 期望 %.6g)\n', name, got, expect);
    else
        fprintf('    FAIL  %s\n', name);
    end
end
end

function T = t_report(T, s, name)
n_this = T.n - s.n;
f_this = T.fail - s.fail;
if f_this == 0
    fprintf('  ==> %s: PASS (%d/%d)\n', name, n_this, n_this);
else
    fprintf('  ==> %s: FAIL (%d/%d 失败)\n', name, f_this, n_this);
end
end

function best = brute_force_gospa(X, Y, c, p, alpha)
%BRUTE_FORCE_GOSPA 独立实现: 枚举所有"部分指派"求最小
%   把部分指派编码为映射 f: {1..n} -> {0,1..m}, 其中 0 表示"该真值点不匹配",
%   并检查非零值互不相同 (单射)。共 (m+1)^n 种组合, n,m<=5 时最多 7776 种,
%   与任何指派算法都无关, 适合做交叉校验。
n = size(X, 2);
m = size(Y, 2);
w = c^p / alpha;

D = zeros(n, m);
for i = 1:n
    for j = 1:m
        D(i, j) = norm(X(:, i) - Y(:, j))^p;
    end
end

best = inf;
n_code = (m + 1)^n;
for code = 0:(n_code - 1)
    f = zeros(1, n);
    x = code;
    for i = 1:n
        f(i) = mod(x, m + 1);
        x = floor(x / (m + 1));
    end

    nz = f(f > 0);
    if numel(unique(nz)) ~= numel(nz)
        continue;                      % 非单射, 非法指派
    end

    cost = w * (n + m - 2 * numel(nz));
    for i = 1:n
        if f(i) > 0
            cost = cost + D(i, f(i));
        end
    end

    if cost < best
        best = cost;
    end
end
end
