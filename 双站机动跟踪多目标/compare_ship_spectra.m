% compare_ship_spectra.m — 对比两艘舰船的 LOFAR+DEMON 频谱
clear; close all;

%% 文件路径
wav1 = 'Cargo_40s.wav';
wav2 = 'Motor_40s.wav';

% 如果 WAV 不存在则回退到仿真模式
params = init_parameters();

%% 生成频谱
fprintf('提取频谱中...\n');
specs(1) = params.target_specs(1);  specs(1).wav_file = wav1;
specs(2) = params.target_specs(2);  specs(2).wav_file = wav2;

[F1, info1] = generate_ship_spectrum(params, specs(1), wav1);
[F2, info2] = generate_ship_spectrum(params, specs(2), wav2);

fprintf('  Cargo:  %s\n', info1.mode);
fprintf('  Motor:  %s\n', info2.mode);

%% 频谱相关性
xm1 = F1 - mean(F1);
ym2 = F2 - mean(F2);
rho_full = (xm1' * ym2) / (norm(xm1) * norm(ym2));

xm1_lofar = F1(1:params.N_lofar) - mean(F1(1:params.N_lofar));
ym2_lofar = F2(1:params.N_lofar) - mean(F2(1:params.N_lofar));
rho_lofar = (xm1_lofar' * ym2_lofar) / (norm(xm1_lofar) * norm(ym2_lofar));

xm1_demon = F1(params.N_lofar+1:end) - mean(F1(params.N_lofar+1:end));
ym2_demon = F2(params.N_lofar+1:end) - mean(F2(params.N_lofar+1:end));
rho_demon = (xm1_demon' * ym2_demon) / (norm(xm1_demon) * norm(ym2_demon));

fprintf('\nPearson 相关系数:\n');
fprintf('  LOFAR+DEMON 联合:  rho = %.4f\n', rho_full);
fprintf('  LOFAR 单独:        rho = %.4f\n', rho_lofar);
fprintf('  DEMON 单独:        rho = %.4f\n', rho_demon);

%% 绘图
lofar_idx = 1:params.N_lofar;
demon_idx = params.N_lofar + (1:params.N_demon);

figure('Name', '舰船频谱对比', 'Position', [100, 100, 1000, 500]);

% LOFAR
subplot(2, 1, 1);
semilogx(params.lofar_freqs, F1(lofar_idx), 'b-', 'LineWidth', 1.5); hold on;
semilogx(params.lofar_freqs, F2(lofar_idx), 'r--', 'LineWidth', 1.5);
xlabel('频率 (Hz)'); ylabel('归一化幅度');
title(sprintf('LOFAR 频谱 (10–500 Hz),  rho = %.4f', rho_lofar));
legend('Cargo', 'Motor', 'Location', 'best'); grid on;

% DEMON
subplot(2, 1, 2);
plot(params.demon_freqs, F1(demon_idx), 'b-', 'LineWidth', 1.5); hold on;
plot(params.demon_freqs, F2(demon_idx), 'r--', 'LineWidth', 1.5);
xlabel('频率 (Hz)'); ylabel('归一化幅度');
title(sprintf('DEMON 频谱 (1–50 Hz),  rho = %.4f', rho_demon));
legend('Cargo', 'Motor', 'Location', 'best'); grid on;

sgtitle(sprintf('两艘舰船辐射噪声频谱对比  (联合 rho = %.4f)', rho_full));
fprintf('完成。\n');
