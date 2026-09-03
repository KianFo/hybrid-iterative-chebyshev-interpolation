% STUDY_CONVERGENCE  Reproduce Figure 2 of the paper: reconstruction SNR as a
% function of the number of iterations for four methods,
%
%     (1) classical iterative method        (plain LPF,      no acceleration)
%     (2) hybrid method                      (1 module,       no acceleration)
%     (3) accelerated hybrid method          (1 module,       Chebyshev)
%     (4) accelerated hybrid method          (2 modules,      Chebyshev)
%
% Results are averaged over several random bandlimited signals, exactly as the
% paper does (it averages over 50 signals).

clear; clc; close all;

%% ---------------------------------------------------------------- settings
N        = 4096;
Nyq      = 64;
L        = N / Nyq;
maxIters = 8;
nSignals = 20;              % averaging ensemble (paper uses 50)

itAxis = 0:maxIters;
snr1 = zeros(1, maxIters+1);   % classical iterative
snr2 = zeros(1, maxIters+1);   % hybrid, no acceleration
snr3 = zeros(1, maxIters+1);   % accelerated hybrid, 1 module
snr4 = zeros(1, maxIters+1);   % accelerated hybrid, 2 modules

%% ------------------------------------------------ average over the ensemble
fprintf('Averaging over %d signals ...\n', nSignals);
for s = 1:nSignals
    x = make_test_signal('bandlimited', N, L, s);

    base = struct('L', L, 'relax', 0.94, 'A', 1.03, 'B', 1.05);

    p1 = base; p1.method='iterative';   p1.modules=0; p1.ref=x;
    p2 = base; p2.method='iterative';   p2.modules=1; p2.ref=x;
    p3 = base; p3.method='accelerated'; p3.modules=1; p3.A=1.03; p3.B=1.05; p3.ref=x;
    p4 = base; p4.method='accelerated'; p4.modules=2; p4.A=1.00; p4.B=1.01; p4.ref=x;

    p1.iters=maxIters; [~,t1]=hybrid_reconstruct(x,p1);
    p2.iters=maxIters; [~,t2]=hybrid_reconstruct(x,p2);
    p3.iters=maxIters; [~,t3]=hybrid_reconstruct(x,p3);
    p4.iters=maxIters; [~,t4]=hybrid_reconstruct(x,p4);

    snr1 = snr1 + t1;  snr2 = snr2 + t2;
    snr3 = snr3 + t3;  snr4 = snr4 + t4;
end
snr1=snr1/nSignals; snr2=snr2/nSignals; snr3=snr3/nSignals; snr4=snr4/nSignals;

%% --------------------------------------------------------- report + figure
fprintf('\n iters |  classical |  hybrid  | accel M=1 | accel M=2\n');
fprintf('-------+------------+----------+-----------+----------\n');
for k = 1:maxIters+1
    fprintf('  %2d   |   %6.1f   |  %6.1f  |  %7.1f  |  %6.1f\n', ...
        itAxis(k), snr1(k), snr2(k), snr3(k), snr4(k));
end

figure('Name','Convergence (paper Fig. 2)','Color','w','Position',[100 100 820 560]);
plot(itAxis, snr1, '-s', 'LineWidth', 1.8); hold on;
plot(itAxis, snr2, '-^', 'LineWidth', 1.8);
plot(itAxis, snr3, '-o', 'LineWidth', 2.0);
plot(itAxis, snr4, '-d', 'LineWidth', 2.0);
grid on;
xlabel('number of iterations'); ylabel('SNR (dB)');
title('Convergence of the interpolation methods (avg. over signals)');
legend('classical iterative', 'hybrid (1 module)', ...
       'accelerated hybrid (1 module)', 'accelerated hybrid (2 modules)', ...
       'Location', 'NorthWest');
