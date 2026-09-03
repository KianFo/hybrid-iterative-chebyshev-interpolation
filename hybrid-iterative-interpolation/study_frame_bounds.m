% STUDY_FRAME_BOUNDS  Search for the frame bounds A and B that maximise the
% SNR of the Chebyshev-accelerated hybrid method.  The paper notes there is no
% closed-form optimum pair (A, B) -- they are chosen experimentally, once,
% before the system is run.  This script performs that grid search and draws
% the resulting SNR surface.

clear; clc; close all;

%% ---------------------------------------------------------------- settings
N   = 4096;
Nyq = 64;
L   = N / Nyq;

modules = 1;
iters   = 2;

x = make_test_signal('bandlimited', N, L, 42);

Avals = 0.90 : 0.005 : 1.10;
Bvals = 0.90 : 0.005 : 1.20;
SNR   = nan(numel(Avals), numel(Bvals));

bestSNR = -Inf; bestA = NaN; bestB = NaN;

%% -------------------------------------------------------------- grid search
fprintf('Searching %d x %d (A,B) pairs ...\n', numel(Avals), numel(Bvals));
p = struct('L',L,'modules',modules,'iters',iters,'method','accelerated');
for i = 1:numel(Avals)
    for j = 1:numel(Bvals)
        if Bvals(j) <= Avals(i), continue; end   % frame bounds need B > A
        p.A = Avals(i); p.B = Bvals(j);
        xr = hybrid_reconstruct(x, p);
        SNR(i,j) = snr_metric(x, xr);
        if SNR(i,j) > bestSNR
            bestSNR = SNR(i,j); bestA = Avals(i); bestB = Bvals(j);
        end
    end
end

fprintf('Best SNR = %.2f dB at A = %.3f, B = %.3f\n', bestSNR, bestA, bestB);

%% -------------------------------------------------------------------- plot
figure('Name','Frame-bound search','Color','w','Position',[120 120 820 600]);
[Bm, Am] = meshgrid(Bvals, Avals);
surf(Bm, Am, SNR, 'EdgeColor', 'none'); colorbar; colormap turbo;
xlabel('B'); ylabel('A'); zlabel('SNR (dB)');
title(sprintf('SNR surface  (max %.1f dB at A=%.3f, B=%.3f)', bestSNR, bestA, bestB));
view(-40, 35); hold on;
plot3(bestB, bestA, bestSNR, 'p', 'MarkerSize', 16, ...
      'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
