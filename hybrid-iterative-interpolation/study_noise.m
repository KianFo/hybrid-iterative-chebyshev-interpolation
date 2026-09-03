% STUDY_NOISE  Reproduce Figure 3 of the paper: behaviour under noise.
%
% A small amount of white Gaussian noise (thermal) plus uniform noise
% (quantisation) is injected at every iteration step.  The SNR first climbs as
% the interpolation converges, then rolls off as the injected noise starts to
% propagate.  The accelerated hybrid method reaches its peak sooner and stays
% more robust than the classical iterative method.
%
% The noisy iteration is written out explicitly here (rather than reusing
% hybrid_reconstruct) so the per-step noise injection is completely visible.

clear; clc; close all;

%% ---------------------------------------------------------------- settings
N   = 4096;
Nyq = 64;
L   = N / Nyq;

maxIters = 12;
nRuns    = 20;                 % average over noise realisations
sigmaW   = 2e-3;               % white Gaussian noise std per step
sigmaU   = 2e-3;               % uniform (quantisation-like) noise half-range

itAxis   = 0:maxIters;
snrHybrid    = zeros(1, maxIters+1);
snrClassical = zeros(1, maxIters+1);

%% ------------------------------------------------ operators as local handles
S  = @(v) zoh_operator(v, L);
Phyb = @(v) modular_operator(v, L, 1);     % hybrid: one module
Pcls = @(v) modular_operator(v, L, 0);     % classical: plain lowpass

% Chebyshev constants for the accelerated hybrid method
A = 1.03; B = 1.05;
rho = (B-A)/(B+A); alpha = 2/(A+B);
lambdaRelax = 0.94;                        % classical relaxation factor

%% ----------------------------------------------------- Monte-Carlo averaging
fprintf('Averaging noisy runs (%d realisations) ...\n', nRuns);
for r = 1:nRuns
    x = make_test_signal('bandlimited', N, L, r);
    noise = @() sigmaW*randn(1,N) + sigmaU*(2*rand(1,N)-1);

    % ---- accelerated hybrid, with noise added at every step --------------
    D   = Phyb(S(x));
    y2  = D;
    snrHybrid(1) = snrHybrid(1) + snr_metric(x, y2);
    y1  = D + alpha*(D - Phyb(S(D))) + noise();
    snrHybrid(2) = snrHybrid(2) + snr_metric(x, y1);
    lam = 2; y = y1;
    for n = 2:maxIters
        lam = 1/(1 - (rho^2/4)*lam);
        y = lam*(y1 + alpha*(D - Phyb(S(y1))) - y2) + y2 + noise();
        y2 = y1; y1 = y;
        snrHybrid(n+1) = snrHybrid(n+1) + snr_metric(x, y);
    end

    % ---- classical iterative, with noise added at every step -------------
    Dc = Pcls(S(x));
    xc = Dc;
    snrClassical(1) = snrClassical(1) + snr_metric(x, xc);
    for k = 1:maxIters
        xc = xc + lambdaRelax*(Dc - Pcls(S(xc))) + noise();
        snrClassical(k+1) = snrClassical(k+1) + snr_metric(x, xc);
    end
end
snrHybrid    = snrHybrid    / nRuns;
snrClassical = snrClassical / nRuns;

%% -------------------------------------------------------------------- plot
figure('Name','Noise robustness (paper Fig. 3)','Color','w','Position',[120 120 820 560]);
plot(itAxis, snrHybrid,    '-o', 'LineWidth', 2.0); hold on;
plot(itAxis, snrClassical, '-s', 'LineWidth', 2.0);
grid on;
xlabel('number of iterations'); ylabel('SNR (dB)');
title('SNR versus iterations with noise injected at each step');
legend('accelerated hybrid (1 module)', 'classical iterative', ...
       'Location', 'NorthEast');

[pk, kh] = max(snrHybrid);
fprintf('Hybrid peaks at %.1f dB after %d iterations.\n', pk, itAxis(kh));
