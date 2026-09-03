% DEMO_SINGLE  Reconstruct one signal with the accelerated hybrid method and
% visualise the result in time and frequency.  A good first script to run.
%
% Reference: ParandehGheibi, Akhaee & Marvasti, "Design of Interpolation
% Functions Using Iterative Methods", EUSIPCO 2006.

clear; clc; close all;

%% ------------------------------------------------------------------ set-up
N   = 4096;                 % length of the high-resolution ("analog") signal
Nyq = 64;                   % number of Nyquist-rate samples
L   = N / Nyq;              % hold factor / oversampling ratio (T in samples)

p.L       = L;
p.modules = 1;              % one cosine module is enough (see the paper)
p.iters   = 2;              % two iterations already reach ~100 dB
p.method  = 'accelerated';
p.A       = 1.03;           % frame bounds (see study_frame_bounds.m)
p.B       = 1.05;

signalKind = 'tones';       % try 'bandlimited', 'tones' or 'chirp'
x = make_test_signal(signalKind, N, L, 42);

%% ------------------------------------------------------------- reconstruct
xSampled = zoh_operator(x, L);              % what a plain D/A would output
xRec     = hybrid_reconstruct(x, p);        % proposed method

fprintf('Signal            : %s\n', signalKind);
fprintf('Raw S&H SNR       : %6.2f dB\n', snr_metric(x, xSampled));
fprintf('Reconstructed SNR : %6.2f dB  (%d module, %d iterations)\n', ...
        snr_metric(x, xRec), p.modules, p.iters);

%% -------------------------------------------------------------- plot: time
tAxis = (0:N-1) / N;
figure('Name', 'Reconstruction', 'Color', 'w', 'Position', [80 80 1000 620]);

subplot(2,2,[1 2]);
plot(tAxis, x,        'LineWidth', 1.8); hold on;
stairs(tAxis, xSampled, 'LineWidth', 0.8);
plot(tAxis, xRec, '--', 'LineWidth', 1.6);
grid on; xlim([0.15 0.35]);
xlabel('time'); ylabel('amplitude');
legend('original', 'sample-and-hold', 'reconstructed', 'Location', 'best');
title(sprintf('%s  |  reconstruction SNR = %.1f dB', ...
      signalKind, snr_metric(x, xRec)));

%% ---------------------------------------------------------- plot: spectrum
f  = (0:N-1);
Xo = abs(fft(x));   Xr = abs(fft(xRec));
subplot(2,2,3);
plot(f, 20*log10(Xo + eps), 'LineWidth', 1.2); hold on;
plot(f, 20*log10(Xr + eps), '--', 'LineWidth', 1.2);
grid on; xlim([0 2*N/(2*L)]);
xlabel('FFT bin'); ylabel('magnitude (dB)');
legend('original', 'reconstructed'); title('spectrum');

%% ------------------------------------------------------------- plot: error
subplot(2,2,4);
plot(tAxis, x - xRec, 'LineWidth', 1.0); grid on;
xlabel('time'); ylabel('error'); title('reconstruction error x - xRec');
