% =========================================================================
% FINAL VALIDATION: Alignment Check & Paper Figure 2 Reproduction
% =========================================================================
clc; clear; close all;

%% 1. System Parameters
N_modules = 1;              
iters_max = 5; % Test up to 5 iterations to plot the curve
A = 1.03;      % optimized A
B = 1.05;      % optimized B
lambda_opt = 0.94; % Optimal lambda for classical method

sig_len = 4096;             
nyq_samples = 64;           
L = sig_len / nyq_samples;  

%% 2. Generate Symmetric Analog Signal
rng(42); 
x_rand = randn(1, sig_len);
X_f = fft(x_rand);
cutoff_idx = floor(sig_len / (2 * L)); 
X_f((cutoff_idx + 1) : (sig_len - cutoff_idx + 1)) = 0; % Perfect Symmetry
x_analog = real(ifft(X_f));
x_analog = x_analog / max(abs(x_analog)); 

% Valid indices for SNR (10% margin)
ignore_len = floor(0.10 * sig_len);
valid_idx = (ignore_len + 1) : (sig_len - ignore_len);
x_orig_valid = x_analog(valid_idx);
signal_power = sum(x_orig_valid.^2);

%% 3. Run Both Methods Across Multiple Iterations
snr_chebyshev = zeros(1, iters_max + 1);
snr_classical = zeros(1, iters_max + 1);

fprintf('Testing Iterations to Replicate Figure 2...\n');

for k = 0:iters_max
    % Run Chebyshev
    x_cheb = chebyshev_hybrid_interpolator(x_analog, L, N_modules, k, A, B);
    err_cheb = sum((x_orig_valid - x_cheb(valid_idx)).^2);
    snr_chebyshev(k+1) = 10 * log10(signal_power / err_cheb);
    
    % Run Classical
    x_class = classical_hybrid_interpolator(x_analog, L, N_modules, k, lambda_opt);
    err_class = sum((x_orig_valid - x_class(valid_idx)).^2);
    snr_classical(k+1) = 10 * log10(signal_power / err_class);
end

%% 4. Generate the Final Reconstructed Signal for Error Plotting
x_final_recon = chebyshev_hybrid_interpolator(x_analog, L, N_modules, 2, A, B);
error_signal = x_analog - x_final_recon;

%% 5. Plotting the Proofs
figure('Name', 'Final Validation & Alignment Proof', 'Position', [100, 100, 1000, 500]);

% --- Plot 1: Replicating Figure 2 of the Paper ---
subplot(1, 1, 1);
iterations_vec = 0:iters_max;
plot(iterations_vec, snr_chebyshev, '-og', 'LineWidth', 2, 'MarkerFaceColor', 'g'); hold on;
plot(iterations_vec, snr_classical, '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
title('Reproduction of Paper Figure 2');
xlabel('Number of Iterations');
ylabel('SNR (dB)');
legend('Accelerated Hybrid (1 Module)', 'Classical Hybrid (1 Module)', 'Location', 'SouthEast');
ylim([10, 180]);
