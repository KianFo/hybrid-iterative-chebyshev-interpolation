% =========================================================================
% 3D PLOT: SNR vs Number of Modules vs Iterations
% =========================================================================
clc; clear; close all;

%% 1. Parameters & Optimal Bounds
% Index 1 corresponds to N=1, Index 2 to N=2, Index 3 to N=3
A_opt = [1.03,  0.98,  1.00 ]; 
B_opt = [1.05,  0.99,  1.01 ]; 

N_vals = 1:3;           % Modules to test: 1, 2, and 3
iters_vals = 1:5;       % Iterations to test: 1 to 5

sig_len = 4096;             
nyq_samples = 64;           
L = sig_len / nyq_samples;  

%% 2. Generate Symmetric Analog Signal
rng(42); 
x_rand = randn(1, sig_len);
X_f = fft(x_rand);
cutoff_idx = floor(sig_len / (2 * L)); 
X_f((cutoff_idx + 1) : (sig_len - cutoff_idx + 1)) = 0; 
x_analog = real(ifft(X_f));
x_analog = x_analog / max(abs(x_analog)); 

ignore_len = floor(0.10 * sig_len);
valid_idx = (ignore_len + 1) : (sig_len - ignore_len);
x_orig_valid = x_analog(valid_idx);
signal_power = sum(x_orig_valid.^2);

%% 3. Calculate SNR for all combinations
SNR_matrix = zeros(length(iters_vals), length(N_vals));

fprintf('Generating 3D Surface Data...\n');
for n_idx = 1:length(N_vals)
    N = N_vals(n_idx);
    A = A_opt(n_idx);
    B = B_opt(n_idx);
    
    for iter_idx = 1:length(iters_vals)
        k = iters_vals(iter_idx);
        
        x_recon = chebyshev_hybrid_interpolator(x_analog, L, N, k, A, B);
        
        err = sum((x_orig_valid - x_recon(valid_idx)).^2);
        if err > 0
            SNR_matrix(iter_idx, n_idx) = 10 * log10(signal_power / err);
        else
            SNR_matrix(iter_idx, n_idx) = 120; % Cap at 120dB if error is practically zero
        end
    end
end
fprintf('Done!\n');

%% 4. Plot the 3D Surface
figure('Name', 'System Performance 3D', 'Position', [150, 150, 800, 600]);

[N_mesh, It_mesh] = meshgrid(N_vals, iters_vals);

surf(N_mesh, It_mesh, SNR_matrix, 'FaceAlpha', 0.85, 'EdgeColor', 'k');
colormap('parula');
colorbar;

% Adjusting axes
xlabel('Number of Modules (N)', 'FontWeight', 'bold');
ylabel('Iterations', 'FontWeight', 'bold');
zlabel('SNR (dB)', 'FontWeight', 'bold');
title('Chebyshev Hybrid Interpolator Performance');

% Set axes limits and ticks so they only show integers
set(gca, 'XTick', 1:3);
set(gca, 'YTick', 1:5);
zlim([0, max(SNR_matrix(:)) + 10]);

view(-40, 30); 
grid on;

hold on;
for n_idx = 1:length(N_vals)
    for iter_idx = 1:length(iters_vals)
        snr_val = SNR_matrix(iter_idx, n_idx);
        % Only label SNR values greater than 60dB to avoid clutter
        if snr_val > 60
            text(N_vals(n_idx), iters_vals(iter_idx), snr_val + 5, ...
                sprintf('%.1f', snr_val), 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', 'r');
        end
    end
end
hold off;