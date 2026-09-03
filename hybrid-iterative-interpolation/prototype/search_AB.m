% =========================================================================
% GRID SEARCH: Finding Optimum A and B for Chebyshev Hybrid Interpolator
% =========================================================================
clc; clear; close all;

%% 1. System Parameters
% N: Number of cosine modules (Set to 1 as recommended for best efficiency)
N_modules = 3;              

% iters: Number of Chebyshev iterations (Set to 2 as it yields great results)
iters = 2;                  

% Signal Generation Parameters
sig_len = 4096;             % Total length of the high-resolution analog signal
nyq_samples = 64;           % Number of discrete samples (Nyquist rate)
L = sig_len / nyq_samples;  % Hold Factor (Ratio of analog to sampled rate, T)

%% 2. Generate the "Analog" Baseband Signal
rng(42); 
x_rand = randn(1, sig_len);

X_f = fft(x_rand);
cutoff_idx = floor(sig_len / (2 * L)); 

% SYMMETRIC ZEROING FOR SIGNAL GENERATION:
X_f((cutoff_idx + 1) : (sig_len - cutoff_idx + 1)) = 0; 

x_analog = real(ifft(X_f));
x_analog = x_analog / max(abs(x_analog)); % Normalize

%% 3. Define the Search Space for A and B
% We will test A from 0.1 to 1.5, and B from 0.5 to 3.0
% In frame theory, B is typically strictly greater than A (B > A)
A_vals = 0.810 : 0.002 : 1.100;
B_vals = 0.850 : 0.002 : 1.150;

% Matrices to store the results for our 3D plot
SNR_matrix = zeros(length(A_vals), length(B_vals));

% Variables to keep track of the absolute best results
best_SNR = -inf;
best_A = 0;
best_B = 0;

%% 4. Run the Grid Search
fprintf('Starting Grid Search for Optimum A and B...\n');
fprintf('Testing %d combinations...\n', length(A_vals) * length(B_vals));

% Define the 10% margin to ignore during SNR calculation (to avoid edge transients)
ignore_len = floor(0.10 * sig_len);
valid_idx = (ignore_len + 1) : (sig_len - ignore_len);
x_orig_valid = x_analog(valid_idx);
signal_power = sum(x_orig_valid.^2);

tic; % Start a timer
for i = 1:length(A_vals)
    for j = 1:length(B_vals)
        
        A = A_vals(i);
        B = B_vals(j);
        
        % Mathematical constraint for Chebyshev: B must be greater than A
        if B <= A
            SNR_matrix(i, j) = NaN; % Ignore invalid combinations
            continue; 
        end
        
        % Run the system with current A and B
        x_recon = chebyshev_hybrid_interpolator(x_analog, L, N_modules, iters, A, B);
        
        % Calculate SNR
        x_recon_valid = x_recon(valid_idx);
        error_power = sum((x_orig_valid - x_recon_valid).^2);
        
        % Prevent log(0) if error is magically zero, otherwise calculate normally
        if error_power == 0
            current_SNR = 200; % Arbitrary high limit
        else
            current_SNR = 10 * log10(signal_power / error_power);
        end
        
        % Store in matrix for plotting
        SNR_matrix(i, j) = current_SNR;
        
        % Update best parameters if this SNR is higher than the previous best
        if current_SNR > best_SNR
            best_SNR = current_SNR;
            best_A = A;
            best_B = B;
        end
    end
end
search_time = toc; % Stop timer

%% 5. Display Results
fprintf('\nSearch Completed in %.2f seconds!\n', search_time);
fprintf('----------------------------------------\n');
fprintf('OPTIMUM PARAMETERS FOUND:\n');
fprintf('Best A : %.2f\n', best_A);
fprintf('Best B : %.2f\n', best_B);
fprintf('Max SNR: %.2f dB\n', best_SNR);
fprintf('----------------------------------------\n');

%% 6. Plot the Performance Surface (3D Graph)
figure('Name', 'Parameter Search Surface', 'Position', [150, 150, 800, 600]);
[B_mesh, A_mesh] = meshgrid(B_vals, A_vals);
surf(B_mesh, A_mesh, SNR_matrix, 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('Parameter B');
ylabel('Parameter A');
zlabel('SNR (dB)');
title(sprintf('SNR Surface (Max: %.2f dB at A=%.2f, B=%.2f)', best_SNR, best_A, best_B));
view(-45, 45); % Set a nice viewing angle

% Add a red star at the maximum point
hold on;
plot3(best_B, best_A, best_SNR + 2, 'p', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
hold off;

% =========================================================================
% LOCAL FUNCTIONS (Paste your previously made functions here)
% =========================================================================
% IMPORTANT: Make sure to paste the updated versions of:
% 1. chebyshev_hybrid_interpolator
% 2. S_operator
% 3. P_operator
% 4. ideal_lowpass_filter
% below this line if you are running this as a single standalone script.