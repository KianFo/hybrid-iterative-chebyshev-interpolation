% =========================================================================
% MAIN SCRIPT: Chebyshev Hybrid Interpolator Testing
% =========================================================================
clc; clear; close all;

%% 1. System Parameters
N_modules = 1;              % Number of cosine modules
iters = 2;                  % Number of Chebyshev iterations
A = 1.04;                   % Optimum Lower Bound
B = 1.06;                   % Optimum Upper Bound

sig_len = 4096;             % Total length of the analog signal
nyq_samples = 64;           % Number of discrete samples (Nyquist rate)
L = sig_len / nyq_samples;  % Hold Factor (T)

%% 2. Signal Generation Selection
% Change this variable (1, 2, or 3) to test different signals!
signal_type = 2; 

t = linspace(0, 1, sig_len); % Time vector for 1 second

switch signal_type
    case 1
        % Signal 1: Bandlimited Random Signal (Paper's Method)
        rng(42); 
        x_rand = randn(1, sig_len);
        X_f = fft(x_rand);
        cutoff_idx = floor(sig_len / (2 * L)); 
        X_f(cutoff_idx+1 : end-cutoff_idx) = 0; 
        x_analog = real(ifft(X_f));
        test_name = 'Bandlimited Random Signal';
        
    case 2
        % Signal 2: Sum of Sinusoids (e.g., a simple audio chord)
        % Frequencies must be below the Nyquist limit (Nyquist is nyq_samples/2 = 32 Hz)
        f1 = 5;  % 5 Hz
        f2 = 12; % 12 Hz
        f3 = 25; % 25 Hz (Close to Nyquist)
        x_analog = sin(2*pi*f1*t) + 0.5*cos(2*pi*f2*t) + 0.25*sin(2*pi*f3*t);
        test_name = 'Sum of Sinusoids';
        
    case 3
        % Signal 3: Chirp Signal (Swept Frequency)
        % Sweeps from 1 Hz to 30 Hz
        x_analog = chirp(t, 1, 1, 30);
        test_name = 'Chirp Signal (1 to 30 Hz)';
end

% Normalize the chosen signal
x_analog = x_analog / max(abs(x_analog)); 

%% 3. Run the Interpolator
fprintf('Running Hybrid Interpolation on: %s\n', test_name);
x_recon = chebyshev_hybrid_interpolator(x_analog, L, N_modules, iters, A, B);

%% 4. Calculate SNR
ignore_len = floor(0.10 * sig_len);
valid_idx = (ignore_len + 1) : (sig_len - ignore_len);

x_orig_valid = x_analog(valid_idx);
x_recon_valid = x_recon(valid_idx);

signal_power = sum(x_orig_valid.^2);
error_power = sum((x_orig_valid - x_recon_valid).^2);
SNR_dB = 10 * log10(signal_power / error_power);

fprintf('Achieved SNR: %.2f dB\n\n', SNR_dB);

%% 5. Plotting
figure('Name', sprintf('Test: %s', test_name), 'Position', [100, 100, 900, 500]);

% Generate the S&H signal just for plotting comparison
x_sh = S_operator(x_analog, L);

plot(t, x_analog, 'b', 'LineWidth', 1.5); hold on;
plot(t, x_sh, 'r', 'LineWidth', 1);
plot(t, x_recon, 'k--', 'LineWidth', 1.5);

title(sprintf('%s | Reconstruction SNR = %.2f dB', test_name, SNR_dB));
legend('Ideal Analog Signal', 'Raw Sample-and-Hold', 'Chebyshev Reconstructed');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0.1, 0.3]); % Zoom in to see the details clearly
