function [x_p] = P_operator(x_in, L, N)
% P_OPERATOR Applies the Modular Method: cosine mixing followed by lowpass filtering.
% Inputs:
%   x_in - The input signal (1D array, usually from S_operator)
%   L    - The hold factor (Ratio of analog fs to sampled fs, representing T)
%   N    - Number of cosine modules (Recommended: 1)
% Output:
%   x_p  - The filtered, interpolated signal

    % Ensure input is a row vector for consistent matrix math
    x_in = x_in(:).'; 
    sig_len = length(x_in);
    
    % Step 1: Generate the discrete time index vector
    t = 0 : (sig_len - 1);
    
    % Step 2: Calculate the Modular Harmonics Multiplier
    modulator = ones(1, sig_len);
    for s = 1:N
        modulator = modulator + 2 * cos(2 * pi * s * t / L);
    end
    
    % Step 3: Mix (multiply) the input signal with the modulator
    x_mixed = x_in .* modulator;
    
    % Step 4: Pass the mixed signal through the standalone lowpass filter
    x_p = ideal_lowpass_filter(x_mixed, L);
    
end