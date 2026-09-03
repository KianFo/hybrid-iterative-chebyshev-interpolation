function [x_filtered] = ideal_lowpass_filter(x_in, L)
% IDEAL_LOWPASS_FILTER Applies an ideal brick-wall lowpass filter
% with PERFECT symmetry to prevent imaginary frequency leakage.

    sig_len = length(x_in);
    X_f = fft(x_in);
    cutoff_idx = floor(sig_len / (2 * L));
    
    % SYMMETRIC ZEROING:
    % Keeps exactly DC + (cutoff_idx) positive freqs and matching negative freqs.
    % Zeroes out everything in between cleanly.
    X_f((cutoff_idx + 1) : (sig_len - cutoff_idx + 1)) = 0;
    
    x_filtered = real(ifft(X_f));
end