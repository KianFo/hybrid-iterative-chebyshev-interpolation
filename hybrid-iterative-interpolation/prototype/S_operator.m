function [x_sh] = S_operator(x_analog, L)
% S_OPERATOR Simulates a ZERO-PHASE Sample-and-Hold process.

    % Step 1: Sampling
    x_discrete = x_analog(1:L:end);
    
    % Step 2: Holding (Causal)
    x_sh = repelem(x_discrete, L);
    
    % Step 3: Phase Correction (The Magic Fix)
    % Shift the signal left by half a period to center the pulses
    x_sh = circshift(x_sh, -floor(L/2));
    
    % Step 4: Dimension Matching
    x_sh = x_sh(1:length(x_analog));
end