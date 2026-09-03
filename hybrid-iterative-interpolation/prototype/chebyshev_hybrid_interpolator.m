function [x_final] = chebyshev_hybrid_interpolator(x_sampled, L, N, iters, A, B)
% CHEBYSHEV_HYBRID_INTERPOLATOR Reconstructs an analog signal from discrete 
% samples using the mathematically corrected Chebyshev acceleration.

    % 1. Calculate parameters
    rho = (B - A) / (B + A);
    alpha = 2 / (A + B);

    % 2. Pre-calculate Chebyshev Lambda Sequence
    lam = zeros(1, iters);
    lam(1) = 2; 
    for n = 2:iters
        lam(n) = 1 / (1 - (rho^2 / 4) * lam(n-1));
    end

    % 3. The Constant Driving Term (D)
    % This is the base modular estimation that drives the whole system
    x_sh_initial = S_operator(x_sampled, L);
    D = P_operator(x_sh_initial, L, N); 

    % 4. Iteration 0 (Initial Guess)
    y_0 = D;

    % 5. Iteration 1 (Standard Gradient Descent Step)
    PS_y0 = P_operator(S_operator(y_0, L), L, N);
    y_1 = y_0 + alpha * (D - PS_y0);

    % Initialize memory variables for the loop
    y_prev2 = y_0; % Represents y_{n-2}
    y_prev1 = y_1; % Represents y_{n-1}

    % 6. The Main Chebyshev Iteration Loop
    for n = 2:iters
        
        % Apply operators to the previous guess
        PS_y_prev1 = P_operator(S_operator(y_prev1, L), L, N);
        
        % Execute the TRUE Chebyshev formula
        % Notice how cleanly it calculates the error: (D - PS_y_prev1)
        step = y_prev1 + alpha * (D - PS_y_prev1) - y_prev2;
        y_n = lam(n) * step + y_prev2;
        
        % Update variables for the next loop iteration
        y_prev2 = y_prev1;
        y_prev1 = y_n;
        
    end

    % 7. Output 
    if iters == 0
        x_final = y_0;
    elseif iters == 1
        x_final = y_1;
    else
        x_final = y_n;
    end
end