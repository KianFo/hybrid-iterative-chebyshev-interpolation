function [x_final] = classical_hybrid_interpolator(x_sampled, L, N, iters, lambda)
    % The standard iteration method without Chebyshev (Equation 4)
    x_sh_initial = S_operator(x_sampled, L);
    D = P_operator(x_sh_initial, L, N); 
    x_k = D;
    
    for k = 1:iters
        PS_xk = P_operator(S_operator(x_k, L), L, N);
        x_k = lambda * D + x_k - lambda * PS_xk;
    end
    x_final = x_k;
end