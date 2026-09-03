function [xRec, snrTrace] = hybrid_reconstruct(x, p)
%HYBRID_RECONSTRUCT  Reconstruct a bandlimited signal from its Nyquist-rate
%   sample-and-hold version using the hybrid iterative interpolation method of
%   A. ParandehGheibi, M. A. Akhaee and F. Marvasti, "Design of Interpolation
%   Functions Using Iterative Methods", EUSIPCO 2006.
%
%   xRec = HYBRID_RECONSTRUCT(x, p) returns the reconstructed signal.  The
%   parameter struct p has the following fields:
%
%       p.L        hold factor / oversampling ratio  (the period T in samples)
%       p.modules  number of cosine modules M in the P operator
%                  (0 -> classical iterative method, >=1 -> hybrid method)
%       p.iters    number of correction iterations
%       p.method   'iterative'   : relaxation update (Eq. 4 of the paper)
%                  'accelerated' : Chebyshev acceleration (Eq. 26)
%       p.relax    relaxation factor lambda      (used by 'iterative')
%       p.A, p.B   frame bounds                  (used by 'accelerated')
%
%   [xRec, snrTrace] = HYBRID_RECONSTRUCT(x, p) additionally returns, when
%   p.ref is supplied, the SNR (dB) after every iteration so convergence can
%   be plotted.  p.ref is the ground-truth signal and is optional.
%
%   IDEA.  Sampling then band-limiting the true signal gives the "measurement"
%   D = P*S*x.  Because P*S is close to the identity on the baseband, every
%   method below refines an estimate y until P*S*y matches D, at which point
%   y has converged to the original x.  The Chebyshev variant reaches the same
%   fixed point in far fewer iterations by reusing the two previous estimates.
%
%   See also ZOH_OPERATOR, MODULAR_OPERATOR, SNR_METRIC.

    % Operator handles: S = sample-and-hold, P = modular band-limiter.
    S  = @(v) zoh_operator(v, p.L);
    P  = @(v) modular_operator(v, p.L, p.modules);
    PS = @(v) P(S(v));

    wantTrace = (nargout > 1) && isfield(p, 'ref');
    snrTrace  = [];

    D = PS(x);                 % initial modular estimate  x_hat = P*S*x

    switch lower(p.method)
        %------------------------------------------------------------------
        case 'iterative'       % relaxation:  y <- y + lambda (D - P S y)
        %------------------------------------------------------------------
            y = D;
            if wantTrace, snrTrace(1) = snr_metric(p.ref, y); end
            for k = 1:p.iters
                y = y + p.relax * (D - PS(y));
                if wantTrace, snrTrace(k+1) = snr_metric(p.ref, y); end
            end
            xRec = y;

        %------------------------------------------------------------------
        case 'accelerated'     % Chebyshev three-term acceleration
        %------------------------------------------------------------------
            rho   = (p.B - p.A) / (p.B + p.A);
            alpha = 2 / (p.A + p.B);

            yTwoBack = D;                          % y_0
            if wantTrace, snrTrace(1) = snr_metric(p.ref, yTwoBack); end
            if p.iters == 0, xRec = yTwoBack; return; end

            yOneBack = D + alpha * (D - PS(D));    % y_1  (single relaxed step)
            if wantTrace, snrTrace(2) = snr_metric(p.ref, yOneBack); end

            lambda = 2;                            % Chebyshev weight, lambda_1
            y = yOneBack;
            for n = 2:p.iters
                lambda = 1 / (1 - (rho^2/4) * lambda);          % Eq. 27
                y = lambda * (yOneBack + alpha*(D - PS(yOneBack)) - yTwoBack) ...
                    + yTwoBack;                                  % Eq. 26
                yTwoBack = yOneBack;
                yOneBack = y;
                if wantTrace, snrTrace(n+1) = snr_metric(p.ref, y); end
            end
            xRec = y;

        otherwise
            error('hybrid_reconstruct:badMethod', ...
                  'Unknown method "%s" (use ''iterative'' or ''accelerated'').', ...
                  p.method);
    end
end
