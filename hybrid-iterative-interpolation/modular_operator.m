function xProj = modular_operator(x, L, M)
%MODULAR_OPERATOR  Band-limiting operator  P  (the "modular" interpolator).
%
%   xProj = MODULAR_OPERATOR(x, L, M) compensates the sample-and-hold
%   distortion using the modular method of Marvasti (1985): the signal is
%   first mixed with the first M harmonics of the sampling impulse train,
%
%       g(t) = 1 + 2 * sum_{s=1}^{M} cos(2*pi*s*t / L) ,
%
%   and the product is then passed through an ideal brick-wall lowpass filter
%   onto the baseband.  g(t) is the truncated Fourier series of the Dirac
%   comb of period L; mixing with it folds the spectral replicas created by
%   the hold back onto the baseband so that they cancel the hold distortion.
%
%   With M = 1 the composite operator P*S is already very close to the
%   identity on the baseband (the paper reports ||I - P*S|| ~ 0.06), which is
%   what gives the iterative scheme its fast convergence.  Setting M = 0 turns
%   P into a plain ideal lowpass filter (the classical, non-hybrid case).
%
%   See also ZOH_OPERATOR, HYBRID_RECONSTRUCT.

    x = x(:).';
    n = numel(x);
    t = 0:n-1;

    % --- modular mixing: multiply by the truncated sampling comb ------------
    g = ones(1, n);
    for s = 1:M
        g = g + 2*cos(2*pi*s*t / L);
    end
    mixed = x .* g;

    % --- ideal lowpass filter onto the baseband ----------------------------
    % Keep only |frequency| < Nyquist.  The zeroed band is chosen so that the
    % surviving spectrum stays conjugate-symmetric; otherwise a real signal
    % would lose its Nyquist-bin energy on every pass and the iteration would
    % hit an artificial SNR floor.
    K = floor(n/(2*L));
    X = fft(mixed);
    X(K+1 : n-K+1) = 0;                  % symmetric about DC in MATLAB indexing
    xProj = real(ifft(X));
end
