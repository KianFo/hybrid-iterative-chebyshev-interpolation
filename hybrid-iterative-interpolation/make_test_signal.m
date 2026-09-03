function x = make_test_signal(kind, N, L, seed)
%MAKE_TEST_SIGNAL  Build a length-N, unit-normalised test signal that is
%   bandlimited to the baseband defined by the oversampling ratio L (so that
%   sampling one point in every L is exactly Nyquist-rate sampling).
%
%   x = MAKE_TEST_SIGNAL(kind, N, L, seed) supports:
%       'bandlimited' - lowpass-filtered pseudo-random noise (paper's signal)
%       'tones'       - a sum of three sinusoids below the Nyquist frequency
%       'chirp'       - a linear frequency sweep staying below Nyquist
%
%   seed seeds the random generator for the 'bandlimited' case (default 0).

    if nargin < 4, seed = 0; end
    t     = (0:N-1) / N;          % normalised time, one "second" of data
    fNyq  = N / (2*L);            % Nyquist frequency in cycles over the record

    switch lower(kind)
        case 'bandlimited'
            rng(seed);
            X = fft(randn(1, N));
            K = floor(N/(2*L));
            X(K+1 : N-K+1) = 0;          % same symmetric baseband as the LPF
            x = real(ifft(X));

        case 'tones'
            % three tones placed well inside the band (0.15, 0.4, 0.8 * Nyq)
            x = 1.00*sin(2*pi*round(0.15*fNyq)*t) ...
              + 0.50*cos(2*pi*round(0.40*fNyq)*t) ...
              + 0.25*sin(2*pi*round(0.80*fNyq)*t);

        case 'chirp'
            x = chirp(t, 1, t(end), floor(0.95*fNyq));

        otherwise
            error('make_test_signal:badKind', 'Unknown signal kind "%s".', kind);
    end

    x = x / max(abs(x));         % normalise to unit peak amplitude
end
