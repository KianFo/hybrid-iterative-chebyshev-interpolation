function xHold = zoh_operator(x, L)
%ZOH_OPERATOR  Sample-and-hold (zero-order-hold) sampling operator  S.
%
%   xHold = ZOH_OPERATOR(x, L) keeps the Nyquist-rate samples of x -- one
%   sample every L points, i.e. the ideal samples x(nT) -- and holds each of
%   them constant until the next sample.  The result is the familiar
%   "staircase" that a plain D/A converter would produce.
%
%   A raw held staircase is delayed by half a hold interval, because the value
%   x(nT) is held over the whole interval [nT, (n+1)T) instead of being
%   centred on nT.  This T/2 group delay is exactly the ZOH phase distortion.
%   We cancel it with a circular shift of L/2 samples so that S behaves as a
%   zero-phase operator and lines up with the ideal lowpass filter that the
%   band-limiting operator P applies afterwards.
%
%   See also MODULAR_OPERATOR, HYBRID_RECONSTRUCT.

    x = x(:).';                          % work with a row vector
    samples = x(1:L:end);                % ideal Nyquist samples  x(nT)
    xHold   = repelem(samples, L);       % hold each sample for L points
    xHold   = circshift(xHold, -floor(L/2));   % remove the ZOH T/2 delay
    xHold   = xHold(1:numel(x));         % trim to the original length
end
