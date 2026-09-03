function snr = snr_metric(ref, est, marginFrac)
%SNR_METRIC  Reconstruction signal-to-noise ratio in dB.
%
%   snr = SNR_METRIC(ref, est) compares the estimate est against the reference
%   ref and returns 10*log10( signal power / error power ).  The first and
%   last 10% of the record are ignored so that filter/transient edge effects
%   do not corrupt the measurement (this matches the paper's methodology).
%
%   snr = SNR_METRIC(ref, est, marginFrac) uses a custom edge margin fraction.

    if nargin < 3, marginFrac = 0.10; end
    ref = ref(:).';
    est = est(:).';

    n   = numel(ref);
    m   = floor(marginFrac * n);
    idx = (m+1):(n-m);

    err = ref(idx) - est(idx);
    snr = 10 * log10(sum(ref(idx).^2) / sum(err.^2));
end
