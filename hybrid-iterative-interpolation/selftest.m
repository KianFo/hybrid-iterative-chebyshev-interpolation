% SELFTEST  Headless numerical check of the implementation (no figures).
function selftest()
    N = 4096; Nyq = 64; L = N/Nyq;
    x = make_test_signal('bandlimited', N, L, 42);

    fprintf('Raw S&H            : %6.2f dB\n', snr_metric(x, zoh_operator(x,L)));

    p = struct('L',L,'modules',1,'method','accelerated','A',1.03,'B',1.05);
    for it = 0:4
        p.iters = it;
        fprintf('Accel hybrid it=%d : %7.2f dB\n', it, snr_metric(x, hybrid_reconstruct(x,p)));
    end

    q = struct('L',L,'modules',0,'method','iterative','relax',0.94,'iters',2);
    fprintf('Classical it=2     : %6.2f dB\n', snr_metric(x, hybrid_reconstruct(x,q)));
    disp('SELFTEST OK');
end
