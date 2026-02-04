function m = bandpower_map(Pxx, f, rngHz)
    idx = f >= rngHz(1) & f <= rngHz(2);
    m = trapz(f(idx), Pxx(:,idx), 2);
    if isrow(m), m = m.'; end % enforce [S x 1]
end

