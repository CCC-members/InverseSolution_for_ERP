% ------- legacy helpers kept for multi-trial fallback --------
function SSA = mssa_original(X, L, r)
    [T,S,K] = size(X);
    Xcat = reshape(permute(X,[1 3 2]), T*K, S);
    N = size(Xcat,1);
    L = min(L, N-1);
    Kwin = N - L + 1;
    H = zeros(L*S, Kwin);
    for s=1:S
        H((s-1)*L+(1:L), :) = hankel(Xcat(1:L,s), Xcat(L:end,s));
    end
    [U,Sv,V] = svd(H, 'econ');
    r = min([r, size(U,2)]);
    Ur = U(:,1:r); Sr = diag(Sv(1:r,1:r)); Vr = V(:,1:r);
    U_time = zeros(N, r);
    for i=1:r
        H_i = Ur(:,i) * (Sr(i) * Vr(i))';
        U_time(:,i) = hankel_antiavg(H_i, L);
    end
    V_space = zeros(S, r);
    for i=1:r
        V_space(:,i) = Xcat \ U_time(:,i);
    end
    SSA = struct('U', U_time, 'V', V_space, 's', Sr);
end

