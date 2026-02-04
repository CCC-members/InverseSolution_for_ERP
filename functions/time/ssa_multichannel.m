function [Y, singvals] = ssa_multichannel(T_scores, L, rssa, opts)
    [T, r0] = size(T_scores);
    L = min(L, T-1);
    K = T - L + 1;
    H = zeros(L*r0, K);
    for i=1:r0
        xi = T_scores(:,i);
        H((i-1)*L+(1:L), :) = hankel(xi(1:L), xi(L:end));
    end
    [U,S,~] = svd(H, 'econ');
    sing = diag(S);
    tol = max(size(H)) * eps(max(sing));
    r_num = sum(sing > tol);
    rssa = min([rssa, r_num]);
    U = U(:,1:rssa); singvals = sing(1:rssa);
    Y = zeros(T, rssa);
    for k=1:rssa
        Hk = U(:,k) * (U(:,k)' * H);    % rank-1 approximation
        Y(:,k) = hankel_antiavg(Hk, L);
    end
end


function x = hankel_antiavg(H, L)
    [~,K] = size(H); N = L + K - 1;
    x = zeros(N,1); c = zeros(N,1);
    for i=1:L
        for j=1:K
            idx = i + j - 1;
            x(idx) = x(idx) + H(i,j);
            c(idx) = c(idx) + 1;
        end
    end
    x = x ./ max(c,1);
end
