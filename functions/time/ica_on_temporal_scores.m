function [IC_time, IC_space] = ica_on_temporal_scores(T_scores, S_loads, nIC, opts)
    [T, r] = size(T_scores);
    % whiten with numerical rank guard
    C = (T_scores.'*T_scores)/(T-1);
    [E,D] = eig((C+C')/2); 
    [d,ord]=sort(diag(D),'descend'); E=E(:,ord); d=max(d,1e-12);
    r_numerical = sum(d > max(size(C))*eps(max(d)));
    r_use = min([nIC, r, r_numerical]);

    Zwh = T_scores * E(:,1:r_use) * diag(1./sqrt(d(1:r_use))); % [T x r_use]
    % FastICA or fixed-point
    if exist('fastica','file') == 2
        [Sic, ~, W] = fastica(Zwh.', 'verbose','off', 'displayMode','off', 'numOfIC', r_use);
        IC_time = Sic.';         % [T x r_use]
    else
        W = randn(r_use, r_use);
        for it=1:200
            G = tanh(Zwh*W.'); dG = 1 - G.^2;
            Wn = (G.'*Zwh)/T - diag(mean(dG,1))*W;
            [Uu,Us,~] = svd(Wn, 'econ');
            Wn = Uu*diag(1./sqrt(diag(Us)))*Uu'*Wn;
            if norm(abs(diag(Wn*W'))-1,'fro') < 1e-6, W = Wn; break; end
            W = Wn;
        end
        IC_time = Zwh*W.';       % [T x r_use]
        W = W.';                 % make W consistent below
    end
    % map to sources via ridge/pinv (temporal->IC transform)
    M_temporal_to_IC = E(:,1:r_use) * diag(1./sqrt(d(1:r_use))) * safe_inv(W, opts.rcond);
    IC_space = S_loads * M_temporal_to_IC;    % [S x r_use]
end


function Winv = safe_inv(W, rcondv)
% Pseudo-inverse with rcond guard (for non-square or ill-conditioned)
    if nargin < 2 || isempty(rcondv), rcondv = 1e-10; end
    if size(W,1) == size(W,2)
        if rcond(W) < rcondv
            Winv = pinv(W, rcondv);
        else
            Winv = inv(W);
        end
    else
        Winv = pinv(W, rcondv);
    end
end