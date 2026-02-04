function [J, W_slor, lambda_opt, gcv_vals, lam_grid] = sloreta_auto(L, V, C)
% sLORETA with automatic lambda selection by GCV
% Inputs
%   L : [nChan x nSrc] leadfield (fixed orientation)
%   V : [nChan x nTime] data matrix
%   C : [nChan x nChan] noise covariance (optional, defaults to I)
% Outputs
%   J          : [nSrc x nTime] sLORETA estimates
%   W_slor     : [nSrc x nChan] sLORETA inverse operator (maps ORIGINAL V -> J)
%   lambda_opt : chosen regularization parameter (scalar)
%   gcv_vals   : GCV values over lam_grid (for diagnostics)
%   lam_grid   : lambda grid used

    if nargin < 3 || isempty(C), C = eye(size(L,1)); end
    nChan = size(L,1);

    % --- Whitening (correct orientation) ---
    % C ~ U*U', lower chol; whitening is U^{-1} * X  -> MATLAB: U \ X
    U = chol(C + 1e-12*eye(nChan), 'lower');
    Lw = U \ L;         % [nChan x nSrc] correccion aqui
    Vw = U \ V;         % [nChan x nTime] correccion aqui

    % --- SVD of whitened leadfield ---
    [U_s, svals, ~] = svd(Lw, 'econ');
    s = diag(svals);    % singular values, length r = rank(Lw)

    % --- GCV scan on a log grid ---
    lam_grid = logspace(-6, 2, 60);
    I_ch = eye(nChan);
    gcv_vals = zeros(size(lam_grid));

    for k = 1:numel(lam_grid)
        lam = lam_grid(k);
        % Hat matrix in whitened sensor space: H = U_s * diag(s.^2./(s.^2+lam)) * U_s'
        f = (s.^2) ./ (s.^2 + lam);
        H = U_s * (diag(f)) * U_s';
        resid = Vw - H*Vw;
        num = sum(resid(:).^2);
        den = (trace(I_ch - H))^2;
        gcv_vals(k) = num / max(den, 1e-30);
    end

    [~, kopt] = min(gcv_vals);
    lambda_opt = lam_grid(kopt);

    % --- Tikhonov inverse in whitened space with optimal lambda ---
    LtL = Lw.' * Lw;
    W   = (LtL + lambda_opt * eye(size(LtL))) \ (Lw.');  % [nSrc x nChan]

    % --- sLORETA row-wise standardization ---
    % Resolution matrix in source space (whitened): R = W*Lw
    R = W * Lw;
    d = sqrt(max(diag(R), 1e-15));   % avoid zero/neg roundoff
    W_std = W ./ d;                  % divide each row i by sqrt(R_ii)

    % --- Unwhiten operator so it acts on ORIGINAL V ---
    % We currently have J = W_std * (U \ V). To get J = (W_std / U) * V:
    W_slor = W_std / U;

    % --- Final estimates on original data ---
    J = W_slor * V;
end
