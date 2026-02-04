function D = dmd_projected(T_scores, fs, r, opts)
    Z  = T_scores;
    X1 = Z(1:end-1,:).';
    X2 = Z(2:end,  :).';

    % Truncated SVD of X1
    [U,S,V] = svd(X1, 'econ');
    sing = diag(S);
    tol  = max(size(X1)) * eps(max(sing));
    r_num = sum(sing > tol);
    r = min([r, r_num]);

    Ur = U(:,1:r); Sr = S(1:r,1:r); Vr = V(:,1:r);

    % Regularized inverse to avoid warnings
    lambda = max(opts.ridge, 1e-8);
    invSr = diag(1./(diag(Sr) + lambda));

    % Low-order operator
    Atil = Ur' * X2 * Vr * invSr;

    % Eigen-decomposition
    [W, Dlam] = eig(Atil);
    lambda_d  = diag(Dlam);
    dt = 1/fs;

    % Guard log on tiny/zero eigenvalues
    lambda_d(abs(lambda_d) < 1e-12) = 1e-12;
    omega   = log(lambda_d)/dt;
    freq_hz = abs(imag(omega))/(2*pi);

    % Mode amplitudes (pinv for safety)
    a0 = Ur' * Z(1,:).';
    b  = pinv(W, opts.rcond) * a0;

    % Time dynamics
    t = (0:size(Z,1)-1)' * dt;
    Phi_t = exp(t * omega.');           % [T x r]
    Phi_t = Phi_t .* (ones(size(Phi_t,1),1) * b.');  % scale by b

    D = struct('modes_time', Phi_t, 'lambda', lambda_d, 'omega', omega, ...
               'freq_hz', freq_hz, 'amplitudes', abs(b));
end


