function out = time_series_analysis(X, fs, Cortex, opts)
% PSD (0.1–50 Hz), ICA, MSSA, DMD for inverse solution.
% X: [T x S] or [T x S x K]  (single-trial -> K==1)
% fs: Hz, Cortex: Brainstorm tess struct
% opts: nICA, ssa_L, ssa_r, dmd_r, r_pod, show_plots, bands

    if nargin < 4, opts = struct; end
    if ndims(X) == 2, X = reshape(X, size(X,1), size(X,2), 1); end
    [T,S,K] = size(X);

    % ---------- defaults ----------
    def.ssa_L = max(10, round(2*fs));
    def.ssa_r = 20;
    def.dmd_r = min(40, T-2);
    def.r_pod = min(80, T-2);
    def.nICA  = min(30, def.r_pod);
    def.show_plots = true;
    def.bands = struct('slow',[0.1 0.9],'delta',[1 4],'theta',[4 7], ...
                       'alpha',[8 12],'beta',[13 30],'gammaL',[30 50]);
    fns = fieldnames(def);
    for i=1:numel(fns)
        if ~isfield(opts, fns{i}), opts.(fns{i}) = def.(fns{i}); end
    end

    % ---------- detrend/center ----------
    X = detrend(X, 0);
    X = bsxfun(@minus, X, mean(X,1));

    % ---------- 1) PSD via Welch ----------
    nfft = max(2^nextpow2(round(2*fs)), 512);
    wlen = round(2*fs);                         % 2-s Hamming
    nover = round(0.5*wlen);
    [Pxx, f] = welch_per_source_singletrial(X(:,:,1), fs, wlen, nover, nfft); % [S x F]
    bandNames = fieldnames(opts.bands);
    BandMaps = struct();
    for i=1:numel(bandNames)
        rngHz = opts.bands.(bandNames{i});
        BandMaps.(bandNames{i}) = bandpower_map(Pxx, f, rngHz);  % [S x 1]
    end

    % =================== FAST PATH FOR ONE TRIAL ===================
    if K == 1
        % ---- Temporal POD (economy SVD) ----
        [U_t, Sigma, V_s] = svd(X(:,:,1), 'econ');  % X ≈ U_t * Sigma * V_s'
        sing = diag(Sigma);
        r_full = numel(sing);
        r = min([opts.r_pod, r_full]);
        U_t = U_t(:,1:r);
        V_s = V_s(:,1:r);
        Sigma_r = Sigma(1:r,1:r);
        T_scores = U_t * Sigma_r;   % [T x r]
        S_loads  = V_s;             % [S x r]

        % ---- ICA on temporal PCs, map back to space ----
        [IC_time, IC_space] = ica_on_temporal_scores(T_scores, S_loads, opts.nICA);

        % ---- MSSA on temporal PCs, map back to space ----
        L = min(opts.ssa_L, T-1);
        [SSA_time, SSA_sing] = ssa_multichannel(T_scores, L, min(opts.ssa_r, r));
        SSA_space = S_loads * (T_scores \ SSA_time);

        % ---- Projected DMD (on temporal PCs), map back ----
        DMD = dmd_projected(T_scores, fs, min(opts.dmd_r, r));
        keep = DMD.freq_hz >= 0.1 & DMD.freq_hz <= 50 & isfinite(DMD.freq_hz);
        DMD_time  = DMD.modes_time(:, keep);
        DMD_modes = S_loads * (T_scores \ DMD_time);
        DMD_freqs = DMD.freq_hz(keep);
        DMD_amps  = DMD.amplitudes(keep);
    else
        % ---------- 2) ICA (multi-trial fallback) ----------
        Xica = reshape(permute(X, [1 3 2]), T*K, S);
        Xica = bsxfun(@minus, Xica, mean(Xica,1));
        [~, A_ica, S_ica] = run_ica(Xica, min(opts.nICA, T-2));
        IC_time = reshape(S_ica, T, K, []);
        IC_space = A_ica;

        % ---------- 3) MSSA ----------
        L = min(opts.ssa_L, T-1);
        SSA = mssa_original(X, L, opts.ssa_r);
        SSA_time = reshape(SSA.U, T, []);
        SSA_space = SSA.V;
        SSA_sing = SSA.s;

        % ---------- 4) DMD ----------
        DMD = dmd_concat_trials(X, fs, opts.dmd_r);
        keep = DMD.freq_hz >= 0.1 & DMD.freq_hz <= 50 & imag(DMD.omega)==0;
        DMD_modes = DMD.modes(:, keep);
        DMD_freqs = DMD.freq_hz(keep);
        DMD_amps  = DMD.amplitudes(keep);
        DMD_time  = DMD.reconstruct_time(:, keep);
    end

    % ---------- output ----------
    out = struct('f', f, 'Pxx', Pxx, 'BandMaps', BandMaps, ...
                 'IC_space', IC_space, 'IC_time', squeeze(IC_time), ...
                 'SSA_space', SSA_space, 'SSA_time', SSA_time, 'SSA_sing', SSA_sing, ...
                 'DMD_modes', DMD_modes, 'DMD_freqs', DMD_freqs, 'DMD_amps', DMD_amps, ...
                 'DMD_time', DMD_time);

    % ---------- plotting (your renderer) ----------
    if opts.show_plots && nargin >= 3 && ~isempty(Cortex)
        cmap = load("tools/mycolormap_brain_basic_conn.mat");
        nV = size(Cortex.Vertices,1);
        fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);

        % Alpha & Delta bandpower maps
        PlotSourceTopography_safe(Cortex, fixmap(BandMaps.alpha), cmap, 'Alpha band power (8–12 Hz)');
        PlotSourceTopography_safe(Cortex, fixmap(BandMaps.delta), cmap, 'Delta band power (1–4 Hz)');

        % ICA #1
        if ~isempty(IC_space),  PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(IC_space(:,1))), cmap, 'ICA #1 spatial'); end
        % SSA #1
        if ~isempty(SSA_space), PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(SSA_space(:,1))), cmap, 'SSA #1 spatial'); end
        % DMD top amplitude
        if ~isempty(DMD_modes)
            [~, iTop] = max(DMD_amps);
            PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(real(DMD_modes(:,iTop)))), cmap, ...
                sprintf('DMD mode @ %.2f Hz (max amp)', DMD_freqs(iTop)));
        end

        % Temporal summary figure
        figure('Name','Temporal components');
        subplot(3,1,1); plot(out.f, 10*log10(mean(out.Pxx,1))); xlim([0.1 50]); grid on
        xlabel('Hz'); ylabel('dB'); title('Mean PSD (across sources)');
        subplot(3,1,2); if ~isempty(out.IC_time), plot(out.IC_time(:,1)); title('ICA #1'); grid on; end
        subplot(3,1,3); if ~isempty(out.SSA_time), plot(out.SSA_time(:,1)); title('SSA #1 temporal PC'); grid on; end
    end
end

% ===================== helpers =====================

function [Pxx, f] = welch_per_source_singletrial(X, fs, wlen, nover, nfft)
% X: [T x S], returns Pxx as [S x F]
    [P, f] = pwelch(X, hamming(wlen), nover, nfft, fs, 'psd'); % P: [F x S]
    Pxx = P.'; % [S x F]
end

function m = bandpower_map(Pxx, f, rngHz)
    idx = f >= rngHz(1) & f <= rngHz(2);
    m = trapz(f(idx), Pxx(:,idx), 2);
    if isrow(m), m = m.'; end % enforce [S x 1]
end

function [IC_time, IC_space] = ica_on_temporal_scores(T_scores, S_loads, nIC)
    [T, r] = size(T_scores);
    r_use = min(nIC, r);
    % whiten
    C = (T_scores.'*T_scores)/(T-1);
    [E,D] = eig((C+C')/2); [d,ord]=sort(diag(D),'descend'); E=E(:,ord); d=max(d,1e-12);
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
    end
    % map to sources
    M_temporal_to_IC = E(:,1:r_use) * diag(1./sqrt(d(1:r_use))) * (W'\eye(r_use));
    IC_space = S_loads * M_temporal_to_IC;    % [S x r_use]
end

function [Y, singvals] = ssa_multichannel(T_scores, L, rssa)
    [T, r0] = size(T_scores);
    L = min(L, T-1);
    K = T - L + 1;
    H = zeros(L*r0, K);
    for i=1:r0
        xi = T_scores(:,i);
        H((i-1)*L+(1:L), :) = hankel(xi(1:L), xi(L:end));
    end
    [U,S,~] = svd(H, 'econ');
    rssa = min(rssa, size(U,2));
    U = U(:,1:rssa); singvals = diag(S(1:rssa,1:rssa));
    Y = zeros(T, rssa);
    for k=1:rssa
        Hk = U(:,k) * (U(:,k)' * H);
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

function D = dmd_projected(T_scores, fs, r)
    Z = T_scores;
    X1 = Z(1:end-1,:).';
    X2 = Z(2:end,  :).';
    [U,S,V] = svd(X1, 'econ');
    r = min([r, size(U,2)]);
    Ur = U(:,1:r); Sr = S(1:r,1:r); Vr = V(:,1:r);
    Atil = Ur' * X2 * Vr / Sr;
    [W, Dlam] = eig(Atil);
    lambda = diag(Dlam);
    dt = 1/fs;
    omega = log(lambda)/dt;
    freq_hz = abs(imag(omega))/(2*pi);
    b = W \ (Ur' * Z(1,:).');
    t = (0:size(Z,1)-1)' * dt;
    Phi_t = zeros(size(Z,1), r);
    for i=1:r
        Phi_t(:,i) = exp(omega(i)*t) * b(i);
    end
    D = struct('modes_time', Phi_t, 'lambda', lambda, 'omega', omega, ...
               'freq_hz', freq_hz, 'amplitudes', abs(b));
end

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
        H_i = Ur(:,i) * (Sr(i) * Vr(:,i))';
        U_time(:,i) = hankel_antiavg(H_i, L);
    end
    V_space = zeros(S, r);
    for i=1:r
        V_space(:,i) = Xcat \ U_time(:,i);
    end
    SSA = struct('U', U_time, 'V', V_space, 's', Sr);
end

function D = dmd_concat_trials(~, ~, ~) %#ok<STOUT,INUSD>
    error('dmd_concat_trials not used in single-trial fast path.');
end

% ===================== PLOTTING (your renderer) =====================

function PlotSourceTopography_safe(Cortex,J,cmap_struct,ttl)
% Wrapper around user's PlotSourceTopography with guards and a title.
    if ~isfield(Cortex,'Faces') || ~isfield(Cortex,'Vertices')
        warning('Cortex lacks Faces/Vertices; skipping PlotSourceTopography.');
        return;
    end
    % Ensure fields required by your function exist
    if ~isfield(Cortex,'SulciMap'), Cortex.SulciMap = zeros(size(Cortex.Vertices,1),1); end
    if ~isfield(Cortex,'VertConn')
        if exist('tess_vertconn','file')==2
            Cortex.VertConn = tess_vertconn(Cortex.Faces, size(Cortex.Vertices,1));
        else
            % fallback: no smoothing if we cannot build connectivity
            Cortex.VertConn = [];
        end
    end
    % Match length
    nV = size(Cortex.Vertices,1);
    if numel(J) ~= nV
        warning('Map length (%d) != #vertices (%d). Truncating/padding.', numel(J), nV);
        J = J(:);
        if numel(J) > nV, J = J(1:nV);
        else, J(end+1:nV,1) = 0; end
    end
    % Plot
    PlotSourceTopography(Cortex, J(:), cmap_struct);
    try, sgtitle(ttl); catch, title(ttl); end
end

function cmap = make_default_colormap()
% Build a default colorMap struct compatible with your PlotSourceTopography
    if exist('turbo','file')==2
        base = turbo(256);
    else
        base = parula(256);
    end
    cmap = struct('cmap_a', base);
end

% === Your function (verbatim, with tiny safety tweaks) ===
function PlotSourceTopography(Cortex,J,colorMap)
    fig = figure;
    % Try to load template axes; if missing, create new axes
    try
        template = load('axes.mat');
        if isfield(template,'axes') && isvalid(template.axes)
            currentAxes = template.axes;
            set(currentAxes,"Parent",fig);
        else
            currentAxes = axes('Parent',fig);
        end
    catch
        currentAxes = axes('Parent',fig);
    end

    sources_iv = sqrt(abs(J));
    if max(sources_iv(:))>0
        sources_iv = sources_iv / max(sources_iv(:));
    end
    smoothValue             = 0.66;
    SurfSmoothIterations    = 10;

    % Smooth vertices if we have connectivity and tess_smooth
    if isfield(Cortex,'VertConn') && ~isempty(Cortex.VertConn) && exist('tess_smooth','file')==2
        Vertices = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);
    else
        Vertices = Cortex.Vertices;
    end
    if ~isfield(Cortex,'SulciMap') || isempty(Cortex.SulciMap)
        sulci = 0;
    else
        sulci = Cortex.SulciMap*0.06;
        if numel(sulci) ~= size(Vertices,1), sulci = 0; end
    end

    patch(currentAxes, ...
        'Faces',Cortex.Faces, ...
        'Vertices',Vertices, ...
        'FaceVertexCData',sulci + log(1+sources_iv), ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'AlphaDataMapping', 'none', ...
        'EdgeColor',        'none', ...
        'EdgeAlpha',        1, ...
        'BackfaceLighting', 'lit', ...
        'AmbientStrength',  0.5, ...
        'DiffuseStrength',  0.5, ...
        'SpecularStrength', 0.2, ...
        'SpecularExponent', 1, ...
        'SpecularColorReflectance', 0.5, ...
        'FaceLighting',     'gouraud', ...
        'EdgeLighting',     'gouraud', ...
        'FaceAlpha',1);
    set(currentAxes,'xcolor','w','ycolor','w','zcolor','w');
    view(currentAxes,0,0);
    if isstruct(colorMap) && isfield(colorMap,'cmap')
        colormap(currentAxes,colorMap.cmap);
    else
        colormap(currentAxes,parula(256));
    end
    rotate3d(currentAxes,'on');
end

function z = zscore_safe(x)
    mu = mean(x,'omitnan'); sd = std(x,0,'omitnan'); sd(sd==0)=1;
    z = (x - mu) ./ sd;
end
