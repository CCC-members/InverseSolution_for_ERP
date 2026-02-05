function out = TimeSeriesAnalysis(EEG, J, fs, Cortex, opts)
% TIME_SERIES_ANALYSIS
% PSD (0.1–50 Hz), ICA, MSSA, DMD on inverse solution.
% X: [T x S] or [T x S x K]  (single-trial -> K==1)
% fs: Hz, Cortex: Brainstorm tess struct
% opts: nICA, ssa_L, ssa_r, dmd_r, r_pod, show_plots, bands, pod_energy, ridge, rcond
%       spectra_fmax, spectra_overlayN
    
    X = J';
    clear J;

    if nargin < 5, opts = struct; end
    if ndims(X) == 2, X = reshape(X, size(X,1), size(X,2), 1); end
    [T,S,K] = size(X);

    % ---------- defaults ----------
    def.ssa_L = max(10, round(2*fs));
    def.ssa_r = 20;
    def.dmd_r = min(40, T-2);
    def.r_pod = min(80, T-2);
    def.nICA  = min(30, def.r_pod);
    def.show_plots = true;
    def.bands = struct('delta',[0.5 4],'theta',[4 7], ...
                       'alpha',[8 14],'beta',[14 30],'gammaL',[30 128]);
    % numerical robustness knobs
    def.pod_energy = 0.999; % keep 99.9% temporal energy
    def.ridge      = 1e-6;  % ridge λ for LS
    def.rcond      = 1e-10; % rcond for pinv/rank
    % spectra display (line plots)
    def.spectra_fmax     = 38;  % show up to 48 Hz
    def.spectra_overlayN = 300; % number of individual source lines (0 = none)

    fns = fieldnames(def);
    for i=1:numel(fns)
        if ~isfield(opts, fns{i}), opts.(fns{i}) = def.(fns{i}); end
    end

    % ---------- detrend/center ----------
    X = detrend(X, 0);
    X = bsxfun(@minus, X, mean(X,1));

    % ---------- 1) PSD via Welch ----------
    disp(strcat("---->> PSD via Welch"));
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
        % ---- Temporal POD (economy SVD) with adaptive rank ----
        disp(strcat("---->> Temporal POD (economy SVD) with adaptive rank"));
        [U_t, Sigma, V_s] = svd(X(:,:,1), 'econ');  % X ≈ U_t * Sigma * V_s'
        sing = diag(Sigma);

        tol = max(size(X(:,:,1))) * eps(max(sing));
        r_num = sum(sing > tol);
        cs = cumsum(sing.^2);
        r_energy = find(cs >= opts.pod_energy * cs(end), 1, 'first');
        if isempty(r_energy), r_energy = r_num; end
        r = min([opts.r_pod, r_num, r_energy]);

        U_t = U_t(:,1:r);
        V_s = V_s(:,1:r);
        Sigma_r = Sigma(1:r,1:r);

        T_scores = U_t * Sigma_r;   % [T x r]
        S_loads  = V_s;             % [S x r]

        % ---- ICA on temporal PCs, map back to space ----
        disp(strcat("---->> ICA on temporal PCs, map back to space"));
        [IC_time, IC_space] = ica_on_temporal_scores(T_scores, S_loads, opts.nICA, opts);

        % ---- MSSA on temporal PCs, map back to space ----
        disp(strcat("---->> MSSA on temporal PCs, map back to space"));
        L = min(opts.ssa_L, T-1);
        [SSA_time, SSA_sing] = ssa_multichannel(T_scores, L, min(opts.ssa_r, r), opts);
        W_ssa     = safe_regress(T_scores, SSA_time, opts.ridge, opts.rcond); % [r x rssa]
        SSA_space = S_loads * W_ssa;

        % ---- Projected DMD (on temporal PCs), map back ----
        disp(strcat("---->> Projected DMD (on temporal PCs), map back"));
        DMD = dmd_projected(T_scores, fs, min(opts.dmd_r, r), opts);
        keep = DMD.freq_hz >= 0.1 & DMD.freq_hz <= 50 & isfinite(DMD.freq_hz);
        DMD_time  = DMD.modes_time(:, keep);
        W_dmd     = safe_regress(T_scores, DMD_time, opts.ridge, opts.rcond);
        DMD_modes = S_loads * W_dmd;
        DMD_freqs = DMD.freq_hz(keep);
        DMD_amps  = DMD.amplitudes(keep);
    else
        % ---------- 2) ICA (multi-trial fallback) ----------
        disp(strcat("---->> ICA (multi-trial fallback)"));
        Xica = reshape(permute(X, [1 3 2]), T*K, S);
        Xica = bsxfun(@minus, Xica, mean(Xica,1));
        [~, A_ica, S_ica] = run_ica(Xica, min(opts.nICA, T-2));
        IC_time = reshape(S_ica, T, K, []);
        IC_space = A_ica;

        % ---------- 3) MSSA ----------
        disp(strcat("---->> MSSA"));
        L = min(opts.ssa_L, T-1);
        SSA = mssa_original(X, L, opts.ssa_r);
        SSA_time = reshape(SSA.U, T, []);
        SSA_space = SSA.V;
        SSA_sing = SSA.s;

        % ---------- 4) DMD ----------
        disp(strcat("---->> DMD"));
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

    % ---------- plotting (maps + spectra lines 0–48 Hz) ----------
    % if opts.show_plots
    %     if nargin >= 3 && ~isempty(Cortex)
    %         cmap = load("tools/mycolormap_brain_basic_conn.mat");
    %         nV = size(Cortex.Vertices,1);
    %         fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);
    % 
    %         % Alpha & Delta bandpower maps
    %         PlotSourceTopography_safe(Cortex, fixmap(BandMaps.delta), cmap, strcat(EEG.SubID,'-',EEG.Condition,'-Delta band power (0.1–4 Hz)'));
    % 
    % 
    %         PlotSourceTopography_safe(Cortex, fixmap(BandMaps.alpha), cmap, strcat(EEG.SubID,'-',EEG.Condition,'-Alpha band power (8–14 Hz)'));
    % 
    %         % ICA #1
    %         if ~isempty(IC_space),  PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(IC_space(:,1))), cmap, strcat(EEG.SubID,'-',EEG.Condition,'-ICA #1 spatial')); end
    %         % SSA #1
    %         if ~isempty(SSA_space), PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(SSA_space(:,1))), cmap, strcat(EEG.SubID,'-',EEG.Condition,'-SSA #1 spatial')); end
    %         % DMD top amplitude
    %         if ~isempty(DMD_modes)
    %             [~, iTop] = max(DMD_amps);
    %             PlotSourceTopography_safe(Cortex, fixmap(zscore_safe(real(DMD_modes(:,iTop)))), cmap, ...
    %                 sprintf('DMD mode @ %.2f Hz (max amp)', DMD_freqs(iTop)));
    %         end
    %     end
    % 
    %     % === Spectra: overlay line plots up to 48 Hz (mean + multiple sources) ===
    %     plot_spectra_lines(out.Pxx, out.f, opts.spectra_fmax, opts.spectra_overlayN, true);
    % end
end

% ===================== helpers =====================

function [Pxx, f] = welch_per_source_singletrial(X, fs, wlen, nover, nfft)
% X: [T x S], returns Pxx as [S x F]
    [P, f] = pwelch(X, hamming(wlen), nover, nfft, fs, 'psd'); % P: [F x S]
    Pxx = P.'; % [S x F]
end



function B = safe_regress(X, Y, lambda, rcondv)
% Ridge-regularized least squares: (X'X + λI) \ (X'Y), with pinv fallback
    if nargin < 3 || isempty(lambda), lambda = 0; end
    if nargin < 4 || isempty(rcondv), rcondv = 1e-10; end
    [~, p] = size(X);
    if lambda > 0
        B = (X.'*X + lambda*eye(p)) \ (X.'*Y);
    else
        B = pinv(X, rcondv) * Y;
    end
end

function D = dmd_concat_trials(~, ~, ~) %#ok<STOUT,INUSD>
    error('dmd_concat_trials not used in single-trial fast path.');
end

% ===================== PLOTTING (your renderer + spectra lines) =====================



function plot_spectra_lines(Pxx, f, fmax, overlayN, showMean)
% Plot mean spectrum + N individual source spectra as overlaid line plots.
% Pxx: [S x F], f: [F x 1]
% fmax: max frequency to display (Hz)
% overlayN: number of sources to overlay (e.g., 100–500). If 0, just mean.
% showMean: bool, plot mean (thicker) on top.
    if nargin < 3 || isempty(fmax),     fmax = 48; end
    if nargin < 4 || isempty(overlayN), overlayN = 300; end
    if nargin < 5 || isempty(showMean), showMean = true; end

    idx = (f > 0 & f <= fmax);
    f2  = f(idx);
    S   = size(Pxx,1);

    % choose evenly spaced sources to avoid spatial bias
    if overlayN > 0
        src_idx = unique(round(linspace(1, S, min(overlayN,S))));
    else
        src_idx = [];
    end

    % dB conversion
    Pdb = 10*log10(max(Pxx(:,idx), eps));

    figure('Name', sprintf('Source spectra (lines, 0–%.0f Hz)', fmax));
    hold on;

    % individual sources as thin lines
    if ~isempty(src_idx)
        lw = 0.5;
        for k = 1:numel(src_idx)
            plot(f2, Pdb(src_idx(k),:), '-', 'LineWidth', lw);
        end
    end

    % mean spectrum on top
    if showMean
        mu = mean(Pdb, 1);
        plot(f2, mu, '-', 'LineWidth', 2);
        if isempty(src_idx)
            legend({'Mean'}, 'Location','northeast');
        else
            legend({'Sources','Mean'}, 'Location','northeast');
        end
    end

    hold off; grid on;
    xlim([0 fmax]); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    title(sprintf('Source spectra (N=%d overlays) — 0–%.0f Hz', numel(src_idx), fmax));
end

function z = zscore_safe(x)
    mu = mean(x,'omitnan'); sd = std(x,0,'omitnan'); sd(sd==0)=1;
    z = (x - mu) ./ sd;
end
