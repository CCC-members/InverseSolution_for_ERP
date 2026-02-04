function SSA = ssa_decompose(x, L, nRec, outpref)
% SSA = ssa_decompose_and_plot(x, L, nRec, outpref)
% Singular Spectrum Analysis for a single time series x (1 x T) or matrix (nComp x T) (works per-row)
% - L: embedding window length
% - nRec: number of leading components to keep for reconstruction (scalar or array)
% Plots eigen-spectrum, selected reconstructed components, and residual.
%
if nargin<4, outpref='ssa'; end
if size(x,1) > 1
    % apply SSA to each row and return cell
    for r=1:size(x,1)
        SSA{r} = ssa_decompose_and_plot(x(r,:), L, nRec, [outpref '_row' num2str(r)]);
    end
    return;
end

x = x(:)';
T = length(x);
if L >= T, error('L must be < length(x)'); end
K = T - L + 1;

% Build trajectory matrix
X = zeros(L, K);
for i=1:K
    X(:,i) = x(i:i+L-1)';
end

[U,S,V] = svd(X, 'econ');
d = diag(S);
% eigen-spectrum
figure('Name','SSA eigenspectrum','Color','w');
plot(d,'o-','LineWidth',1.2); xlabel('Component'); ylabel('Singular value'); title('SSA singular values');
if ~isempty(outpref), saveas(gcf,[outpref '_eigs.png']); end

% reconstruct first nRec components (grouping)
if isempty(nRec), nRec = 2; end
if numel(nRec)==1
    groups = {1:nRec};
else
    groups = nRec;
end

recon = zeros(length(groups), T);
for g=1:numel(groups)
    inds = groups{g};
    Xg = zeros(L,K);
    for kidx = inds
        Xg = Xg + d(kidx) * U(:,kidx) * V(:,kidx)';
    end
    % diagonal averaging to reconstruct series
    r = zeros(1,T);
    counts = zeros(1,T);
    for i=1:K
        for j=1:L
            r(i+j-1) = r(i+j-1) + Xg(j,i);
            counts(i+j-1) = counts(i+j-1) + 1;
        end
    end
    r = r ./ counts;
    recon(g,:) = r;
end

% plot reconstructions and residual
figure('Name','SSA Reconstructions','Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.4]);
nplots = size(recon,1) + 1;
subplot(nplots,1,1);
plot(x); title('Original'); ylabel('Amp');
for g=1:size(recon,1)
    subplot(nplots,1,g+1);
    plot(recon(g,:)); title(['SSA recon group ' num2str(g)]);
end
if ~isempty(outpref), saveas(gcf,[outpref '_recon.png']); end

SSA.U = U; SSA.S = S; SSA.V = V; SSA.recon = recon;

end
