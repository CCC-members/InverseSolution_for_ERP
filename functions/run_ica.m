function ICA = run_ica(J, nComp, cortex, outpref)
% ICA = run_ica_and_plot(J, nComp, cortex, outpref)
% Runs ICA on J (nVox x nTime). Plots spatial maps for each component and component timecourses.
% Returns structure with .mixing (nVox x nComp), .sources (nComp x nTime)
%
if nargin<3, cortex = []; end
if nargin<4, outpref='ica'; end

[nVox, nT] = size(J);
% center and whiten via PCA (common preproc)
X = double(J);
X = X - mean(X,2);
% PCA
[U,Sigma,~] = svd((X*X')/nT);
% choose nComp
if isempty(nComp), nComp = min(20, size(U,2)); end
Ured = U(:,1:nComp);
% whiten
whitenMat = diag(1./sqrt(diag(Sigma(1:nComp,1:nComp)))) * Ured';
Xw = whitenMat * X;

% Try to use fastica if available
use_fastica = exist('fastica','file')==2;
if use_fastica
    try
        [icasig, A, W] = fastica(X, 'numOfIC', nComp, 'verbose','off');
        sources = icasig; % nComp x nT
        mixing = A;       % nVox x nComp
    catch
        use_fastica = false;
    end
end

% fallback: simple symmetric FastICA implementation
if ~use_fastica
    % symmetric decorrelation algorithm
    rng(0);
    W = randn(nComp);
    maxIter = 1000; tol = 1e-6;
    for it=1:maxIter
        Wold = W;
        % g(u)=tanh(u)
        gwtx = tanh(W * Xw);
        g_wtx = 1 - gwtx.^2;
        W = (gwtx * Xw')/nT - diag(mean(g_wtx,2))*W;
        % symmetric decorrelation
        [Uw,Sw] = svd(W);
        W = Uw * eye(nComp) * Uw' * W;
        if max(abs(abs(diag(W*Wold'))-1)) < tol, break; end
    end
    sources = W * Xw; % nComp x nT
    % mixing matrix approximate inverse
    mixing = pinv(W*whitenMat);
end

ICA.sources = sources;
ICA.mixing = mixing;

% Plot topographies (mixing columns) on cortex
for c=1:size(mixing,2)
    map = mixing(:,c);
    figure('Name',['ICA comp ' num2str(c)],'Color','w','Units','normalized','Position',[0.05 0.05 0.25 0.4]);
    if ~isempty(cortex)
        patch('Vertices',cortex.Vertices,'Faces',cortex.Faces,...
              'FaceVertexCData',map,'FaceColor','interp','EdgeColor','none');
        axis equal; axis off; camlight headlight; material dull;
    else
        subplot(2,1,1);
        plot(map); title(['ICA mixing component ' num2str(c)]);
    end
    colorbar;
    % timecourse
    subplot(2,1,2);
    plot(sources(c,:));
    xlabel('Time (samples)'); ylabel('Amplitude');
    title(['ICA source ' num2str(c)]);
    if ~isempty(outpref), saveas(gcf,[outpref '_comp' num2str(c) '.png']); end
end

end
