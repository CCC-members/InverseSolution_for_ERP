function D = dmd_decompose(J, dt, r, cortex, outpref)
% D = dmd_decompose_and_plot(J, dt, r, cortex, outpref)
% Exact DMD on J (nVox x nTime). Returns modes, eigenvalues, frequencies, amplitudes.
% Plots spatial maps of leading modes and DMD spectrum (frequency vs growth).
%
% Inputs:
%  - J: nVox x nTime
%  - dt: timestep (seconds) = 1/Fs
%  - r: truncation rank (if empty, choose via energy)
%  - cortex: surface struct (optional) for plotting
%  - outpref: prefix for saving plots (optional)
%
if nargin<4, cortex=[]; end
if nargin<5, outpref='dmd'; end

X1 = double(J(:,1:end-1));
X2 = double(J(:,2:end));
[U,Sigma,V] = svd(X1,'econ');

if isempty(r)
    svals = diag(Sigma);
    cum = cumsum(svals.^2) / sum(svals.^2);
    r = find(cum>0.99,1); % keep 99% energy
    r = max(r,10);
end

Ur = U(:,1:r); Sr = Sigma(1:r,1:r); Vr = V(:,1:r);

Atilde = Ur' * X2 * Vr / Sr;
[W,D_eigs] = eig(Atilde);
Phi = X2 * Vr / Sr * W; % DMD modes (nVox x r)
lambda = diag(D_eigs);
omega = log(lambda)/dt;
freqs = imag(omega)/(2*pi); % Hz
growth = real(omega);

% amplitudes (initial condition least-squares)
x1 = X1(:,1);
b = Phi \ x1;

D.Phi = Phi; D.lambda = lambda; D.freqs = freqs; D.growth = growth; D.b = b;
D.r = r;

% plot spectrum (freq vs growth)
figure('Name','DMD Spectrum','Color','w');
scatter(freqs, growth, 40, abs(b),'filled');
xlabel('Frequency (Hz)'); ylabel('Growth rate (1/s)');
title('DMD modes: frequency vs growth (color=amplitude)');
colorbar; grid on;
if ~isempty(outpref), saveas(gcf,[outpref '_spectrum.png']); end

% Plot spatial maps of leading modes (by amplitude)
[~,ord] = sort(abs(b),'descend');
nplot = min(8, numel(ord));
for k=1:nplot
    idx = ord(k);
    mode_map = real(Phi(:,idx));
    figure('Name',['DMD mode ' num2str(k)],'Color','w','Units','normalized','Position',[0.4 0.05 0.5 0.4]);
    if ~isempty(cortex)
        patch('Vertices',cortex.Vertices,'Faces',cortex.Faces,...
              'FaceVertexCData',mode_map,'FaceColor','interp','EdgeColor','none');
        axis equal; axis off; camlight headlight; material dull;
    else
        plot(mode_map); title(['DMD mode ' num2str(idx)]);
    end
    colorbar;
    subtitle = sprintf('freq=%.3f Hz, growth=%.4f, amp=%.3e', freqs(idx), growth(idx), abs(b(idx)));
    title(subtitle);
    if ~isempty(outpref), saveas(gcf,[outpref '_mode' num2str(k) '.png']); end
end

end
