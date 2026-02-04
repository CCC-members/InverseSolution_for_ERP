function S = compute_spectra(J, Fs, cortex, outpref, freqrange)
% S = compute_spectra_and_plot(J, Fs, cortex, outpref, freqrange)
% Computes power spectral density (Welch) for each voxel (row of J).
% Plots:
%  1) mean PSD across voxels
%  2) spatial map of band power (user-specified bands) on cortex
%
% Inputs:
%  - J: [nVox x nTime]
%  - Fs: sampling frequency (Hz)
%  - cortex: struct with vertices and faces (Brainstorm format)
%  - outpref: string prefix for saving figs (optional)
%  - freqrange: [min max] to display (optional)
%
% Output:
%  - S: struct with fields .Pxx (nVox x nFreq), .F (1 x nFreq)
%
if nargin<4, outpref='spec'; end
if nargin<5, freqrange=[]; end

[nVox, nT] = size(J);
win = round(2*Fs); % 2-second window
noverlap = round(0.5*win);
nfft = max(2^nextpow2(win), 512);

% Preallocate
Pxx = [];
for v=1:nVox
    x = double(J(v,:));
    [pxx, f] = pwelch(x, win, noverlap, nfft, Fs);
    if v==1
        Pxx = zeros(nVox, numel(pxx));
    end
    Pxx(v,:) = pxx(:)';
end

S.Pxx = Pxx;
S.F = f;

% Plot mean PSD
figure('Name','Mean PSD across voxels','Color','w','Units','normalized','Position',[0.05 0.6 0.4 0.3]);
plot(f,10*log10(mean(Pxx,1)),'LineWidth',1.2);
xlabel('Frequency (Hz)'); ylabel('Power (dB/Hz)');
title('Mean PSD across voxels');
grid on;
if ~isempty(freqrange)
    xlim(freqrange);
end
if ~isempty(outpref), saveas(gcf,[outpref '_meanPSD.png']); end

% Define standard bands and compute band-power maps
bands = struct('name',{'Delta','Theta','Alpha','Beta','Gamma'}, ...
               'range',{[1 4],[4 8],[8 13],[13 30],[30 60]});
bandmap = zeros(nVox, numel(bands));
for b=1:numel(bands)
    idx = f>=bands(b).range(1) & f<=bands(b).range(2);
    bandmap(:,b) = trapz(f(idx), Pxx(:,idx),2);
end

% Plot each band on cortex
for b=1:numel(bands)
    figure('Name',['Band ' bands(b).name],'Color','w','Units','normalized','Position',[0.5 0.6 0.4 0.3]);
    patch('Vertices',cortex.Vertices,'Faces',cortex.Faces,...
          'FaceVertexCData',bandmap(:,b),'FaceColor','interp','EdgeColor','none');
    axis equal; axis off; camlight headlight; material dull;
    title(['Band power: ' bands(b).name ' (' num2str(bands(b).range) ' Hz)']);
    colorbar;
    if ~isempty(outpref), saveas(gcf,[outpref '_band_' bands(b).name '.png']); end
end

end

