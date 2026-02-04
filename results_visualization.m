%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Information Theoretical Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
%


%%
%% Preparing WorkSpace
%%
clc;
close all;
restoredefaultpath;
clearvars -except varargin;
addpath(genpath('functions'));
addpath(genpath('tools'));

disp('==========================================================================');
disp('==========================================================================');
disp(" - ->> Starting process");
disp('==========================================================================');
disp('==========================================================================');

output_path = "Outputs";
Cortex = load(fullfile(output_path,'structural','Cortex.mat'));
CortexHigh = load(fullfile(output_path,'structural','Cortex_high.mat'));

%%
%%  Parcels correlation
%%
Scouts = Cortex.Atlas(Cortex.iAtlas).Scouts;
ScoutsHigh = CortexHigh.Atlas(CortexHigh.iAtlas).Scouts;


cmap = load("tools/mycolormap_brain_basic_conn.mat");
nV = size(Cortex.Vertices,1);
fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);

bands = struct('delta',[0.1 4],'theta',[4 7], ...
    'alpha',[8 14],'alpha10',[9.5 10.5],'beta',[14 30],'gammaL',[30 50]);
bandNames = fieldnames(bands);

%%
%%  Subjects Visualization
%%
subjects = dir(output_path);
subjects = subjects([subjects.isdir]);
subjects(ismember({subjects.name},{'.','..','structural'})) = [];

for i=1:length(subjects)
    subject = subjects(i);

    disp(" - -------------------------------------------------------------------------");
    disp(strcat("-->> Processing subject: ",subject.name));
    disp(" - -------------------------------------------------------------------------");
    EEG_eo = load(fullfile(subject.folder,subject.name,"eyes_open","EEG.mat"));
    EEG_ec = load(fullfile(subject.folder,subject.name,"eyes_closed","EEG.mat"));
    
    %%
    %% Frequencies average by Attlas Scouts (Eyes open)
    %%
    BandMaps_eo = EEG_eo.BandMaps;
    % Source map by subject
    PlotSourceTopography_safe(Cortex, fixmap(BandMaps_eo.delta), cmap, strcat(EEG_eo.SubID," - ",EEG_eo.Condition," - Delta band power (0.1–4 Hz)"));
    PlotSourceTopography_safe(Cortex, fixmap(BandMaps_eo.alpha), cmap, strcat(EEG_eo.SubID," - ",EEG_eo.Condition," - Alpha band power (8–14 Hz)"));

    % Scout map by subject 
    ScoutHighCDataDelta_eo = zeros(length(CortexHigh.Vertices),1);
    ScoutHighCDataAlpha_eo = zeros(length(CortexHigh.Vertices),1);
    for s=1:length(Scouts)
        scout = Scouts(s);
        scoutHigh = ScoutsHigh(s);
        ScoutHighCDataDelta_eo(scoutHigh.Vertices) = median(BandMaps_eo.delta(scout.Vertices),1);
        ScoutHighCDataAlpha_eo(scoutHigh.Vertices) = median(BandMaps_eo.alpha(scout.Vertices),1);
    end
    nV = size(CortexHigh.Vertices,1);
    fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);
    PlotScoutTopography(CortexHigh,ScoutHighCDataDelta_eo,cmap,strcat(EEG_eo.SubID," - ",EEG_eo.Condition," - Delta - Scouts map (0.1–4 Hz)"));
    PlotScoutTopography(CortexHigh,ScoutHighCDataAlpha_eo,cmap,strcat(EEG_eo.SubID," - ",EEG_eo.Condition," - Alpha - Scouts map (8–14 Hz)"));


    %%
    %% Frequencies average by Attlas Scouts (Eyes closed)
    %%
    BandMaps_ec = EEG_ec.BandMaps;
    % Source map by subject
    PlotSourceTopography_safe(Cortex, fixmap(BandMaps_ec.delta), cmap, strcat(EEG_ec.SubID," - ",EEG_ec.Condition," - Delta band power (0.1–4 Hz)"));
    PlotSourceTopography_safe(Cortex, fixmap(BandMaps_ec.alpha), cmap, strcat(EEG_ec.SubID," - ",EEG_ec.Condition," - Alpha band power (8–14 Hz)"));

    % Scout map by subject 
    ScoutHighCDataDelta_ec = zeros(length(CortexHigh.Vertices),1);
    ScoutHighCDataAlpha_ec = zeros(length(CortexHigh.Vertices),1);
    for s=1:length(Scouts)
        scout = Scouts(s);
        scoutHigh = ScoutsHigh(s);
        ScoutHighCDataDelta_ec(scoutHigh.Vertices) = median(BandMaps_ec.delta(scout.Vertices),1);
        ScoutHighCDataAlpha_ec(scoutHigh.Vertices) = median(BandMaps_ec.alpha(scout.Vertices),1);
    end
    nV = size(CortexHigh.Vertices,1);
    fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);
    PlotScoutTopography(CortexHigh,ScoutHighCDataDelta_ec,cmap,strcat(EEG_ec.SubID," - ",EEG_ec.Condition," - Delta - Scouts map (0.1–4 Hz)"));
    PlotScoutTopography(CortexHigh,ScoutHighCDataAlpha_ec,cmap,strcat(EEG_ec.SubID," - ",EEG_ec.Condition," - Alpha - Scouts map (8–14 Hz)"));    
end

disp('==========================================================================');
disp('==========================================================================');
disp(" - ->> Process finished");
disp('==========================================================================');
disp('==========================================================================');