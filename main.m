%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Inverse Solution for ERP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dependencies
% EEGLAB (2024) or higher
%
%
%
% Authors
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
%
%    November 15, 2025

%%
%% Preparing WorkSpace
%%
clc;
close all;
restoredefaultpath;
clearvars -except varargin;

disp('==========================================================================');
disp('==========================================================================');
disp('-->> Starting process');
disp('==========================================================================');
disp('==========================================================================');

%%
%% Preparing dependencies
%%
addpath(genpath('functions'));
addpath(genpath('tools'));
addpath(genpath('external'));
addpath(genpath('templates'));
addpath(genpath('data'));


%%
%%  Defining properties
%%
properties.activation_params = jsondecode(fileread("app/activation_params.json"));
properties.sensor_params = jsondecode(fileread("app/sensor_params.json"));
properties.general_params = jsondecode(fileread("app/general_params.json"));
properties.anat_params = jsondecode(fileread("app/anat_params.json"));
properties.event_params = jsondecode(fileread("app/event_params.json"));


%%
%%  Statics parameters
%%
cmap                    = load("tools/mycolormap_brain_basic_conn.mat");
Fs                      = properties.sensor_params.samp_freq.value;
wname                   = 'db4';
levels                  = 7;
band_names              = {'32-64Hz', '16-32Hz', '8-16Hz', '4-8Hz', '2-4Hz', '1-2Hz', '0.5-1Hz','<0.5Hz'};
base_path               = properties.general_params.workspace.input_path;
output_path             = properties.general_params.workspace.output_path;



if(isequal(properties.anat_params.type,"template"))
    anat_params = properties.anat_params.type_list(1);
    disp('-->> Loading structural data for template anatomy');
    Channels                    = load(fullfile(anat_params.base_path,anat_params.channels));
    Leadfield                   = load(fullfile(anat_params.base_path,anat_params.leadfield));
    Cortex                      = load(fullfile(anat_params.base_path,anat_params.cortex_low));
    CortexHigh                  = load(fullfile(anat_params.base_path,anat_params.cortex_high));    
    Scouts                      = Cortex.Atlas(Cortex.iAtlas).Scouts;
    ScoutsHigh                  = CortexHigh.Atlas(CortexHigh.iAtlas).Scouts;
    ref_labels                  = jsondecode(fileread(fullfile(anat_params.base_path,anat_params.labels)));
end

subjects = dir(fullfile(base_path));
subjects(ismember({subjects.name},{'.','..','derivatives','.datalad','labels','structural'})) = [];
subjects = subjects([subjects.isdir]);
%%
%%  Propcessing EEG Data
%%
participants_file = fullfile(base_path,properties.general_params.workspace.part_file);
isPartInfo = false;
[~,~,part_ext] = fileparts(participants_file);
if(isfile(participants_file))
    if(isequal(part_ext,'.tsv') || isequal(part_ext,'.csv'))
        isPartInfo = true;
        opts = detectImportOptions(participants_file, FileType="text");
        participants = table2struct(readtable(fullfile(participants_file),opts));
    end
    if(isequal(part_ext,'.json'))
        isPartInfo = true;
        participants = jsondecode(fileread(participants_file));
    end
end

addpath(properties.general_params.dependencies.eeglab.base_path);
eeglab nogui;



for i=1:length(subjects)
    subject = subjects(i);
    subID = subject.name;
    disp('--------------------------------------------------------------------------');
    disp(strcat("-->> Processing subject: ",subID));
    disp('--------------------------------------------------------------------------');

    subject_path = fullfile(output_path,subID);
    if(~isfolder(subject_path))
        mkdir(subject_path);
    end

    %%
    %% Loading Structural data
    %%
    if(isequal(properties.anat_params.type,'individual'))
        anat_params = properties.anat_params.type_list(2);
        disp('-->> Loading structural data for individual anatomy');
        Channels                    = load(fullfile(anat_params.base_path,subID,anat_params.channels));
        Leadfield                   = load(fullfile(anat_params.base_path,subID,anat_params.leadfield));
        Cortex                      = load(fullfile(anat_params.base_path,subID,anat_params.cortex_low));
        CortexHigh                  = load(fullfile(anat_params.base_path,subID,anat_params.cortex_high));
        Scouts                      = Cortex.Atlas(Cortex.iAtlas).Scouts;
        ScoutsHigh                  = CortexHigh.Atlas(CortexHigh.iAtlas).Scouts;
        ref_labels                  = jsondecode(fileread(fullfile(anat_params.base_path,anat_params.labels)));
    end

    %%
    %% Loading Functional data
    %%
    EEG                     = ImportEEG(properties, subID);    
    EEG                     = remove_eeg_channels_by_labels(ref_labels, EEG);
    labels                  = {EEG.chanlocs.labels};

    %%
    %% Import events
    %%
    EEG                     = ImportEvents(properties,EEG);       

    %%
    %% Select conditions
    %%
    EEGs                    = SelectConditions(properties,EEG); 

    %%
    %% Structural integration
    %%
    disp('--------------------------------------------------------------------------');
    disp(strcat("-->> Structural preprocessing"));
    disp('--------------------------------------------------------------------------');
    % Gain = bst_gain_orient(leadfield.Gain, leadfield.GridOrient, leadfield.GridAtlas);

    %%
    %% Loadding Structural data
    %%
    

    nV                          = size(Cortex.Vertices,1);
    fixmap                      = @(v) reshape(v(1:min(numel(v),nV)), [], 1);

    [rChannels,rLeadfield]      = remove_channels_by_preproc_data(labels,Channels,Leadfield);
    [sChannels,sLeadfield]      = sort_channels_by_preproc_data(labels,rChannels,rLeadfield);
    subject.modality            = 'EEG';
    subject.Headmodel           = sLeadfield;
    subject.Scortex             = Cortex;
    [subject,properties]        = get_activation_priors(subject,properties);

    struct_path = fullfile(subject_path,'structural');
    if(~isfolder(struct_path))
        mkdir(struct_path);
    end
    save(fullfile(struct_path,'Cortex.mat'),'-struct','Cortex');
    save(fullfile(struct_path,'CortexHigh.mat'),'-struct','CortexHigh');
    save(fullfile(struct_path,'Leadfield.mat'),'-struct','Leadfield');
    save(fullfile(struct_path,'Channels.mat'),'-struct','Channels');
      
    for e=1:length(EEGs)
        segmentEEG                          = EEGs(e);
        seg_path = fullfile(subject_path,'eeg');
        if(~isfolder(seg_path))
            mkdir(seg_path);
        end
        save(fullfile(seg_path,strcat('EEG_cond-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'.mat')),'-struct','segmentEEG');
        subject.MEEG                        = segmentEEG;
        data                                = subject.MEEG.data;
        V                                   = modwt(data.', wname, levels);
        V                                   = permute(V,[3 2 1]);
        J                                   = zeros(length(subject.Scortex.Vertices),size(V,2),size(V,3));
        for lev = 1:(levels+1)
            Svv                             = cov(squeeze(V(:,:,lev)).');
            % PlotSensorTopography(Svv,sChannels,strcat("Topo plot ",band_names{lev}));

            sensor_level_out.Svv            = Svv;
            sensor_level_out.band.str_band  = band_names{lev};
            subject.sensor_level_out        = sensor_level_out;
            [subject,properties,outputs]    = activation_level_sssblpp(subject,properties);
            Tjv                             = outputs.T;
            Tjv                             = bst_gain_orient(Tjv', Leadfield.GridOrient, Leadfield.GridAtlas);
            Tjv                             = Tjv';
            J(:,:,lev)                      = Tjv*squeeze(V(:,:,lev));
            % PlotSourceTopography_safe(Cortex, squeeze(var(J(:,:,lev),0,2)), cmap, strcat(subID,'-','Eyes closed',band_names{lev}));
        end

        J                                   = permute(J,[3 2 1]);
        J                                   = imodwt(J, wname).';
        segmentEEG.J                        = J;
        out                                 = TimeSeriesAnalysis(segmentEEG,Fs,Cortex);
        out.J                               = J;

        fig_path = fullfile(subject_path,'Figures');
        if(~isfolder(fig_path))
            mkdir(fig_path);
        end
        disp(strcat("---->> Saving data"));
        fig_delta                           =  PlotSourceTopography(Cortex, fixmap(out.BandMaps.delta), cmap, strcat(subID,'-',segmentEEG.condition,'-seg-',num2str(segmentEEG.segment),'-Delta band power (0.5–4 Hz)'));
        saveas(fig_delta,fullfile(fig_path,strcat(subID,'-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'-Delta band power (0.5–4 Hz)','.fig')));
        close(fig_delta);
        fig_alpha                            = PlotSourceTopography(Cortex, fixmap(out.BandMaps.alpha), cmap, strcat(subID,'-',segmentEEG.condition,'-seg-',num2str(segmentEEG.segment),'-Alpha band power (8–14 Hz)'));
        saveas(fig_alpha,fullfile(fig_path,strcat(subID,'-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'-Alpha band power (8–14 Hz)','.fig')));
        close(fig_alpha);


        %%
        %% Scout map by subject
        %%
        ScoutHighCDataDelta              = zeros(length(CortexHigh.Vertices),1);
        ScoutHighCDataAlpha              = zeros(length(CortexHigh.Vertices),1);
        for s=1:length(Scouts)
            scout                           = Scouts(s);
            scoutHigh                       = ScoutsHigh(s);
            ScoutHighCDataDelta(scoutHigh.Vertices) = median(out.BandMaps.delta(scout.Vertices),1);
            ScoutHighCDataAlpha(scoutHigh.Vertices) = median(out.BandMaps.alpha(scout.Vertices),1);
        end
        fig_delta = PlotScoutTopography(CortexHigh,ScoutHighCDataDelta,cmap,strcat(subID," - ",segmentEEG.condition,'-seg-',num2str(segmentEEG.segment)," - Delta - Scouts map (0.1–4 Hz)"));
        saveas(fig_delta,fullfile(fig_path,strcat(subID,'-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'- Delta - Scouts map (0.1–4 Hz)','.fig')));
        close(fig_delta);
        fig_alpha = PlotScoutTopography(CortexHigh,ScoutHighCDataAlpha,cmap,strcat(subID," - ",segmentEEG.condition,'-seg-',num2str(segmentEEG.segment)," - Alpha - Scouts map (8–14 Hz)"));
        saveas(fig_alpha,fullfile(fig_path,strcat(subID,'-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'- Alpha - Scouts map (8–14 Hz)','.fig')));
        close(fig_alpha);

        %%
        %%  Time series by Scouts (Eyes open)
        %%
        J = segmentEEG.J;
        JScout = zeros(length(Scouts),size(J,2));
        for s=1:length(Scouts)
            scout = Scouts(s);
            sVertices = scout.Vertices;
            JScout(s,:) = median(J(scout.Vertices,:),1);
        end
        out.JScouts = JScout;

        out_path = fullfile(subject_path,'time');
        if(~isfolder(out_path))
            mkdir(out_path);
        end
        save(fullfile(out_path,strcat('Times_cond-',segmentEEG.condition,'_seg-',num2str(segmentEEG.segment),'.mat')),'-struct','out','-v7.3');

        disp('--------------------------------------------------------------------------');
        disp('--------------------------------------------------------------------------');
    end

end
cmap = load("tools/mycolormap_brain_basic_conn.mat");
nV = size(Cortex.Vertices,1);
fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);

%%
%%  EO average (Delta & Alpha)
%%
BandMapsDelta = [BandMaps.delta];
BandMapsDelta = median(BandMapsDelta,2);
BandMapsAlpha = [BandMaps.alpha];
BandMapsAlpha = median(BandMapsAlpha,2);
fig_delta_ave = PlotSourceTopography_safe(Cortex, fixmap(BandMapsDelta), cmap, strcat('Eyes closed','-Delta band average power (0.1–4 Hz)'));
saveas(fig_delta_ave,fullfile(output_path,strcat('Eyes closed','-Delta band average power (0.1–4 Hz)','.fig')));
close(fig_delta_ave);
fig_alpha_ave = PlotSourceTopography_safe(Cortex, fixmap(BandMapsAlpha), cmap, strcat('Eyes closed','-Alpha band average power (8–14 Hz)'));
saveas(fig_alpha_ave,fullfile(output_path,strcat('Eyes closed','-Alpha band average power (8–14 Hz)','.fig')));
close(fig_alpha_ave);

% EO Scouts
ScoutHighCDataDelta = zeros(length(CortexHigh.Vertices),1);
ScoutHighCDataAlpha = zeros(length(CortexHigh.Vertices),1);
for s=1:length(Scouts)
    scout = Scouts(s);
    scoutHigh = ScoutsHigh(s);
    ScoutHighCDataDelta(scoutHigh.Vertices) = median(BandMapsDelta(scout.Vertices),1);
    ScoutHighCDataAlpha(scoutHigh.Vertices) = median(BandMapsAlpha(scout.Vertices),1);
end
nV = size(CortexHigh.Vertices,1);
fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);
fig_scout_delta_ave = PlotScoutTopography(CortexHigh,ScoutHighCDataDelta,cmap,strcat(subID," - ",'Eyes closed',"- Average - Delta - Scouts map (0.1–4 Hz)"));
saveas(fig_scout_delta_ave,fullfile(output_path,strcat('Eyes closed','-Average-Delta-Scouts map (0.1–4 Hz)','.fig')));
close(fig_scout_delta_ave);
fig_scout_alpha_ave = PlotScoutTopography(CortexHigh,ScoutHighCDataAlpha,cmap,strcat(subID," - ",'Eyes closed',"- Average - Alpha - Scouts map (8–14 Hz)"));
saveas(fig_scout_alpha_ave,fullfile(output_path,strcat('Eyes closed','-Average-Alpha-Scouts map (8–14 Hz)','.fig')));
close(fig_scout_alpha_ave);


% Pxx_ave = median(Pxx,3);
% bands = struct('delta',[0.1 4],'theta',[4 7], ...
%     'alpha',[8 14],'alpha10',[9.5 10.5],'beta',[14 30],'gammaL',[30 50]);
%
% bandNames = fieldnames(bands);
% BandMaps = struct();
% for i=1:numel(bandNames)
%     rngHz = bands.(bandNames{i});
%     BandMaps.(bandNames{i}) = bandpower_map(Pxx_ave, out.f, rngHz);  % [S x 1]
% end
% PlotSourceTopography_safe(Cortex, fixmap(BandMaps.delta), cmap, strcat(subID,'-','Eyes closed','-Delta band power (0.1–4 Hz)'));
% PlotSourceTopography_safe(Cortex, fixmap(BandMaps.alpha), cmap, strcat(subID,'-','Eyes closed','-Alpha band power (8–14 Hz)'));
%
% Pxx_ave = median(Pxx,3);
% for i=1:numel(bandNames)
%     rngHz = bands.(bandNames{i});
%     BandMaps.(bandNames{i}) = bandpower_map(Pxx_ave, out.f, rngHz);  % [S x 1]
% end
% cmap = load("tools/mycolormap_brain_basic_conn.mat");
% nV = size(Cortex.Vertices,1);
% fixmap = @(v) reshape(v(1:min(numel(v),nV)), [], 1);
% PlotSourceTopography_safe(Cortex, fixmap(BandMaps.delta), cmap, strcat(EEG_ec.SubID,'-',EEG_ec.Condition,'-Delta band power (0.1–4 Hz)'));
% PlotSourceTopography_safe(Cortex, fixmap(BandMaps.alpha), cmap, strcat(EEG_ec.SubID,'-',EEG_ec.Condition,'-Alpha band power (8–14 Hz)'));
%
%
%

%%
%%  Saving the Time Series
%%
J_Ave = median(TemiSeriesTensor,3);
J_Scout_Ave = median(Scout_tensor,3);
EEG_Ave.J = J_Ave;
EEG_Ave.JScout = J_Scout_Ave;
EEG_Ave.SubID = 'Ave';
EEG_Ave.Condition = 'Eyes open';
% TimeSeriesAnalysis(EEG_Ave,Fs,Cortex);

J_Ave_ec = median(tensor,3);
J_Scout_Ave_ec = median(scout_tensor,3);
EEG_Ave.J = J_Ave_ec;
EEG_Ave.JScout = J_Scout_Ave_ec;
EEG_Ave.SubID = 'Ave';
EEG_Ave.Condition = 'Eyes closed';
% TimeSeriesAnalysis(EEG_Ave,Fs,Cortex);

disp(strcat("---->> Saving data"));
save(fullfile(output_path,'Time_Ave.mat'),'J_Ave');
save(fullfile(output_path,'Time_Ave_Scouts.mat'),'J_Scout_Ave');
save(fullfile(output_path,'TimeSeries.mat'),'TemiSeriesTensor');

save(fullfile(output_path,'Time_Ave.mat'),'J_Ave_ec');
save(fullfile(output_path,'Time_Ave_Scouts.mat'),'J_Scout_Ave_ec');
save(fullfile(output_path,'TimeSeries_ec.mat'),'tensor');


%%
%%  Export figures
%%
%export_fig(fig,fullfile(figure.folder,figure.name),'-transparent','-png');


disp('==========================================================================');
disp('==========================================================================');
disp('-->> Process finished');
disp('==========================================================================');
disp('==========================================================================');




