function EEG  = ImportEEG(properties, subID)
%% Example of batch code to reject bad channels...
%
%


ref_file                = properties.general_params.workspace.ref_file;
dataset_path            = properties.general_params.workspace.input_path;
file_name               = fullfile(dataset_path,subID,strrep(ref_file,'SubID',subID));
[base_path,filename,ext]    = fileparts(file_name);
freq_sample             = properties.sensor_params.samp_freq.value;



%% Step1: Import data.
switch lower(ext)
    case '.edf'
        EEG                     = pop_biosig(file_name);        
        % Convert 10-10 to 10-20 equivalents
        new_labels              = replace(replace(replace(replace(replace({EEG.chanlocs.labels}','-FPz',''),'-REF',''),'-LE',''),'-A12',''),' ','');
        % new_labels              = replace(replace(replace(replace(replace(replace(replace(replace(replace( ...
        %     {EEG.chanlocs.labels}','-FPz',''),'-REF',''),'-LE',''),'-A12',''),' ',''), ...
        %     'T7','T3'),'T8','T4'),'P7','T5'),'P8','T6');
        [EEG.chanlocs.labels]   = new_labels{:};
    case '.set'
        EEG                     = pop_loadset(file_name);
        indx = ismember({EEG.chanlocs.labels},{'T7','T8','P7','P8'});
        labels = {EEG.chanlocs.labels};
        labels(indx) = {'T3','T4','T5','T6'};
        [EEG.chanlocs.labels]   = labels{:};
     
    case '.matrix'
        load(file_name);
        EEG                     = eeg_emptyset;
        EEG.srate               = 128;
        EEG.data                = data;
        EEG.nbchan              = size(EEG.data,1);
        EEG.pnts                = size(EEG.data,2);
        EEG.xmin                = 0;
        EEG.xmax                = EEG.xmin+(EEG.pnts-1)*(1/EEG.srate);
        EEG.times               = (0:EEG.pnts-1)/EEG.srate.*1000;
        if(exist('labels','var'))
            EEG.chanlocs(length(labels)+1:end,:)    = [];
            new_labels                              = labels;
            [EEG.chanlocs.labels]                   = new_labels{:};
        end
    case '.dat'
        EEG                     = pop_loadBCI2000(file_name);   
        
    
    case '.txt'
        EEG                     = eeg_emptyset;
        [~,filename,~]          = fileparts(file_name);
        EEG.filename            = filename;
        EEG.filepath            = filepath;
        EEG.subject             = subID;
        data                    = readmatrix(file_name);
        data                    = data';
        EEG.data                = data;
        EEG.nbchan              = length(EEG.chanlocs);
        EEG.pnts                = size(data,2);
        EEG.srate               = 200;
        EEG.min                 = 0;
        EEG.max                 = EEG.xmin+(EEG.pnts-1)*(1/EEG.srate);
        EEG.times               = (0:EEG.pnts-1)/EEG.srate.*1000;
    case '.mff'
        EEG                     = pop_readegimff( file_name );
end
EEG.SubID                       = subID;
EEG.subject                     = subID;
EEG.setname                     = subID;
EEG.filename                    = filename;
EEG.filepath                    = base_path;



%% Step 6: Downsample the data.
if(EEG.srate > freq_sample)
    EEG                         = pop_resample(EEG, freq_sample);
end

end
