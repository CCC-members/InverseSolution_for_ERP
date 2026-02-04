function [OutEEG, errMessage] = ImportEvents(properties,EEG)


disp(strcat("-->> Importing EEG Events"));
errMessage = [];

events_filename    = strrep(EEG.filename,'eeg','events');
events_file = fullfile(EEG.filepath,strcat(events_filename,'.tsv'));
if(isfile(events_file))
    EEG.event = struct;
    opts = detectImportOptions(events_file, 'FileType', 'text', 'Delimiter', '\t');
    opts.VariableNamesLine = 1;
    opts.ExtraColumnsRule = 'ignore';
    opts.MissingRule = 'fill';
    opts = setvartype(opts, 'char');  % treat all columns as text
    T = readtable(events_file, opts);
    events = table2struct(T);
    for i = 1:length(events)
        EEG.event(i).latency   = str2double(events(i).onset);
        EEG.event(i).duration  = str2double(events(i).duration);
        if(isfield(events(i),'label'))
            EEG.event(i).type     = normalize_event(properties,events(i).label);
        end
        if(isfield(events(i),'trial_type'))
            EEG.event(i).type     = normalize_event(properties,events(i).trial_type);
        end
        if(isfield(events(i),'sample'))
            EEG.event(i).sample   = str2double(events(i).sample);
        end
    end
end
OutEEG = EEG;

end