function [OutEEGs, errMessage] = SelectConditions(properties,EEG)

disp(strcat("-->> Selecting conditions from EEG data"));
disp('--------------------------------------------------------------------------');
errMessage          = [];
select_list         = properties.event_params.select;
selection           = properties.event_params.selection;
exp_descrip         = properties.event_params.experiment;
type                = selection.type;
trialsPerSet        = selection.nTrial;

count               = 1;

subID               = EEG.subject;
events              = [EEG.event];
if(isempty(fields(events)))
    return;
end
if(~isempty(select_list))
    for l=1:length(select_list)
        event_elem                          = select_list(l);
        if(event_elem.get)
            event_id                        = event_elem.id;
            disp(strcat("---->> Selecting: ",event_id));            
            switch type
                case 'segment'
                    regions                 = get_regions(events,event_id,EEG.pnts,EEG.srate);
                    times                   = [regions.start; regions.end]';
                    for j=1:size(times,1)
                        time                = times(j,:);
                        newEEG              = pop_select(EEG, 'time', time);
                        newEEG.condition    = event_id;
                        newEEG.segment      = num2str(j);
                        OutEEGs(count)      = eeg_checkset(newEEG);
                        count               = count + 1;
                    end

                case 'continue'
                    regions                 = get_regions(events,event_id,EEG.pnts,EEG.srate);
                    times                   = [regions.start; regions.end]';
                    for j=1:size(times,1)
                        time                = times(j,:);
                        newEEG              = pop_select(EEG, 'time', time);
                        newEEG.condition    = event_id;
                        newEEGs(j)          = newEEG;
                    end
                    mergedEEG = pop_mergeset(newEEGs, [1:length(newEEGs)]);
                    mergedEEG.segment   = 1;
                    OutEEGs(count)      = eeg_checkset(mergedEEG);
                    count               = count + 1;
                    
                case 'trials'
                    seg_time = str2double(strrep(selection.time,'s',''));
                    if(isequal(lower(exp_descrip.marks.method),'border'))
                        regions = get_regions_by_borders(properties,events,event_id,select_list,EEG.pnts,EEG.srate);
                    else
                        if(EEG.annotated)
                            regions = get_regions_by_duration(events,event_id);
                        else
                            regions = get_regions(events,event_id,EEG.pnts,EEG.srate);
                        end
                    end
                    if(~isempty(fields(regions)))
                        EEG.Regions = regions;
                        EEG.urevent = [];
                        EEG.event   = [];
                        % try

                        times  = [regions.start; regions.end]';
                        seg = 1;
                        for t=1:size(times,1)
                            time = times(t,:);
                            if (time(2) - time(1)) < seg_time
                                time(2) = time(1) + seg_time;
                            end
                            newEEG = pop_select(EEG, 'time', time);
                            newEEG = eeg_regepochs(newEEG, 'recurrence', seg_time, 'limits', [0 seg_time], 'rmbase', NaN);
                            ntrials = newEEG.trials;
                            for nt = 1:ntrials
                                tmpEEG = pop_select(newEEG, 'trial', nt);
                                tmpEEG.trials = 1;
                                tmpEEG.event   = [];
                                tmpEEG.epoch   = [];
                                tmpEEG.urevent = [];
                                tmpEEG = eeg_checkset(tmpEEG);
                                newEEGs(seg) = tmpEEG;
                                seg = seg + 1;
                            end
                        end
                        if(length(newEEGs) < trialsPerSet)
                            continue;
                        end
                        mergedEEG = pop_mergeset(newEEGs, [1:length(newEEGs)]);
                        mergedEEG.event = [];
                        trialsEEG = eeg_regepochs(mergedEEG, 'recurrence', seg_time, 'limits', [0 seg_time], 'rmbase', NaN);
                        OutEEGs(count) = eeg_checkset(trialsEEG);
                        count = count + 1;
                        clear newEEGs;
                        % catch Ex
                        %     disp('--------------------------------------------------------------------------');
                        %     disp("-->> ERROR");
                        %     disp(Ex.message);
                        %     disp('--------------------------------------------------------------------------');
                        %     continue;
                        % end
                    end
            end
        end
    end
end

if(exist("OutEEGs","var"))
    OutEEGs(cellfun(@isempty, {OutEEGs.filename})) = [];    
    errMessage = [];
else
    OutEEGs = [];
    errMessage = "The subject do not contain any of the selected events";
end

end


%%
%%  Get regions from segments
%%
function regions = get_regions(events,event_id,pnts,srate)
count = 1;
regions = struct;
for j=1:length(events)-1
    if(isequal(events(j).type, event_id))
        regions(count).type     = events(j).type;
        regions(count).start    = events(j).latency;
        regions(count).end      = events(j+1).latency;
        count = count + 1;
    end
end
if(isequal(events(end).type, event_id))
    regions(count).type     = events(end).type;
    regions(count).start    = events(end).latency;
    regions(count).end      = pnts(end)/srate;
end
end

function regions = get_regions_by_duration(events,event_id)
count = 1;
regions = struct;
for j=1:length(events)
    if(isequal(events(j).type, event_id))
        regions(count).type     = events(j).type;
        regions(count).start    = events(j).latency;
        regions(count).end      = events(j).latency + events(j).duration;
        count = count + 1;
    end
end
end


function regions = get_regions_by_borders(properties,events,event_id,select_list,pnts,srate)
unit = properties.event_params.experiment.latency.unit;
count = 1;
regions = struct;
for i=1:length(events)
    if(isequal(events(i).type, event_id))
        regions(count).type     = events(i).type;
        regions(count).start    = events(i).latency;
        if(isequal(lower(unit),'point'))
            regions(count).start = events(i).latency/srate;
        end
        for j=i+1:length(events)
            next_event = events(j);
            names = {select_list.name};
            isMatch = any(cellfun(@(x) contains(next_event.type, x, 'IgnoreCase', true), names ));
            if(isMatch)
                regions(count).end      = next_event.latency;
                if(isequal(lower(unit),'point'))
                    regions(count).end = events(j).latency/srate;
                end
                count = count + 1;
                break;
            end
            if(isequal(j,length(events)))
                regions(count).end      = pnts/srate;
                count = count + 1;
                break;
            end
        end
    end
end
end