function eventOut = normalize_event(properties,event_label)

eventOut = [];
event_list = properties.event_params.select;

for e=1:length(event_list)
    event = event_list(e);
    if(any(contains(event.name, event_label, 'IgnoreCase', true)))
        eventOut = event.id;
        break;
    end
end
if(isempty(eventOut))
    eventOut = event_label;
    missing_vents = jsondecode(fileread(strcat('missing_events.json')));
    if(isempty(find(ismember(missing_vents.events,event_label),1)))
        if(isempty(missing_vents.events))
            missing_vents.events{1} = event_label;
        else
            missing_vents.events{end+1} = event_label;
        end
        saveJSON(missing_vents,strcat('functions/preprocessing/missing_events.json'));
    end
else

end

