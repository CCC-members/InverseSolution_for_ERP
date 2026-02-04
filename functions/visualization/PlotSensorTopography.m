function fig = PlotSensorTopography(Svv, Cdata, titleStr)

    % Load colormap
    cmap = load("tools/mycolormap_brain_basic_conn.mat");

    % Create figure
    fig = figure;
    template = load('axes.mat');
    currentAxes = template.axes;

    % Build electrode structure from Cdata
    Channel = Cdata.Channel;
    Channel(ismember({Channel.Name},{'TP7','T9','FTT9h','T10','FTT10h','TP8'})) = [];
    elec_data = [];
    elec_data.pos = zeros(length(Channel),3);
    for ii = 1:length(Channel)
        elec_data.lbl{ii}   = Channel(ii).Name;
        loc                 = Channel(ii).Loc;
        elec_data.pos(ii,:) = mean(loc,2);
    end
    elec_data.label    = elec_data.lbl;
    elec_data.elecpos  = elec_data.pos;
    elec_data.unit     = 'mm';

    % Normalize diagonal of Svv
    SvvDiag = diag(Svv);
    SvvDiag = abs(SvvDiag) / max(abs(SvvDiag(:)));

    % FieldTrip topo data structure
    topo              = [];
    topo.label        = elec_data.label;
    topo.elec         = elec_data;
    topo.dimord       = 'chan_time';
    topo.time         = 0;              % single snapshot
    topo.avg          = SvvDiag(:);     % one value per channel

    % Config for topoplot
    cfg = [];
    cfg.marker        = 'labels';
    cfg.layout        = []; %'standard_1020.elc';
    cfg.channel       = elec_data.label;
    cfg.markersymbol  = '.';
    cfg.colormap      = cmap.cmap_a;
    cfg.markersize    = 20;
    cfg.markercolor   = 'g';
    cfg.comment       = 'no';
    cfg.interactive   = 'no';
    cfg.colorbar      = 'East';
    cfg.rotate          = 90;
    cfg.parent        = currentAxes;   % plot into loaded axes

    % Topoplot
    ft_topoplotER(cfg, topo);

    % Add title
    sgtitle(titleStr);

end
