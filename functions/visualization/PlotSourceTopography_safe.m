function fig = PlotSourceTopography_safe(Cortex,J,cmap_struct,ttl)
    if ~isfield(Cortex,'Faces') || ~isfield(Cortex,'Vertices')
        warning('Cortex lacks Faces/Vertices; skipping PlotSourceTopography.');
        return;
    end
    if ~isfield(Cortex,'SulciMap'), Cortex.SulciMap = zeros(size(Cortex.Vertices,1),1); end
    if ~isfield(Cortex,'VertConn')
        if exist('tess_vertconn','file')==2
            Cortex.VertConn = tess_vertconn(Cortex.Faces, size(Cortex.Vertices,1));
        else
            Cortex.VertConn = []; % no smoothing fallback
        end
    end
    nV = size(Cortex.Vertices,1);
    J = J(:);
    if numel(J) ~= nV
        warning('Map length (%d) != #vertices (%d). Truncating/padding.', numel(J), nV);
        if numel(J) > nV, J = J(1:nV);
        else, J(end+1:nV,1) = 0; end
    end
    fig = PlotSourceTopography(Cortex, J, cmap_struct);
    try, sgtitle(ttl); catch, title(ttl); end
end


function fig = PlotSourceTopography(Cortex,J,colorMap)
    fig = figure;
    try
        template = load('axes.mat');
        if isfield(template,'axes') && isvalid(template.axes)
            currentAxes = template.axes;
            set(currentAxes,"Parent",fig);
        else
            currentAxes = axes('Parent',fig);
        end
    catch
        currentAxes = axes('Parent',fig);
    end

    sources_iv = sqrt(abs(J));
    if max(sources_iv(:))>0, sources_iv = sources_iv / max(sources_iv(:)); end
    smoothValue = 0.66; 
    SurfSmoothIterations = 10;

    if isfield(Cortex,'VertConn') && ~isempty(Cortex.VertConn) && exist('tess_smooth','file')==2
        Vertices = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);
    else
        Vertices = Cortex.Vertices;
    end
    if ~isfield(Cortex,'SulciMap') || isempty(Cortex.SulciMap)
        sulci = 0;
    else
        sulci = Cortex.SulciMap*0.06;
        if numel(sulci) ~= size(Vertices,1), sulci = 0; end
    end

    patch(currentAxes, ...
        'Faces',Cortex.Faces, ...
        'Vertices',Vertices, ...
        'FaceVertexCData',sulci + log(1+sources_iv), ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'AlphaDataMapping', 'none', ...
        'EdgeAlpha',        1, ...
        'BackfaceLighting', 'lit', ...
        'AmbientStrength',  0.5, ...
        'DiffuseStrength',  0.5, ...
        'SpecularStrength', 0.2, ...
        'SpecularExponent', 1, ...
        'SpecularColorReflectance', 0.5, ...
        'FaceLighting',     'gouraud', ...
        'EdgeLighting',     'gouraud', ...
        'FaceAlpha',1);
    set(currentAxes,'xcolor','w','ycolor','w','zcolor','w');
    view(currentAxes,0,0);
    if isstruct(colorMap) && isfield(colorMap,'cmap_a')
        colormap(currentAxes,colorMap.cmap_a);
    else
        colormap(currentAxes,parula(256));
    end
    rotate3d(currentAxes,'on');
end
