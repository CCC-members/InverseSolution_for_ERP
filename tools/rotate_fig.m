function rotate_fig(rotate)
switch lower(rotate)
    case {'front','f'}
        view(90,0)
    case {'back','b'}
        view(270,0)
    case {'left','l'}
        view(180,0)
    case {'right','r'}
        view(0,0)
    case {'frontright','fr'}
        view(135,30)
    case {'backright','br'}
        view(45,30)
    case {'frontleft','fl'}
        view(-135,30)
    case {'backleft','bl'}
        view(-45,30)
    case 'top'
        view(180,90)
    case 'bottom'    % undocumented option!
        view(0,-90)
end
end

