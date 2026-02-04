function [Cortex,iHideVert]  = SplitHemisphere(Cortex,hemi)
            [rH, lH, isConnected, iStruct, iRightScout, iLeftScout] = tess_hemisplit(Cortex);
            controls = app.PanelControls.Children;            
            for i=1:length(controls)
                control = controls(i);
                if(isequal('matlab.ui.control.CheckBox',class(control)) && isequal(control.Text,'Left Hemisphere'))
                    hemi(1) = control;
                end
                if(isequal('matlab.ui.control.CheckBox',class(control)) && isequal(control.Text,'Right Hemisphere'))
                    hemi(2) = control;
                end
                if(isequal('matlab.ui.control.CheckBox',class(control)) && isequal(control.Text,'Show Scouts'))
                    hemi(3) = control;
                end
            end
            if(isequal(hemi,'right'))
                iHideVert = lH;
            end
            if(isequal(hemi,'left'))
                iHideVert = rH;
            end
            if(isempty(hemi))
                iHideVert = [];
            end
            
            Cortex.Vertices(iHideVert,:)    = [];
            Cortex.VertConn(iHideVert,:)    = [];
            Cortex.VertConn(:,iHideVert)    = [];
            Cortex.SulciMap(iHideVert,:)    = [];
            Cortex.VertNormals(iHideVert,:) = [];
            Cortex.Curvature(iHideVert,:)   = [];
            isHideFaces                     = any(ismember(Cortex.Faces, iHideVert), 2);
            Cortex.Faces(isHideFaces,:)     = [];
            if(~hemi(1).Value)
                Cortex.Faces                = Cortex.Faces(:,:) - length(iHideVert);
            end

            if(isequal(length(hemi),3) && ~isempty(hemi(3)) && hemi(3).Value)
                sAtlas = Cortex.Atlas(Cortex.iAtlas);
                if isempty(sAtlas)
                    return;
                end
                sScouts = sAtlas.Scouts;
                if isempty(sScouts)
                    return;
                end
                sScouts(ismember({sScouts.Label},{'Thalamus L','Thalamus R','Unknown L','Unknown R','VentralDC L L','VentralDC L R','VentralDC R',...
                    'Ventricle inf-lat L','Ventricle inf-lat R','Ventricle lat L L','Ventricle lat L R','Ventricle lat R',...
                    'White L','White R'})) = [];
                if(~hemi(1).Value)
                    sScouts(ismember({sScouts.Region},{'LU','UU'})) = []; 
                end
                if( ~hemi(2).Value)
                    sScouts(ismember({sScouts.Region},{'RU','UU'})) = [];
                end  
                for i=1:length(sScouts)
                    sScouts(i).Vertices(ismember(sScouts(i).Vertices,iHideVert)) = [];
                    if(~hemi(1).Value)
                        sScouts(i).Vertices = sScouts(i).Vertices - length(iHideVert);                        
                    end
                end
                Cortex.Atlas(Cortex.iAtlas).Scouts = sScouts;
            end
        end       