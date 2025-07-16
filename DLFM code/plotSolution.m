clc;
clear;
close all;
addpath(genpath('./src'));

%% Simulation parameters
runparameters_Fig2C;

%% Loop over all cases
for caseind = 1:Ncases

    frameNum         = framemat(caseind);
    parentfolder     = sprintf(parentfolderformat,cellIDmat(caseind));
    framefolder      = sprintf(['%s/',framefolderformat],parentfolder,frameNum);
    phaseimagefile   = sprintf(['%s/',phaseimagefileformat],framefolder,frameNum);
    fluoreImageFile  = sprintf(['%s/',fluoroimagefileformat],framefolder,frameNum);
    defcoordsfile    = sprintf(['%s/',nodesdatafileformat],framefolder,frameNum);
    
    fprintf(1,'Processing frame %d of %d\n',caseind, Ncases);

    % Energy term weights
    wt_nodalerror    = wtNmat(caseind);
    wt_line          = wtLmat(caseind);
    wt_el_pen        = wtEpmat(caseind);
    
    E_actctr         = E_actctrmat(caseind);
    dia_actctr       = dia_actctrmat(caseind);
    
    if ~usemaskbasedlimits
        xmin = xminmat(caseind);
        xmax = xmaxmat(caseind);
        ymin = yminmat(caseind);
        ymax = ymaxmat(caseind);
    end
    
    % Axis position (Maximize Plot area as we don't have axis labels or
    % figure title)
    axisposition = [0.05 0.05 0.9 0.9];

    % Creating Folder for saving tension figures
    TensionFigures      = sprintf('%s/Tension_Fig',parentfolder);
    TensionFluroFigures = sprintf('%s/TensionFluro_Fig',parentfolder);
    DisplacementFigures = sprintf('%s/Displacement_Fig',parentfolder);
    PhaseImageFigures   = sprintf('%s/PhaseImage_Fig',parentfolder);
    if ~isfolder(TensionFigures) 
        mkdir(TensionFigures) 
    end
    if ~isfolder(TensionFluroFigures) 
        mkdir(TensionFluroFigures) 
    end
    if ~isfolder(DisplacementFigures) 
        mkdir(DisplacementFigures) 
    end
    if ~isfolder(PhaseImageFigures) 
        mkdir(PhaseImageFigures) 
    end
    
    for regind = length(regcoeffmat)
        
        regcoeff   = regcoeffmat(regind);
        fprintf(1,'Plotting for reg parameter = %3.2E\n',regcoeff);
        outputfoldername      = sprintf('%s/frame%02d/withAC_regcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,parentfolder,frameNum,regcoeff,E_actctr,dia_actctr,wt_nodalerror,wt_line,wt_el_pen);
        
        cellmaskfile = sprintf('%s/%02d_cellROI.mat',framefolder,frameNum);
        cellmask = load(cellmaskfile);
        maskcoords = cellmask.hroipoly.Position;

        imageresolutionfile  = sprintf('%s/pix2um.txt',framefolder);
        pix2um               = load(imageresolutionfile);
        
        if usemaskbasedlimits == 1
            maskxmin = min(maskcoords(:,1));
            maskxmax = max(maskcoords(:,1));
            maskymin = min(maskcoords(:,2));
            maskymax = max(maskcoords(:,2));
            
            maskpolyin = polyshape(maskcoords(:,1),maskcoords(:,2));
            [cellcenX,cellcenY] = centroid(maskpolyin);
            
            xmin = cellcenX - framesize_x/2;
            xmax = cellcenX + framesize_x/2;
            ymin = cellcenY - framesize_y/2;
            ymax = cellcenY + framesize_y/2;
        end
        
        if forwardproblem == 1
            solfilename = sprintf('%s/ForwardProbSolution.mat',FWDproblemparentfolder);
            rundatafilename = sprintf('%s/rundata.mat',FWDproblemparentfolder);
        else
            solfilename = sprintf('%s/InverseProbSolution_reg%2.4E.mat',outputfoldername,regcoeff);
            beammeshfilename = sprintf('%s/beammesh.mat',framefolder);
            rundatafilename = sprintf('%s/rundata.mat',outputfoldername);
            datax = load(sprintf('%s/n2_uxm.txt',framefolder));
            datay = load(sprintf('%s/n2_uym.txt',framefolder));
        end
        
        if exist(solfilename)
            load(solfilename);
            load(beammeshfilename);
        elseif exist(rundatafilename)
            load(rundatafilename);
            fprintf(1,'Using rundata %s\n',rundatafilename);
        else
            error('Could not find the solution file at\n%s',solfilename);
        end
        
        beammesh      = DOFhandlers(1);
        substratemesh = DOFhandlers(2);

        unodal = reshape(beammesh.u(:),3,beammesh.Nnodes)';
        currpos   = beammesh.coords + dispampfactor*unodal(:,1:2);
        refpos    = beammesh.coords;
        
        %% Initialize plots
        if plotsigma == 1
            fh_sigma   = figure();
            fh_sigmaFl = figure();
        else
            fh_sigma = -1;
            fh_sigmaFl = -1;
        end
        if plotcontractility == 1
            fh_C = figure();
            polarscatter(pi,0,'k.');
        else
            fh_C = -1;
        end
        
        %% Displacement plot
        if plotdisps == 1
            fh_disp = figure(); hold on;
            figure(fh_disp);
            if disp_backimage == 1
                cellimagefile = fluoreImageFile;
            else
                cellimagefile = phaseimagefile;
            end
            if isfile(cellimagefile)
                cellimg = double(imread(cellimagefile));
                hcimg = imshow(cellimagefile);
                hcimg.YData = hcimg.YData*pix2um;
                hcimg.XData = hcimg.XData*pix2um;
                hadisp = findall(fh_disp, 'type', 'axes');
                hadisp.Visible = 'On';
            else
                error('Did not find %s\n', cellimagefile);
            end
            
            dispnodalb  = reshape(beammesh.u(:),3,beammesh.Nnodes)';
            defcoords   = beammesh.coords + dispampfactor*dispnodalb(:,1:2);
            refcoords   = beammesh.coords;
            
            dispnodals   = reshape(substratemesh.u(:),3,substratemesh.Nnodes)';
            dispnodal    = [dispnodalb;dispnodals];
            defcoords   = [defcoords; substratemesh.coords + dispampfactor*dispnodals(:,1:2)];
            refcoords   = [refcoords; substratemesh.coords];
            
            figure(fh_disp);hold on;
            Ninter = beammeshdata.nGridIntersectNodes;
            quiver_tri(refcoords(1:Ninter,1),refcoords(1:Ninter,2),...
                dispnodal(1:Ninter,1)*dispscalefactor,dispnodal(1:Ninter,2)*dispscalefactor...
                ,disp_arrowheadsize,22.5,disp_arrowlinewidth,disp_arrowcolor);
            
            % Quiver plot
            %figure(fh_disp);hold on;
            %Ninter = beammeshdata.nGridIntersectNodes;
            %qut = quiver(refcoords(1:Ninter,1),refcoords(1:Ninter,2),dispnodal(1:Ninter,1)*dispscalefactor,dispnodal(1:Ninter,2)*dispscalefactor);
            %qut.LineWidth       = 2.0;
            %qut.MaxHeadSize     = 3.0;
            %qut.Color           = 'r';
            %qut.AutoScale       = 'off';
            %plot(refcoords(1:Ninter,1),refcoords(1:Ninter,2),'ko');
            %axis equal
            
            Umax = max(sqrt(dispnodal(1:Ninter,1).^2 + dispnodal(1:Ninter,2).^2));
            
        end
        
        %% Background image (phase) and mesh plot
        if plotsigma == 1
            figure(fh_sigma); hold on;

            cellimagefile = phaseimagefile;
            if isfile(cellimagefile)
                cellimg = imread(cellimagefile);

                % Reduce brighness for arrow visibility
                cellimghsv = rgb2hsv(cellimg);
                cellimghsv(:,:,3) = cellimghsv(:,:,3)/(mean(mean(cellimghsv(:,:,3))))*0.55;
                cellimgconv = hsv2rgb(cellimghsv);
                
                hold on;
                hcimg = imshow(cellimgconv);
                hcimg.YData = hcimg.YData*pix2um;
                hcimg.XData = hcimg.XData*pix2um;
            else
                error('Did not find %s\n', cellimagefile);
            end
            ha = findall(fh_sigma, 'type', 'axes');
            ha.Visible = 'On';

            % Set to 'reverse' when plotting image on top
            set(ha, 'YDir','reverse');

            % add mesh
            if plotnanonet ==1
                % Generate smooth deformed configuration using Hermite interpolation
                [xelems_intp,yelems_intp] = getHermiteInterpolatedCoords(beammesh,3);
                figure(fh_sigma); hold on;
                hdef = plot(xelems_intp',yelems_intp','-','LineWidth',1.0,'color','k');
            end
        end
        
        %% Background image (fluoro) and mesh plot
        if plotsigma == 1
            figure(fh_sigmaFl); hold on;

            cellimagefile = fluoreImageFile;
            if isfile(cellimagefile)
                cellimg = imread(cellimagefile);

                hold on;
                hcimg = imshow(cellimg);
                hcimg.YData = hcimg.YData*pix2um;
                hcimg.XData = hcimg.XData*pix2um;
            else
                error('Did not find %s\n', cellimagefile);
            end
            ha = findall(fh_sigma, 'type', 'axes');
            ha.Visible = 'On';

            % Set to 'reverse' when plotting image on top
            set(ha, 'YDir','reverse');

            % add mesh
            if plotnanonet ==1
                % Generate smooth deformed configuration using Hermite interpolation
                [xelems_intp,yelems_intp] = getHermiteInterpolatedCoords(beammesh,3);
                figure(fh_sigmaFl); hold on;
                hdef = plot(xelems_intp',yelems_intp','-','LineWidth',1.0,'color','k');
            end

            % add target nodal intersections
            % Read measured data (dof wise)
            mDOFs   = [];
            % x-displacements
            framefolder      = sprintf('%s/frame%02d',parentfolder,frameNum); % input data
            datax = load(sprintf('%s/n2_uxm.txt',framefolder));
            % y-displacements
            datay = load(sprintf('%s/n2_uym.txt',framefolder));

            measureddefcoords = zeros(length(datax));
            measureddefcoords(:,1) = datax(:,2) + refpos(datax(:,1),1);
            measureddefcoords(:,2) = datay(:,2) + refpos(datay(:,1),2);
            figure(fh_sigmaFl); hold on;
            plot(measureddefcoords(:,1),measureddefcoords(:,2),'ko','MarkerSize',5,'MarkerFaceColor','r');

            %         subdefcoordsfilename = sprintf('%s/subdef_nodes_%02d.txt',datafoldername,frameNum);
            %         if isfile(subdefcoordsfilename)
            %             subdefcoords = load(subdefcoordsfilename);
            %             plot(subdefcoords(:,1),subdefcoords(:,2),'ko','MarkerSize',5,'MarkerFaceColor','r');
            %         end
        end
        %% Get nodal forces
        nodalRHS = reshape(HiperProblem.Fnodal(1:DOFhandlers(1).Ndofs),3,beammesh.Nnodes)';
        nodalForces = nodalRHS(:,1:2);
        
        %% Plot force perunit length at element centers
        if plotsigma == 1 && plotlinetension == 1
            
            sigmamat = zeros(beammesh.Nelems,2);
            ndofs    = beammesh.dofs;
            connecind = zeros(beammesh.Nnodes,1);
            for nodeID = 1:beammesh.Nnodes
                connecind(nodeID) = numel(find(beammesh.connec(:)==nodeID));
            end
            
            for elemID = 1:beammesh.Nelems
                node1 = beammesh.connec(elemID,1);
                node2 = beammesh.connec(elemID,2);
                
                % Get nodal coords and current displacements
                X1 = beammesh.coords(node1,1);
                Y1 = beammesh.coords(node1,2);
                X2 = beammesh.coords(node2,1);
                Y2 = beammesh.coords(node2,2);
                
                u1 = beammesh.u(ndofs*(node1 - 1)  + 1);
                w1 = beammesh.u(ndofs*(node1 - 1)  + 2);
                
                u2 = beammesh.u(ndofs*(node2 - 1)  + 1);
                w2 = beammesh.u(ndofs*(node2 - 1)  + 2);
                
                L0  = sqrt((X2 - X1)^2 + (Y2 - Y1)^2);
                L   = sqrt(((X2 + u2) - (X1 + u1))^2 + ((Y2 + w2) - (Y1 + w1))^2);
                
                sigmamat(elemID,1) = nodalForces(node1,1)/L/connecind(node1) + nodalForces(node2,1)/L/connecind(node2);
                sigmamat(elemID,2) = nodalForces(node1,2)/L/connecind(node1) + nodalForces(node2,2)/L/connecind(node2);
            end
            
            elemcencoords = 0.5*((currpos(beammesh.connec(:,1),1:2)) + (currpos(beammesh.connec(:,2),1:2)));
            
            sigma_mag = sqrt(sigmamat(:,1).^2 + sigmamat(:,2).^2);
            sigmax = max(sigma_mag);
            
            forcenodes = sigma_mag > 1e-2;
            figure(fh_sigma); hold on;
            quiver_tri(elemcencoords(forcenodes,1),elemcencoords(forcenodes,2),...
                sigmamat(forcenodes,1)*sigmascalefactor,sigmamat(forcenodes,2)*sigmascalefactor...
                ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,sigma_arrowcolor,...
                sigma_quifactor);

            figure(fh_sigmaFl); hold on;
            quiver_tri(elemcencoords(forcenodes,1),elemcencoords(forcenodes,2),...
                sigmamat(forcenodes,1)*sigmascalefactor,sigmamat(forcenodes,2)*sigmascalefactor...
                ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,sigma_arrowcolor,...
                sigma_quifactor);
            
            %% Contour plots of displacement and line tensions
            if plotcontourmaps == 1
                xv = xmin:1:xmax;
                yv = ymin:1:ymax;
                [X,Y] = meshgrid(xv,yv);
                sigX = griddata(elemcencoords(:,1),elemcencoords(:,2),sigmamat(:,1),X,Y,'cubic');
                sigY = griddata(elemcencoords(:,1),elemcencoords(:,2),sigmamat(:,2),X,Y,'cubic');
                sigX(isnan(sigX)) = 0;
                sigY(isnan(sigY)) = 0;
                
                redcmap = flip(cbrewer('div','RdBu', 41));
                
                siglim = max([max(max(abs(sigY))),max(max(abs(sigX)))]);
                minlim  = -siglim*1.25;
                maxlim  = siglim*1.25;
                fh_contoursigY = plotcontourdata(xv,yv,sigY,minlim,maxlim,redcmap);
                fh_contoursigX = plotcontourdata(xv,yv,sigX,minlim,maxlim,redcmap);
                Ninter = beammeshdata.nGridIntersectNodes;
                UX = griddata(refcoords(1:Ninter,1),refcoords(1:Ninter,2),dispnodal(1:Ninter,1),X,Y,'cubic');
                UY = griddata(refcoords(1:Ninter,1),refcoords(1:Ninter,2),dispnodal(1:Ninter,2),X,Y,'cubic');
                
                redcmap = flip(cbrewer('div','PiYG', 41));
                Ulim = max([max(max(abs(UY))),max(max(abs(UX)))]);
                minlim  = -Ulim*1.25;
                maxlim  = Ulim*1.25;
                fh_contourUX = plotcontourdata(xv,yv,UX,minlim,maxlim,redcmap);
                
                UYlim = max(max(abs(UY)));
                fh_contourUY = plotcontourdata(xv,yv,UY,minlim,maxlim,redcmap);
            end
        end
     
        %% Calculate contractility value
        % Get epicenter for the force distribution
        [coord_epicen,contractility,thCmax,Cmax,unitvecangle,contractility_theta]...
            = getContractilityEpicenter(nodalForces,currpos);
        fprintf(1,'Contractility = %3.2f nN, P = %2.3E\n',contractility,Cmax/contractility);

        if plotcontractility == 1
            figure(fh_C); hold on;
            polarplot(unitvecangle,contractility_theta,'color',sigma_arrowcolor,'linewidth',4.0);
            % pax = gca;
            % pax.ThetaAxisUnits = 'radians';
            % pax.RLim = [0,ceil(Cmax/contractility)];
            % pax.RAxisLocation = 0;
            % pax.GridLineStyle = '--';
            % pax.ThetaDir = 'clockwise'; % To match with the x-y orientation with flipped Yaxis
            % title(sprintf('$C$=%2.1f nN',contractility))
        end
        
        %% Plot nodal forces
        if plotnodalforce == 1
            
            Fmag = sqrt(nodalForces(:,1).^2 + nodalForces(:,2).^2);
            Fmax = max(sqrt(nodalForces(:,1).^2 + nodalForces(:,2).^2));
            
            figure(fh_sigma);hold on;
            forcenodes = Fmag > 1e-2;
            quiver_tri(currpos(forcenodes,1),currpos(forcenodes,2),...
                nodalForces(forcenodes,1)*fscalefactor,nodalForces(forcenodes,2)*fscalefactor...
                ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,'r',sigma_quifactor);
        end
        
        
        %% Adjust figure size and add scale bars

        %% sigma with phase back image
        if plotsigma == 1
        figure(fh_sigma);
        fh_sigma.Color = 'w';
        set(gcf,'InvertHardCopy','off');% Added to make the scale bar appear white
        ha = findall(fh_sigma, 'type', 'axes');
        box off
        axis equal
        axis off
        ha.Position = axisposition;
        ha.XLim     =  [xmin xmax];
        ha.YLim     =  [ymin ymax];
        
        fh_sigma.Units         = 'centimeters';
        fh_sigma.Position      = [1 1 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
        
        % Add tension and force scale bars
        % TO DO: Save a text file with all scale bars
        if plotlinetension == 1
            quiver_tri(xmax - 35,ymax-5,sigmascalefactor*sigma_scalebarval,0,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,...
                sigma_arrowcolor,sigma_quifactor);
            
            if disp_showarrowtext == 1
                text(xmax - 35,(ymax-5)-6,sprintf('%d mN/m',sigma_scalebarval),'Color',sigma_arrowcolor,'FontSize',sigma_scalebar_textsize);
            end
            
        elseif plotnodalforce == 1
            quiver_tri(xmax - 55,ymax-5,fscalefactor*250,0,arrowsize,22.5,arrowlinewidth,arrowcolor);
        end
        
        % Add lengthscale scale bar for tension figure
        LscaleArrow            = annotation('line') ;
        LscaleArrow.Parent     = ha;  % or any other existing axes or figure
        LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
        LscaleArrow.Color      = lenscalebar_color;
        LscaleArrow.LineWidth  =  lenscalebar_linewidth;

        %% sigma with fluoro back image
        figure(fh_sigmaFl);
        fh_sigmaFl.Color = 'w';
        set(gcf,'InvertHardCopy','off');% Added to make the scale bar appear white
        ha = findall(fh_sigmaFl, 'type', 'axes');
        box off
        axis equal
        axis off
        ha.Position = axisposition;
        ha.XLim     =  [xmin xmax];
        ha.YLim     =  [ymin ymax];
        
        fh_sigmaFl.Units         = 'centimeters';
        fh_sigmaFl.Position      = [1 1 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
        
        % Add tension and force scale bars
        % TO DO: Save a text file with all scale bars
        if plotlinetension == 1
            quiver_tri(xmax - 35,ymax-5,sigmascalefactor*sigma_scalebarval,0,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,...
                sigma_arrowcolor,sigma_quifactor);
            
            if disp_showarrowtext == 1
                text(xmax - 35,(ymax-5)-6,sprintf('%d mN/m',sigma_scalebarval),'Color',sigma_arrowcolor,'FontSize',sigma_scalebar_textsize);
            end
            
        elseif plotnodalforce == 1
            quiver_tri(xmax - 55,ymax-5,fscalefactor*250,0,arrowsize,22.5,arrowlinewidth,arrowcolor);
        end
        
        % Add lengthscale scale bar for tension figure
        LscaleArrow            = annotation('line') ;
        LscaleArrow.Parent     = ha;  % or any other existing axes or figure
        LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
        LscaleArrow.Color      = lenscalebar_color;
        LscaleArrow.LineWidth  =  lenscalebar_linewidth;
        end
        %% Contractility plot
        if plotcontractility == 1
            figure(fh_C);
            fh_C.Color         = 'w';
            fh_C.Units         = 'centimeters';
            fh_C.Position      = [16 2 12 12];

            Cmax = 3000; % [nN]
            pax = gca;

            pax.RLim = [0, Cmax];
            pax.RTick = 0:500:Cmax;

            pax.ThetaAxisUnits = 'radians';
            pax.ThetaDir = 'clockwise'; % To match with the x-y orientation with flipped Yaxis

            pax.ThetaTick = [0 pi/2 pi 3*pi/2];

            pax.RAxisLocation = 3*pi/2;
            pax.GridLineStyle = '-';

            pax.ThetaTickLabel = {};
            pax.FontSize = 18;
            pax.Layer = 'top';
            
        end
        
        %% Displacements plot
        if plotdisps == 1
            figure(fh_disp);
            fh_disp.Color = 'w';
            set(gcf,'InvertHardCopy','off');
            ha = findall(fh_disp, 'type', 'axes');
            box off
            axis equal
            axis off
            ha.Position = axisposition;
            ha.XLim     =  [xmin xmax];
            ha.YLim     =  [ymin ymax];
            
            fh_disp.Units         = 'centimeters';
            fh_disp.Position      = [9 1 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
            
            % Add scale bar
            quiver_tri(xmax - 35,ymax-5,dispscalefactor*disp_scalebarval,0,...
                disp_arrowheadsize,22.5,disp_arrowlinewidth,disp_arrowcolor);
            
            if disp_showarrowtext == 1
                text(xmax - 35,(ymax-5)-8,sprintf('$%3.1f \\mu m$',disp_scalebarval),'Color',disp_arrowcolor,'FontSize',disp_scalebar_textsize)
            end
            
            LscaleArrow            = annotation('line') ;
            LscaleArrow.Parent     = ha;  % or any other existing axes or figure
            LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
            LscaleArrow.Color      = lenscalebar_color;
            LscaleArrow.LineWidth  =  lenscalebar_linewidth;
            
        end
        
        %% Plot Contour plots of tractions and disps
        if plotcontourmaps == 1
            figurehmat = [fh_contoursigX,fh_contoursigY,fh_contourUX,fh_contourUY];
            for fh = figurehmat
                figure(fh);
                ha = findall(fh, 'type', 'axes');
                box off
                axis equal
                axis off
                ha.Position = axisposition;
                ha.XLim     =  [xmin xmax];
                ha.YLim     =  [ymin ymax];
                
                fh.Units         = 'centimeters';
                fh.Position      = [9 9 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
                
                cb = colorbar('SouthOutside');
                cb.Position = [0.7 0.075 0.25 0.02];
            end
        end
        
        %% plot and save phase and fluorescent image in the same format
        if plotphaseimage == 1
            fh_phase = figure();
            figure(fh_phase);
            set(gcf,'InvertHardCopy','off');
            cellimagefile = phaseimagefile;
            if isfile(cellimagefile)
                cellimg = double(imread(cellimagefile));
                hcimg = imshow(cellimagefile);
                hcimg.YData = hcimg.YData*pix2um;
                hcimg.XData = hcimg.XData*pix2um;
            end
            fh_phase.Color = 'w';
            ha = findall(fh_phase, 'type', 'axes');
            box off
            axis equal
            axis off
            ha.Position = axisposition;
            ha.XLim     = [xmin xmax];
            ha.YLim     = [ymin ymax];
            
            fh_phase.Units         = 'centimeters';
            fh_phase.Position      = [12 2 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
            
            LscaleArrow            = annotation('line') ;
            LscaleArrow.Parent     = ha;  % or any other existing axes or figure
            LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
            LscaleArrow.Color      = lenscalebar_color;
            LscaleArrow.LineWidth  = lenscalebar_linewidth;
        end

        %% Plot contractile epicenter
        %         figure(fh_comb);hold on;
        %         plot(coord_epicen(1),coord_epicen(2),'ko','MarkerSize',5,'MarkerFaceColor','g');
        
        %% Save images
        fprintf(1,'Saving figures...\n');
        imresolution = '-r300'; % resolution at which to save figures
        saveFiguresFolder = sprintf('%s/Figures/frame%d',parentfolder,frameNum);
        if exist(saveFiguresFolder,'dir')~=7
            mkdir(saveFiguresFolder) ;
        end
        
        if plotsigma == 1
            if plotnodalforce == 1
                print(fh_sigma,sprintf('%s/%02d_Forces_frame%d_Reg%2.3E.png',saveFiguresFolder,regind, frameNum,regcoeff),'-dpng',imresolution);
            elseif plotlinetension == 1
                print(fh_sigma,sprintf('%s/frame%04d.png',TensionFigures,frameNum),'-dpng',imresolution);
                %print(fh_sigmaFl,sprintf('%s/frame%04d.png',TensionFluroFigures,frameNum),'-dpng',imresolution);
                print(fh_sigmaFl,sprintf('%s/frame%04d.jpg',TensionFluroFigures,frameNum),'-djpeg',imresolution);
            else
                error('Specify correct plot option\n');
            end
        end
        
        if plotcontractility == 1
            print(fh_C,sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.png',saveFiguresFolder,regind,frameNum,regcoeff),'-dpng',imresolution);
            print(fh_C,sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.eps',saveFiguresFolder,regind,frameNum,regcoeff),'-depsc');
        end
        % save textfile to plot data seperately
        fp_cdata = fopen(sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.txt',saveFiguresFolder,regind,frameNum,regcoeff),'w+');
        fprintf(fp_cdata,'Theta  C_theta  C\n');
        fprintf(fp_cdata,'%2.3E %2.3E %2.3E\n',[unitvecangle(:),contractility_theta(:),contractility*ones(length(unitvecangle),1)]');
        fclose(fp_cdata);
        
        if plotdisps == 1
           print(fh_disp,sprintf('%s/%02d_dispplot_frame%04d.png',saveFiguresFolder,regind,frameNum),'-dpng',imresolution); 
           print(fh_disp,sprintf('%s/frame%04d.png',DisplacementFigures,frameNum),'-dpng',imresolution);
        end
        
        if plotphaseimage == 1
           print(fh_phase,sprintf('%s/%02d_phaseimage_frame%04d.png',saveFiguresFolder,regind,frameNum),'-dpng',imresolution); 
           print(fh_phase,sprintf('%s/frame%04d.png',PhaseImageFigures,frameNum),'-dpng',imresolution);
        end
        
        if plotcontourmaps == 1
            figurehmat = [fh_contoursigX,fh_contoursigY,fh_contourUX,fh_contourUY];
            figurename = {'sigX','sigY','Ux','Uy'};
            for i = 1:length(figurehmat)
                print(figurehmat(i),sprintf('%s/%02d_%s_frame%d.png',saveFiguresFolder,regind,figurename{i},frameNum),'-dpng',imresolution); 
            end
        end
        
        localframenum = 1;
        %% For calculating local force
        if plotlocalforce ==1
            fprintf(1, 'Select a polygonal region to evaluate local forces...\n' )
            % [al1,al2,al3,al4] = measure_local_force(fh_comb,beammesh);
            figure(fh_sigma);
            hroirect          = drawpolygon('StripeColor','y'); waitfordoubleclick;
            measurerectcoords = hroirect.Position;
            selectedmNodes = find(inpolygon(beammesh.coords(:,1),beammesh.coords(:,2),...
                measurerectcoords(:,1),measurerectcoords(:,2)));
            xmin_local = min(measurerectcoords(:,1));
            xmax_local = max(measurerectcoords(:,1));
            
            ymin_local = min(measurerectcoords(:,2));
            ymax_local = max(measurerectcoords(:,2));
            
            Fx = nodalForces(selectedmNodes,1);
            Fy = nodalForces(selectedmNodes,2);
            ForceX = sum(Fx);
            ForceY = sum(Fy);
            Force = sqrt(ForceX^2+ ForceY^2);
            
            polyin = polyshape(measurerectcoords(:,1),measurerectcoords(:,2));
            [x,y] = centroid(polyin);
            
            fscalefactor_local = 8/(Force);
            
            % Plot the local force
            fh_localforce = figure(); hold on;
            set(gcf,'InvertHardCopy','off');
            cellimagefile = phaseimagefile;
            
            cellimg = imread(cellimagefile);
            hold on;
            hcimg = imshow(cellimg);
            hcimg.YData = hcimg.YData*pix2um;
            hcimg.XData = hcimg.XData*pix2um;
            
            fh_localforce.Color = 'w';
            ha = findall(fh_localforce, 'type', 'axes');
            box off
            axis equal
            axis off
            ha.Position = axisposition;
            ha.XLim     =  [xmin_local-0.05*xmin_local xmax_local+0.05*xmax_local];
            ha.YLim     =  [ymin_local-0.05*ymin_local ymax_local+0.05*ymax_local];
            
            
            
            localF_arrowlinewidth = 5.0;
            localF_max            = 250;    
            localF_scalefactor    = 10/localF_max;   % um/(nN)
            localF_arrowheadsize  = 0.3*localF_scalefactor*localF_max;
            localF_arrowcolor     = 'r';
            
            f = 1:length(measurerectcoords);
            patch('Faces',f,'Vertices',measurerectcoords,...
                'EdgeColor',localF_arrowcolor,'FaceColor',[0.0 0.7 1.0],'FaceAlpha',.2,'LineStyle','--','LineWidth',2);
            quiver_tri(x,y,...
                ForceX*fscalefactor_local,ForceY*fscalefactor_local...
                ,localF_arrowheadsize,22.5,localF_arrowlinewidth,'r');
            titlestr = sprintf(' F = %5.3f nN ', Force);
            ht = title(titlestr, 'interpreter', 'latex','FontSize',18);
            
            LscaleArrow            = annotation('line') ;
            LscaleArrow.Parent     = ha;  % or any other existing axes or figure
            LscaleArrow.Position   = [xmin_local-0.04*xmin_local ymax_local+0.04*ymax_local 10 0];
            LscaleArrow.Color      = lenscalebar_color;
            LscaleArrow.LineWidth  =  5;
            
            
            print(fh_localforce,sprintf('%s/frame%d_local_%d.png',TensionFigures,frameNum,localframenum),'-dpng',imresolution);
            localframenum = localframenum + 1;
            
        end
        
       scalefile = sprintf('%s/scale.txt', TensionFigures);
       fileID = fopen(scalefile,'w');
       fprintf(fileID,'Length scale: %d um \n',lenscalebarval);
       fprintf(fileID,'Tension scale: %d mN/m \n',sigma_max);
       fclose(fileID);

       if plotdisps == 1
           scalefile = sprintf('%s/scale.txt', DisplacementFigures);
           fileID = fopen(scalefile,'w');
           fprintf(fileID,'Length scale: %d um \n',lenscalebarval);
           fprintf(fileID,'Disp  scale: %d um \n',disp_max);
           fclose(fileID);
       end
        
            
        %% close figures
        if closefigures == 1
            if plotsigma == 1
                close(fh_sigma);
                close(fh_sigmaFl);
            end
            if plotdisps == 1
                close(fh_disp);
            end
            if plotcontractility == 1
                close(fh_C);
            end
            if plotphaseimage == 1
                close(fh_phase);
            end
        end
        
    end
end


%% Movie command
% ffmpeg -f image2 -r 4 -i ./frame%d.png -crf 25 ./celldiv.wmv
% (Better quality)
% ffmpeg -r 2 -f image2 -i frame%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p combined_movie.mp4

%% Plot fluorescent error function with mesh

% fluroimagefile = sprintf('%s/real_%02d_corrected_image.jpg', outputfoldername, frameNum);
% fluroimagefile = sprintf('%s/real_%02d_corrected_image_nodots.png', datafoldername, frameNum);
% grayimg = rgb2gray(imread(fluroimagefile));
% Idata   = 1 - im2double(grayimg);
%
% I_th  = double(Idata > 0.6);
% Eline = imgaussfilt(I_th,8,'FilterSize',41);
%
% h = fspecial('disk',10);
% % h = fspecial('log',10,1.4);
% Idata_filt = imfilter(Idata,h,'replicate');
% % Idata_filt = Idata_filt.^4;
%
% fh_comb = figure(); hold on;
% figure(fh_comb);
% pix2um = 0.32;
% hcimg = imshow(Eline.^2);
% hcimg.YData = hcimg.YData*pix2um;
% hcimg.XData = hcimg.XData*pix2um;
% ha = findall(fh_comb, 'type', 'axes');
% ha.Visible = 'On';
%
% set(ha, 'YDir','reverse')
%
% hold on;
% hdef = plot(xelems',yelems','-','LineWidth',1.0,'color','b');

%% Generating an artificial force patten
% nodalForces_mag = sqrt(nodalForces(:,1).^2 + nodalForces(:,2).^2);
% coord_epicen = [335,300];
% relpos       = currpos - coord_epicen;
% nodalF_angle = atan2(relpos(:,2),relpos(:,1));
% reldist = vecnorm(relpos,2,2);
% nodalForces_mag_new = (nodalForces_mag>0).*(1 - exp(-reldist/75))*75;
% %nodalForces_mag_new = nodalForces_mag_new .* (reldist < 17.5);
%
% nodalForces_new = -[nodalForces_mag_new.*cos(nodalF_angle),nodalForces_mag_new.*sin(nodalF_angle)];
%
%
% figure(fh_comb);hold on;
% qft = quiver(currpos(:,1),currpos(:,2),nodalForces_new(:,1)*fscalefactor,nodalForces_new(:,2)*fscalefactor);
% qft.LineWidth       = 2.0;
% qft.MaxHeadSize     = 3.0;
% qft.Color           = 'k';
% qft.AutoScale       = 'off';
%
% forcevalX = nodalForces_new(:,1);
% forcevalY = nodalForces_new(:,2);
% forceind = find(forcevalX);
% % Assign forces to forcedofs
% forcedofs = [];
% for i = 1:length(forceind)
% forcedofs = [forcedofs;(forceind(i)-1)*3+1, forcevalX(forceind(i))];
% forcedofs = [forcedofs;(forceind(i)-1)*3+2, forcevalY(forceind(i))];
% end
%
% outputfoldername = sprintf('/home/sohan/ProjectData/nanonettfm/ValidationTesting/cell3-Cpattern');
% fp  = fopen(sprintf('%s/dofExternalForces.txt',outputfoldername),'w+');
% for i = 1:size(forcedofs,1)
% fprintf(fp,'%4d %2.5E\n',forcedofs(i,1),forcedofs(i,2));
% end
% fclose(fp);

% Unused mesh operations
%         xelems  = zeros(beammesh.Nelems,2);
%         yelems  = zeros(beammesh.Nelems,2);
%         xRelems = zeros(beammesh.Nelems,2);
%         yRelems = zeros(beammesh.Nelems,2);
%
%         xelems(:,1) = currpos(beammesh.connec(:,1),1);
%         xelems(:,2) = currpos(beammesh.connec(:,2),1);
%         yelems(:,1) = currpos(beammesh.connec(:,1),2);
%         yelems(:,2) = currpos(beammesh.connec(:,2),2);
%
%         xRelems(:,1) = beammesh.coords(beammesh.connec(:,1),1);
%         xRelems(:,2) = beammesh.coords(beammesh.connec(:,2),1);
%         yRelems(:,1) = beammesh.coords(beammesh.connec(:,1),2);
%         yRelems(:,2) = beammesh.coords(beammesh.connec(:,2),2);