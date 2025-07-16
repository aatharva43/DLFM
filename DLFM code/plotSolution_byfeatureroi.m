% Analyzing forces within a certain feature
% Plot data of all features with different colored arrows
% Plotting both together for a movie or saving them seperately
clc;
clear;
close all;
addpath(genpath('./src'));

set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Simulation parameters
% runparameters_Experimental_Demo;
% runparameters_actin_inhibition;
% runparameters_ROCK_inhibition;
% runparameters_Rectangular_Grid;
% runparameters_Hybrid_Grid;
% runparameters_Cell_Division;
% runparameters_Cell_Spreading;
% runparameters_inversevalidation_cell1;
% runparameters_inversevalidation_cell2;
% runparameters_inversevalidation_cell3;

% runparameters_Fig3C_Inhibition;
% runparameters_Fig3D_Y27632;
% runparameters_Fig3D_Y27632washed;
% runparameters_Fig3E;
% runparameters_force_area_plot;
% runparameters_CNFM_Schematic;
% runparameters_Fig3C_Control_series;
% runparameters_Fig3E_series;
runparameters_Fig3D_Control_series;
% runparameters_Cell_Spreading_series;
% runparameters_Cell_Division_series;

% ID of the feature 
featureroiIDmat = 1;%1:3
featurecolormat = zeros(length(featureroiIDmat),3);

featurecolormat(1,:) = [0.5 1.0 0.0];
featurecolormat(2,:) = [255 69 0]/255;
featurecolormat(3,:) = [0 1 1];

save_featurewise = 0;

plotcombined = 1; % Plot combined for movie
plotcontractility = 1;
showframetime = 1;
Cmax_rlim = 1500;%3000;

% movietext = '';
movietext = 'Control';
% movietext = 'Latrunculin B';
% movietext = 'Y27632';
% movietext = 'After wash';

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
    if save_featurewise == 0
        featurefolder = sprintf('%s/feature_combined',parentfolder);
    else
        featurefolder = sprintf('%s/feature_%02d',parentfolder,featureroiIDmat);
    end

    if ~isfolder(featurefolder)
        mkdir(featurefolder);
    end

    CombinedMovie      = sprintf('%s/CombinedMovie_Fig',featurefolder);
    if ~isfolder(CombinedMovie) 
        mkdir(CombinedMovie) 
    end
    
    TensionFigures      = sprintf('%s/Tension_Fig',featurefolder);
    ContracFigures      = sprintf('%s/Contrac_Fig',featurefolder);
    if ~isfolder(TensionFigures) 
        mkdir(TensionFigures) 
    end

    if ~isfolder(ContracFigures) 
        mkdir(ContracFigures) 
    end
    
    for regind = 1:length(regcoeffmat)
        
        regcoeff   = regcoeffmat(regind);
        
        datafoldername      = sprintf('%s/frame%02d',parentfolder,frameNum);
        outputfoldername      = sprintf('%s/frame%02d/withAC_regcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,parentfolder,frameNum,regcoeff,E_actctr,dia_actctr,wt_nodalerror,wt_line,wt_el_pen);
        
        cellmaskfile = sprintf('%s/%02d_cellROI.mat',datafoldername,frameNum);
        cellmask = load(cellmaskfile);
        maskcoords = cellmask.hroipoly.Position;

        imageresolutionfile  = sprintf('%s/frame%02d/pix2um.txt',parentfolder,frameNum);
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
            beammeshfilename = sprintf('%s/beammesh.mat',datafoldername);
            rundatafilename = sprintf('%s/rundata.mat',outputfoldername);
            datax = load(sprintf('%s/n2_uxm.txt',datafoldername));
            datay = load(sprintf('%s/n2_uym.txt',datafoldername));
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

        unodal    = reshape(beammesh.u(:),3,beammesh.Nnodes)';
        currpos   = beammesh.coords + dispampfactor*unodal(:,1:2);
        refpos    = beammesh.coords;

        %% Initialize plots
        if plotsigma == 1
            fh_sigma   = figure();
        else
            fh_sigma = -1;
        end
        if plotcontractility == 1
            if plotcombined == 1
                fh_C = figure(fh_sigma);
            else
                 fh_C = figure();
            end
            polarscatter(pi,0,'k.');
        else
            fh_C = -1;
        end
        
        %% Background image (phase) and mesh plot
        if plotsigma == 1
            figure(fh_sigma); hold on;
            if plotcombined == 1
                subplot(1,2,1);
            end

            cellimagefile = phaseimagefile;
            if isfile(cellimagefile)
                cellimg = imread(cellimagefile);

                % Reduce brighness for arrow visibility, for RGB images
                if length(size(cellimg)) == 3
                    cellimghsv = rgb2hsv(cellimg);
                    cellimghsv(:,:,3) = cellimghsv(:,:,3)/(mean(mean(cellimghsv(:,:,3))))*0.55;
                    cellimgconv = hsv2rgb(cellimghsv);
                else
                    cellimgconv = cellimg;
                end
                
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
                subplot(1,2,1);
                hdef = plot(xelems_intp',yelems_intp','-','LineWidth',1.0,'color','k');
            end
        end
        
       
        %% Get nodal forces
        nodalRHS = reshape(HiperProblem.Fnodal(1:DOFhandlers(1).Ndofs),3,beammesh.Nnodes)';
        

        for featureroiID = featureroiIDmat
            sigma_arrowcolor_feature = featurecolormat(find(featureroiIDmat==featureroiID),:);
            % Get feature ROI
            featureroifile = sprintf('%s/%02d_feature%02d_ROI.mat',datafoldername,frameNum,featureroiID);
            featureroi     = load(featureroifile);
            featurecoords  = featureroi.hroipoly.Position;

            featureNodes = find(inpolygon(currpos(:,1),currpos(:,2),...
                featurecoords(:,1),featurecoords(:,2)));

            nodalForces = zeros(size(nodalRHS,1),2);
            % only use feature nodes
            nodalForces(featureNodes,1:2) = nodalRHS(featureNodes,1:2);

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
                if plotcombined == 1
                    subplot(1,2,1);
                end
                quiver_tri(elemcencoords(forcenodes,1),elemcencoords(forcenodes,2),...
                    sigmamat(forcenodes,1)*sigmascalefactor,sigmamat(forcenodes,2)*sigmascalefactor...
                    ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,sigma_arrowcolor_feature,...
                    sigma_quifactor);

            end

            %% Calculate contractility value
            % Get epicenter for the force distribution
            [coord_epicen,contractility,thCmax,Cmax,unitvecangle,contractility_theta]...
                = getContractilityEpicenter(nodalForces,currpos,plotcontractility,fh_C);
            fprintf(1,'Contractility = %3.2f nN, P = %2.3E\n',contractility,Cmax/contractility);

            if plotcontractility == 1
                figure(fh_C); hold on;
                if plotcombined == 1
                    subplot(1,2,2);
                end
                polarplot(unitvecangle,contractility_theta,'-','linewidth',4.0,'Color',sigma_arrowcolor_feature);
            end

            % Plot contractile epicenter
            % if plotsigma == 1
            %     figure(fh_sigma);
            %     subplot(1,2,1);hold on;
            %     plot(coord_epicen(1),coord_epicen(2),'ks','MarkerSize',8,'MarkerFaceColor','r');
            % end
        end
        
        %% Adjust figure size and add scale bars

        %% sigma with phase back image
        if plotsigma == 1
        figure(fh_sigma);
        if plotcombined == 1
            subplot(1,2,1);
        end

        fh_sigma.Units         = 'centimeters';
        if plotcombined == 1
            fh_sigma.Position      = [1 1 figurewidth+figurewidth/2  figurewidth*(ymax-ymin)/(xmax-xmin)];
        else
            fh_sigma.Position      = [1 1 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];
        end
        fh_sigma.Color = 'w';
        set(gcf,'InvertHardCopy','off');% Added to make the scale bar appear white
        ha = findall(fh_sigma, 'type', 'axes');
        box off
        axis equal
        axis off
        if plotcombined == 1
            ha.Position =  [0 0 0.5 1.0];
        else
            ha.Position =  [0 0 1.0 1.0];
        end
        ha.XLim     =  [xmin xmax];
        ha.YLim     =  [ymin ymax];
        
        % Add tension and force scale bars
        % TO DO: Save a text file with all scale bars
        if plotlinetension == 1
            quiver_tri(xmax - 35,ymax-5,sigmascalefactor*sigma_scalebarval,0,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,...
                sigma_arrowcolor,sigma_quifactor);
            
            if disp_showarrowtext == 1
                text(xmax - 35,(ymax-5)-6,sprintf('%d mN/m',sigma_scalebarval),'Color',sigma_arrowcolor,'FontSize',sigma_scalebar_textsize);
            end
        end
        
        % Add lengthscale scale bar for tension figure
        LscaleArrow            = annotation('line') ;
        LscaleArrow.Parent     = ha;  % or any other existing axes or figure
        LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
        LscaleArrow.Color      = lenscalebar_color;
        LscaleArrow.LineWidth  =  lenscalebar_linewidth;
        end

        %% Add timestap
        if plotcombined == 1
            if showframetime == 1
                str = sprintf('%d min',round(frametime(caseind)));
                htimest = annotation('textbox',[0.025,0.9,0.1,0.1],'String',str);
                htimest.FontSize = 16;
                htimest.LineStyle = 'none';
            end

            if ~isempty(movietext)
                str = movietext;
                htimest = annotation('textbox',[0.5 - 0.05 0.0 0.1 0.1],'String',str);
                htimest.FontSize = 20;
                htimest.LineStyle = 'none';
                htimest.HorizontalAlignment = 'center';
                htimest.VerticalAlignment = 'middle';
            end
        end

        %% Contractility plot
        if plotcontractility == 1
            figure(fh_C);
            fh_C.Color         = 'w';
            fh_C.Units         = 'centimeters';
            
            if plotcombined == 1
                pax=subplot(1,2,2);
                pax.OuterPosition  = [0.5,0,0.5,1.0];
                pax.InnerPosition  = [0.55, 0.1, 0.4, 0.75];
            else
                pax = gca;
                fh_C.Position      = [16 2 12 12];
            end
            
            pax.ThetaDir       = 'clockwise'; % To match with the x-y orientation with flipped Yaxis
            pax.ThetaAxisUnits = 'radians';
            pax.ThetaTick      = [0, pi/2, pi, 3*pi/2];
            pax.ThetaTickLabel = {};

            
            pax.RLim = [0, Cmax_rlim];
            pax.RTick = linspace(0,Cmax_rlim,3);
            pax.RAxisLocation = 3*pi/2;
            
            pax.GridLineStyle = '-';
            
            pax.FontSize = 18;
            pax.Layer = 'top';

            %             ht = title(sprintf('$C_\\theta$ (nN), $C$ = %3.1f (nN)',contractility));
            if plotcombined == 1
                ht = title(sprintf('$C_\\theta$ (nN)'));
                ht.FontSize = 20;
                ht.Units = "normalized";
                ht.Position = [0.5 1.1 0];
            end

        end
        
        
        %% Save images
        fprintf(1,'Saving figures...\n');
        imresolution = '-r300'; % resolution at which to save figures
        saveFiguresFolder = sprintf('%s/Figures/frame%d',featurefolder,frameNum);
        if exist(saveFiguresFolder,'dir')~=7
            mkdir(saveFiguresFolder) ;
        end
        
        if plotcombined == 1
            %print(fh_sigma,sprintf('%s/frame%04d.png',CombinedMovie,caseind),'-dpng',imresolution);
            print(fh_sigma,sprintf('%s/frame%04d',CombinedMovie,caseind),'-dpng',imresolution);
        else
            if plotsigma == 1
            print(fh_sigma,sprintf('%s/frame%04d.png',TensionFigures,frameNum),'-dpng',imresolution);
            end

            if plotcontractility == 1
                print(fh_C,sprintf('%s/PolarContractility_frame%04d.png',ContracFigures,frameNum),'-dpng',imresolution);
                print(fh_C,sprintf('%s/PolarContractility_frame%04d.eps',ContracFigures,frameNum),'-depsc');

                % save textfile to plot data seperately
                fp_cdata = fopen(sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.txt',saveFiguresFolder,regind,frameNum,regcoeff),'w+');
                fprintf(fp_cdata,'Theta  C_theta  C\n');
                fprintf(fp_cdata,'%2.3E %2.3E %2.3E\n',[unitvecangle(:),contractility_theta(:),contractility*ones(length(unitvecangle),1)]');
                fclose(fp_cdata);
            end
        end
        
        if plotcombined == 1
        scalefile = sprintf('%s/scale.txt', CombinedMovie);
        else
            scalefile = sprintf('%s/scale.txt', TensionFigures);
        end
       fileID = fopen(scalefile,'w');
       fprintf(fileID,'Length scale: %d um \n',lenscalebarval);
       fprintf(fileID,'Tension scale: %d mN/m \n',sigma_max);
       fclose(fileID);
        
            
        %% close figures
        if closefigures == 1
            if plotsigma == 1
                close(fh_sigma);
            end
            if plotcombined == 0
                if plotcontractility == 1
                    close(fh_C);
                end
            end
        end
        
    end
end


%% Movie command
% ffmpeg -f image2 -r 4 -i ./frame%04d.png -crf 10 ./celldiv.wmv
% ffmpeg -r 2 -f image2 -s 1920x1080 -i frame%04d.png -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -crf 25  -pix_fmt yuv420p M6.mov