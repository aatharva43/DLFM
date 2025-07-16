% Plots showing comparison of ground truth and reconstructed force patterns
clc
clear
close all

%% Simulation parameters
runparameters_Fig4A_inversevalidation; % set sigma_e = {0, 0.1, 0.2, 0.5} inside this file
% runparameters_inversevalidation_cell1;
% runparameters_inversevalidation_cell2;
% runparameters_inversevalidation_cell3;
% runparameters_inversevalidation_cell4;
% runparameters_inversevalidation_cell5;
% runparameters_inversevalidation_cell6;

%% Loop over all cases
for caseind = 1:Ncases

    fwdsolfoldername = sprintf(fwdsolfolderformat,cellIDmat(caseind));
    parentfolder     = sprintf(parentfolderformat,cellIDmat(caseind));
    frameNum         = framemat(caseind);

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
    saveFiguresFolder = sprintf('%s/Figures/frame%02d/wN%2.3E_wAC%2.3E/',parentfolder,frameNum,wt_nodalerror,wt_line);
    if ~isfolder(saveFiguresFolder)
        mkdir(saveFiguresFolder)
    end

    for regind = length(regcoeffmat)

        regcoeff   = regcoeffmat(regind);

        datafoldername      = sprintf('%s/frame%02d',parentfolder,frameNum);
        outputfoldername      = sprintf('%s/frame%02d/withAC_regcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,parentfolder,frameNum,regcoeff,E_actctr,dia_actctr,wt_nodalerror,wt_line,wt_el_pen);

        phaseimagefilename     = sprintf('%s/%02d_rot.jpg', datafoldername, frameNum);
        fakefluorimagefilename = sprintf('%s/fake_%02d_corrected_image_rot.jpg', datafoldername, frameNum);

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

        fwdsolfilename = sprintf('%s/ForwardProbSolution.mat',fwdsolfoldername);

        solfilename = sprintf('%s/InverseProbSolution_reg%2.4E.mat',outputfoldername,regcoeff);
        beammeshfilename = sprintf('%s/beammesh.mat',datafoldername);
        datax = load(sprintf('%s/n2_uxm.txt',datafoldername));
        datay = load(sprintf('%s/n2_uym.txt',datafoldername));

        if exist(solfilename)
            load(solfilename);
            load(beammeshfilename);
        else
            error('Could not find the solution file at\n%s',solfilename);
        end

        %% Reconstructed force and deformation data
        beammesh      = DOFhandlers(1);
        substratemesh = DOFhandlers(2);

        unodal = reshape(beammesh.u(:),3,beammesh.Nnodes)';
        currpos   = beammesh.coords + unodal(:,1:2);
        refpos    = beammesh.coords;

        % Get nodal forces
        nodalRHS = reshape(HiperProblem.Fnodal(1:DOFhandlers(1).Ndofs),3,beammesh.Nnodes)';
        nodalForces = nodalRHS(:,1:2);

        %% Ground truth data
        % Load the applied nodal forces in the forward problem
        fwdsolution      = load(sprintf('%s/ForwardProbSolution.mat',fwdsolfoldername));
        % To account for shift made to align the inverse problem mesh with
        % the generated nanonet image
        fwdsolution.DOFhandlers(1).coords = beammesh.coords;
        fwdsolution.DOFhandlers(2).coords = substratemesh.coords;

        fwdunodal    = reshape(fwdsolution.DOFhandlers(1).u(:),3,fwdsolution.DOFhandlers(1).Nnodes)';
        fwdcurrpos   = fwdsolution.DOFhandlers(1).coords + fwdunodal(:,1:2);

        nodalRHS_truth = reshape(fwdsolution.HiperProblem.Fnodal(1:fwdsolution.DOFhandlers(1).Ndofs),3,fwdsolution.DOFhandlers(1).Nnodes)';
        nodalForces_truth = nodalRHS_truth(:,1:2);

        %% Initialize plots
        fh_sigma   = figure();
        if plotnodalforce == 1
            fh_force   = figure();
        else
            fh_force   = -1;
        end
        if plotcontractility == 1
            fh_C = figure();
            polarscatter(pi,0,'k.');
        else
            fh_C = -1;
        end

        %% Plot force perunit length at element centers
        figure(fh_sigma); hold on;
        ha = findall(fh_sigma, 'type', 'axes');
        ha.Visible = 'On';

        % Set to 'reverse' when plotting image on top
        set(ha, 'YDir','reverse');

        % plot mesh
        figure(fh_sigma); hold on;

        % Fwd solution
        [xelems_intp_fwd,yelems_intp_fwd] = getHermiteInterpolatedCoords(fwdsolution.DOFhandlers(1),3);
        hdef_fwd = plot(xelems_intp_fwd',yelems_intp_fwd','-','LineWidth',1.0,'color',[0.75 0.75 0.75]);

        % Reg solution
        [xelems_intp,yelems_intp] = getHermiteInterpolatedCoords(beammesh,3);
        hdef = plot(xelems_intp',yelems_intp','-','LineWidth',2.0,'color','k');
        
        % Fwd solution
        [sigmamat_fwd,elemcencoords_fwd] = getLineTensionfromNodalForces...
            (fwdsolution.DOFhandlers(1),nodalForces_truth);

        sigma_mag_fwd = sqrt(sigmamat_fwd(:,1).^2 + sigmamat_fwd(:,2).^2);
        sigmax_fwd = max(sigma_mag_fwd);

        forcenodes_fwd = sigma_mag_fwd > 1e-2;
        figure(fh_sigma); hold on;
        quiver_tri(elemcencoords_fwd(forcenodes_fwd,1),elemcencoords_fwd(forcenodes_fwd,2),...
            sigmamat_fwd(forcenodes_fwd,1)*sigmascalefactor,sigmamat_fwd(forcenodes_fwd,2)*sigmascalefactor...
            ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,sigma_fwd_arrowcolor,...
            sigma_quifactor);

        % Reg solution
        [sigmamat,elemcencoords] = getLineTensionfromNodalForces(beammesh,nodalForces);

        sigma_mag = sqrt(sigmamat(:,1).^2 + sigmamat(:,2).^2);
        sigmax = max(sigma_mag);

        forcenodes = sigma_mag > 1e-2;
        figure(fh_sigma); hold on;
        quiver_tri(elemcencoords(forcenodes,1),elemcencoords(forcenodes,2),...
            sigmamat(forcenodes,1)*sigmascalefactor,sigmamat(forcenodes,2)*sigmascalefactor...
            ,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,sigma_arrowcolor,...
            sigma_quifactor);

        %% Plot nodal forces
        if plotnodalforce == 1
            figure(fh_force); hold on;
            ha = findall(fh_force, 'type', 'axes');            

            % Set to 'reverse' when plotting image on top
            set(ha, 'YDir','reverse');

            hdef_fwd = plot(xelems_intp_fwd',yelems_intp_fwd','-','LineWidth',1.0,'color',[0.75 0.75 0.75]);
            hdef     = plot(xelems_intp',yelems_intp','-','LineWidth',2.0,'color','k');


            forcenodes_fwd = vecnorm(nodalForces_truth,2,2) > 1e-2;
            quiver_tri(fwdcurrpos(forcenodes_fwd,1),fwdcurrpos(forcenodes_fwd,2),...
                nodalForces_truth(forcenodes_fwd,1)*forcescalefactor,nodalForces_truth(forcenodes_fwd,2)*forcescalefactor...
                ,force_arrowheadsize,22.5,force_arrowlinewidth,sigma_fwd_arrowcolor,...
                force_quifactor);

            forcenodes = vecnorm(nodalForces,2,2) > 1e-2;
            figure(fh_force); hold on;
            quiver_tri(currpos(forcenodes,1),currpos(forcenodes,2),...
                nodalForces(forcenodes,1)*forcescalefactor,nodalForces(forcenodes,2)*forcescalefactor...
                ,force_arrowheadsize,22.5,force_arrowlinewidth,force_arrowcolor,...
                force_quifactor);
        end

        %% Calculate contractility value
        if plotcontractility == 1
            % Get epicenter for the force distribution
            [coord_epicen,contractility,thCmax,Cmax,unitvecangle,contractility_theta]...
                = getContractilityEpicenter(nodalForces,currpos);
            fprintf(1,'Contractility = %3.2f nN, P = %2.3E\n',contractility,Cmax/contractility);

            % Get epicenter for the force distribution of fwd problem
            [~,contractility_fwd,~,Cmax_fwd,unitvecangle_fwd,contractility_theta_fwd]...
                = getContractilityEpicenter(nodalForces_truth,fwdcurrpos);

            figure(fh_C); hold on;
            polarplot(unitvecangle,contractility_theta/contractility_fwd,'linewidth',2.0,'Color',sigma_arrowcolor);
            polarplot(unitvecangle_fwd,contractility_theta_fwd/contractility_fwd,'linewidth',2.0,'Color','k','LineStyle','--');

            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            pax.RLim = [0,ceil(Cmax_fwd/contractility_fwd)];
            pax.RAxisLocation = 0;
            pax.GridLineStyle = '--';
            %     pax.ThetaTickLabel = '';
            pax.ThetaDir = 'clockwise'; % To match with the x-y orientation with flipped Yaxis
            title(sprintf('$C$=%2.1f nN',contractility_fwd))
        end

        %% Adjust figure size and add scale bars

        %% Line tension plot
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
        quiver_tri(xmax - 35,ymax-5,sigmascalefactor*sigma_scalebarval,0,sigma_arrowheadsize,22.5,sigma_arrowlinewidth,...
            sigma_arrowcolor,sigma_quifactor);

        if disp_showarrowtext == 1
            text(xmax - 35,(ymax-5)-6,sprintf('%d mN/m',sigma_scalebarval),'Color',sigma_arrowcolor,'FontSize',sigma_scalebar_textsize);
        end

        % Add lengthscale scale bar for tension figure
        LscaleArrow            = annotation('line') ;
        LscaleArrow.Parent     = ha;  % or any other existing axes or figure
        LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
        LscaleArrow.Color      = lenscalebar_color;
        LscaleArrow.LineWidth  = lenscalebar_linewidth;

        %% Nodal force plot
        if plotnodalforce == 1
            figure(fh_force);
            fh_force.Color = 'w';
            set(gcf,'InvertHardCopy','off');% Added to make the scale bar appear white
            ha = findall(fh_force, 'type', 'axes');
            box off
            axis equal
            axis off
            ha.Position = axisposition;
            ha.XLim     =  [xmin xmax];
            ha.YLim     =  [ymin ymax];

            fh_force.Units         = 'centimeters';
            fh_force.Position      = [1 1 figurewidth figurewidth*(ymax-ymin)/(xmax-xmin)];

            % Add tension and force scale bars
            % TO DO: Save a text file with all scale bars
            quiver_tri(xmax - 35,ymax-5,forcescalefactor*force_scalebarval,0,force_arrowheadsize,22.5,force_arrowlinewidth,...
                force_arrowcolor,force_quifactor);

            if disp_showarrowtext == 1
                text(xmax - 35,(ymax-5)-6,sprintf('%d mN/m',sigma_scalebarval),'Color',sigma_arrowcolor,'FontSize',sigma_scalebar_textsize);
            end

            % Add lengthscale scale bar for tension figure
            LscaleArrow            = annotation('line') ;
            LscaleArrow.Parent     = ha;  % or any other existing axes or figure
            LscaleArrow.Position   = [xmin + 10 ymax-5 lenscalebarval 0];
            LscaleArrow.Color      = lenscalebar_color;
            LscaleArrow.LineWidth  = lenscalebar_linewidth;
        end
        %% Contractility plot
        if plotcontractility == 1
            figure(fh_C);
            fh_C.Color         = 'w';
            fh_C.Units         = 'centimeters';
            fh_C.Position      = [16 2 12 12];

            hpolarax = gca;
            hpolarax.RTick     = [0,1,2];
            hpolarax.RLim      = [0,1.25];
            hpolarax.ThetaTick = [0.0:0.25:2]*pi;
            hpolarax.Position  = [0.2 0.15 0.6 0.65];
            hpolarax.FontSize  = 16;
        end

        %% Save figures
        fprintf(1,'Saving figures...\n');
        imresolution = '-r300'; % resolution at which to save figures

        print(fh_sigma,sprintf('%s/frame%d_Reg%2.3E.png',saveFiguresFolder,frameNum,regcoeff),'-dpng',imresolution);
        print(fh_sigma,sprintf('%s/frame%d_Reg%2.3E.eps',saveFiguresFolder,frameNum,regcoeff),'-depsc');
        
        if plotnodalforce ==1
            print(fh_force,sprintf('%s/frame%d_Reg%2.3E_Forces.png',saveFiguresFolder,frameNum,regcoeff),'-dpng',imresolution);
            print(fh_force,sprintf('%s/frame%d_Reg%2.3E_Forces.eps',saveFiguresFolder,frameNum,regcoeff),'-depsc');
        end

        if plotcontractility == 1
            print(fh_C,sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.png',saveFiguresFolder,regind,frameNum,regcoeff),'-dpng',imresolution);
            print(fh_C,sprintf('%s/%02d_Contractility_frame%d_Reg%2.3E.eps',saveFiguresFolder,regind,frameNum,regcoeff),'-depsc');
        end

        scalefile = sprintf('%s/scale.txt', saveFiguresFolder);
        fileID = fopen(scalefile,'w');
        fprintf(fileID,'Length scale: %d um \n',lenscalebarval);
        fprintf(fileID,'Tension scale: %d mN/m \n',sigma_max);
        if plotnodalforce == 1
            fprintf(fileID,'Force scale: %d mN/m \n',force_max);
        end
        fclose(fileID);


        %% close figures
        if closefigures == 1
            close(fh_sigma);
            if plotcontractility == 1
                close(fh_C);
            end
            if plotnodalforce == 1
                close(fh_force);
            end
        end

    end
end

function [sigmamat,elemcencoords] = getLineTensionfromNodalForces(beammesh,nodalForces)
% Get line tension from nodal forces
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

unodal = reshape(beammesh.u(:),3,beammesh.Nnodes)';
currpos   = beammesh.coords + unodal(:,1:2);
elemcencoords = 0.5*((currpos(beammesh.connec(:,1),1:2)) + (currpos(beammesh.connec(:,2),1:2)));
end