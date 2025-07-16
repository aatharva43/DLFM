% Compute and plot the error between the ground truth solution of the
% forward problem and the inverse solution
clc
clear
close all

%% Simulation parameters
% runparameters_Fig2C;
runparameters_Fig4A_inversevalidation;

%% Loop over all cases
rpatch   =  20;
plotcontractility = 0;
addlegend = 0;

if plotcontractility == 1
    fh_C = figure(1);
    polarscatter(pi,0,'k.');
end
for caseind = 1:Ncases

    fwdsolfoldername = sprintf(fwdsolfolderformat,cellIDmat(caseind));
    parentfolder = sprintf(parentfolderformat,cellIDmat(caseind));
    frameNum     = framemat(caseind);

    % Energy term weights
    wt_nodalerror    = wtNmat(caseind);
    wt_line          = wtLmat(caseind);
    wt_el_pen        = wtEpmat(caseind);
    
    E_actctr         = E_actctrmat(caseind);
    dia_actctr       = dia_actctrmat(caseind);

    % Creating Folder for saving tension figures
    saveFiguresFolder = sprintf('%s/Figures/frame%02d/wN%2.3E_wAC%2.3E/',parentfolder,frameNum,wt_nodalerror,wt_line);
    if ~isfolder(saveFiguresFolder)
        mkdir(saveFiguresFolder)
    end

    % Read measured data (dof wise)
    mDOFs   = [];
    % x-displacements
    datafoldername      = sprintf('%s/frame%02d',parentfolder,frameNum); % input data
    datax = load(sprintf('%s/n2_uxm.txt',datafoldername));
    mDOFs = [mDOFs;3*(datax(:,1)-1)+1];
    % y-displacements
    datay = load(sprintf('%s/n2_uym.txt',datafoldername));
    mDOFs = [mDOFs;3*(datay(:,1)-1)+2];
    Nm = length(mDOFs); % number of measurements (x and y)

    errormat      = [];
    ferrormat      = [];
    Jtotal_mat    = zeros(size(regcoeffmat));

    linecolors = cbrewer('Dark2', length(regcoeffmat),'linear');

    for regind = 1:length(regcoeffmat)
        regcoeff = regcoeffmat(regind);
        outputfoldername      = sprintf('%s/frame%02d/withAC_regcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,parentfolder,frameNum,regcoeffmat(regind),E_actctr,dia_actctr,wt_nodalerror,wt_line,wt_el_pen); % output and run data

        errordata = load(sprintf('%s/errordata.txt',outputfoldername));
        errormat = [errormat;errordata(end,:)];

        regsolfilename      = sprintf('%s/InverseProbSolution_reg%2.4E.mat',outputfoldername,regcoeffmat(regind));
        load(regsolfilename,'JData');

        forceerrordata      = getForceErrors(fwdsolfoldername,regsolfilename,rpatch);
        ferrormat = [ferrormat;forceerrordata];

        Jtotal_mat(regind) = wt_nodalerror * JData.J_m + JData.J_reg + wt_line * JData.J_m_line + wt_el_pen * JData.J_reg_el;

        %% Plot contractility values
        % reg solution data
        solfilename = sprintf('%s/InverseProbSolution_reg%2.4E.mat',outputfoldername,regcoeff);
        load(solfilename);
        % Reconstructed force and deformation data
        beammesh      = DOFhandlers(1);

        % Ground truth data
        fwdsolfilename = sprintf('%s/ForwardProbSolution.mat',fwdsolfoldername);
        
        % Load the applied nodal forces in the forward problem
        fwdsolution      = load(sprintf('%s/ForwardProbSolution.mat',fwdsolfoldername));
        % To account for shift made to align the inverse problem mesh with
        % the generated nanonet image
        fwdsolution.DOFhandlers(1).coords = beammesh.coords;

        fwdunodal    = reshape(fwdsolution.DOFhandlers(1).u(:),3,fwdsolution.DOFhandlers(1).Nnodes)';
        fwdcurrpos   = fwdsolution.DOFhandlers(1).coords + fwdunodal(:,1:2);

        nodalRHS_truth = reshape(fwdsolution.HiperProblem.Fnodal(1:fwdsolution.DOFhandlers(1).Ndofs),3,fwdsolution.DOFhandlers(1).Nnodes)';
        nodalForces_truth = nodalRHS_truth(:,1:2);

        unodal = reshape(beammesh.u(:),3,beammesh.Nnodes)';
        currpos   = beammesh.coords + unodal(:,1:2);

        % Get nodal forces
        nodalRHS = reshape(HiperProblem.Fnodal(1:DOFhandlers(1).Ndofs),3,beammesh.Nnodes)';
        nodalForces = nodalRHS(:,1:2);

        % Get epicenter for the force distribution
        [coord_epicen,contractility,thCmax,Cmax,unitvecangle,contractility_theta]...
            = getContractilityEpicenter(nodalForces,currpos);

        % Get epicenter for the force distribution of fwd problem
        [~,contractility_fwd,~,Cmax_fwd,unitvecangle_fwd,contractility_theta_fwd]...
            = getContractilityEpicenter(nodalForces_truth,fwdcurrpos);

        % plot and compare
        if plotcontractility == 1
            figure(fh_C); hold on;

            if regind == 1
                polarplot(unitvecangle_fwd,contractility_theta_fwd/contractility_fwd,...
                    'linewidth',2.0,'Color','k','LineStyle','--',...
                    'DisplayName', 'target solution');
            end
            linelabel = sprintf('%1.2E',regcoeff);
            polarplot(unitvecangle,contractility_theta/contractility_fwd,...
                'linewidth',2.0,'Color',linecolors(regind,:),...
                'DisplayName',linelabel);

            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            pax.RLim = [0,ceil(Cmax_fwd/contractility_fwd)];
            pax.RAxisLocation = 0;
            pax.GridLineStyle = '--';
            %     pax.ThetaTickLabel = '';
            pax.ThetaDir = 'clockwise'; % To match with the x-y orientation with flipped Yaxis
            title(sprintf('$C$=%2.1f nN',contractility_fwd))
        end
    end

    lambdamat = errormat(:,1);
    [lambdamat, lamorder] = sort(lambdamat);

    resnorm    = errormat(lamorder,2)/sqrt(Nm);
    solnorm    = errormat(lamorder,3);
    Jtotal_mat = Jtotal_mat(lamorder);

    ferrorpointwise     = ferrormat(lamorder,1)*100; % [%]
    contractility       = ferrormat(lamorder,2);
    contractility_truth = ferrormat(lamorder,3);
    ferrorpatch         = ferrormat(lamorder,4)*100; % [%]

    contractility_error = 100* abs(contractility - contractility_truth)./abs(contractility_truth);

    %% Plot force error data
    fhFpatch = figure(2);
    figure(fhFpatch); hold on;
    hplt = plot(lambdamat,ferrorpatch);
    hplt.DisplayName = sprintf('%2.2f',sigma_e);
    hplt_markers = plot(lambdamat,ferrorpatch);
    hplt_markers.Annotation.LegendInformation.IconDisplayStyle = 'off';

    set(gca, 'XScale', 'log');
    xlabel('$\beta$ [nN/$\mu$m]','interpreter','latex');
    ylabel('Force recovery error [\%]','interpreter','latex');

    fhFpatch.Color = 'w';
    fhFpatch.Units = 'centimeters';
    fhFpatch.Position = [3 3 12 10];
    fhFpatch.Renderer = 'Painters';
    box off;
    grid on;

    axh = gca;
    set(axh,'Pos',[0.2 0.2 0.7 0.7]);
    set(axh,'FontSize',14);
    axh.XLabel.FontSize = 18;
    axh.YLabel.FontSize = 18;
    axh.XScale = 'log';
    axh.XTick  = 10.^[-5:1:5];
    axh.YTick  = 0:10:100;
    axh.XLim   = [10^-3, 10];
    axh.YLim   = [0, 100];
    axh.XMinorGrid = 'off';
    axh.YMinorGrid = 'off';

    hplt.LineWidth       = 2.0;
    hplt.Color           = colorval;

    hplt_markers.Color           = 'none';
    hplt_markers.LineWidth       = 1.0;
    hplt_markers.Marker          = 'square';
    hplt_markers.MarkerFaceColor = colorval;
    hplt_markers.MarkerEdgeColor = 'k';
    hplt_markers.MarkerSize      = 6;

    ylim([0,100]);

    if addlegend == 1
        hleg = legend;
        hleg.Box = 'off';
        hleg.Location = 'northwest';
        hleg.FontSize = 14;
    end

    Cerrfilename = sprintf('%s/Fpatcherr_Frame%02d',saveFiguresFolder,frameNum);
    print(fhFpatch, [Cerrfilename,'.png'],'-dpng','-r300');
    print(fhFpatch, [Cerrfilename,'.eps'],'-depsc');

    fp = fopen(sprintf('%s/patchsize_Frame%02d.txt',saveFiguresFolder,frameNum),'w+');
    fprintf(fp,'rp = %2.3f\n', rpatch);
    fclose(fp);


    %% Plot contractility error data
    fhCerr = figure(3); 
    figure(fhCerr); hold on;
    hplt = plot(lambdamat,contractility_error);
    hplt.DisplayName = sprintf('%2.2f',sigma_e);
    hplt_markers = plot(lambdamat,contractility_error);
    hplt_markers.Annotation.LegendInformation.IconDisplayStyle = 'off';

    set(gca, 'XScale', 'log');
    xlabel('$\beta$ [nN/$\mu$m]','interpreter','latex');
    ylabel('Contractility recovery error [\%]','interpreter','latex');

    fhCerr.Color = 'w';
    fhCerr.Units = 'centimeters';
    fhCerr.Position = [3 3 12 10];
    fhCerr.Renderer = 'Painters';
    box off;
    grid on;

    axh = gca;
    set(axh,'Pos',[0.2 0.2 0.7 0.7]);
    set(axh,'FontSize',14);
    axh.XLabel.FontSize = 18;
    axh.YLabel.FontSize = 18;
    axh.XScale = 'log';
    axh.XTick  = 10.^[-5:1:5];
    axh.YTick  = 0:10:100;
    axh.XLim   = [10^-3, 10];
    axh.YLim   = [0, 100];
    axh.XMinorGrid = 'off';
    axh.YMinorGrid = 'off';

    hplt.LineWidth       = 2.0;
    hplt.Color           = colorval;

    hplt_markers.Color           = 'none';
    hplt_markers.LineWidth       = 1.0;
    hplt_markers.Marker          = 'square';
    hplt_markers.MarkerFaceColor = colorval;
    hplt_markers.MarkerEdgeColor = 'k';
    hplt_markers.MarkerSize      = 6;

    ylim([0,100]);
    
    if addlegend == 1
        hleg = legend;
        hleg.Box = 'off';
        hleg.Location = 'northwest';
        hleg.FontSize = 14;
    end

    Cerrfilename = sprintf('%s/Cerr_Frame%02d',saveFiguresFolder,frameNum);
    print(fhCerr, [Cerrfilename,'.png'],'-dpng','-r300');
    print(fhCerr, [Cerrfilename,'.eps'],'-depsc');

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

end


function forceerrordata = getForceErrors(fwdsolfoldername,regsolfilename,rpatch)
% Calculate error between force patterns

% Load the reconstructed force pattern
load(regsolfilename);
nodalRHS = reshape(HiperProblem.Fnodal(1:DOFhandlers(1).Ndofs),3,DOFhandlers(1).Nnodes)';
nodalForces = nodalRHS(:,1:2);

% Current and Ref positions
unodal    = reshape(DOFhandlers(1).u(:),3,DOFhandlers(1).Nnodes)';
currpos   = DOFhandlers(1).coords + unodal(:,1:2);
refpos    = DOFhandlers(1).coords;

% Get epicenter for the force distribution
[~,contractility] = getContractilityEpicenter(nodalForces,currpos);

% Load the applied nodal forces in the forward problem
bcNeu               = load(sprintf('%s/dofExternalForces.txt',fwdsolfoldername));
RHS_truth = zeros(size(HiperProblem.Fnodal));
RHS_truth(bcNeu(:,1)) = bcNeu(:,2);

nodalRHS_truth = reshape(RHS_truth(1:DOFhandlers(1).Ndofs),3,DOFhandlers(1).Nnodes)';
nodalForces_truth = nodalRHS_truth(:,1:2);

fwdsolution      = load(sprintf('%s/ForwardProbSolution.mat',fwdsolfoldername));
fwdunodal    = reshape(fwdsolution.DOFhandlers(1).u(:),3,fwdsolution.DOFhandlers(1).Nnodes)';
fwdcurrpos   = fwdsolution.DOFhandlers(1).coords + fwdunodal(:,1:2);

[~,contractility_truth] = getContractilityEpicenter(nodalForces_truth,fwdcurrpos);

% Pointwise error
Ferror_pointwise = norm(vecnorm(nodalForces_truth - nodalForces,2,2))...
    /norm(vecnorm(nodalForces_truth,2,2));

% Sum of forces in a small radius around a point
magnodalForces_truth = vecnorm(nodalForces_truth,2,2);
forcenodes = find( magnodalForces_truth > 0.1 * max(magnodalForces_truth) );
%rpatch     = 20;

% visualize patch
% fh = figure(); hold on;
% plot(currpos(:,1),currpos(:,2),'k.');
% 
% pch = plot(currpos(1,1),currpos(1,2),'ks','MarkerSize',12,'MarkerFaceColor','b');
% pathnodesh = plot(currpos(1,1),currpos(1,2),'ko','MarkerSize',6,'MarkerFaceColor','r');

for i = 1:length(forcenodes)
    nodeID = forcenodes(i);
    xnode = refpos(nodeID,1);
    ynode = refpos(nodeID,2);
    rlocal = sqrt( (refpos(:,1) - xnode).^2 + ( refpos(:,2) -ynode).^2 );
    patchnodes = find(rlocal <= rpatch);

    % error within patch
    Fpatch_truth = sum(nodalForces_truth(patchnodes,:),1);
    Fpatch       = sum(nodalForces(patchnodes,:),1);

    Ferror_patch(i) = norm(vecnorm(Fpatch_truth - Fpatch,2,2))...
        /norm(vecnorm(Fpatch_truth,2,2));

    %     pch.XData = xnode;
    %     pch.YData = ynode;
    %
    %     pathnodesh.XData = currpos(patchnodes,1);
    %     pathnodesh.YData = currpos(patchnodes,2);
end

mFerror_patch = mean(Ferror_patch);

forceerrordata = [Ferror_pointwise, contractility, contractility_truth, mFerror_patch];
end
