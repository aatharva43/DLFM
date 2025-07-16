clc
clear
close all

addpath(genpath('./src'));

%% Simulation parameters
runparameters_Fig2C;
% runparameters_Fig4A_inversevalidation;

%% Loop over all cases
if exist("colorval","var")
    linecolor = colorval;
else
    linecolor = 'r';
end
addpointlabels = 0;
for caseind = 1:Ncases

    parentfolder = sprintf(parentfolderformat,cellIDmat(caseind));
    frameNum     = framemat(caseind);

    % Energy term weights
    wt_nodalerror    = wtNmat(caseind);
    wt_line          = wtLmat(caseind);
    wt_el_pen        = wtEpmat(caseind);
    
    E_actctr         = E_actctrmat(caseind);
    dia_actctr       = dia_actctrmat(caseind);

    
    % Read measured data (dof wise)
    mDOFs   = [];
    % x-displacements
    datafoldername      = sprintf('%s/frame%02d',parentfolder,frameNum); % input data
    datax = load(sprintf('%s/n2_uxm.txt',datafoldername));
    mDOFs = [mDOFs;3*(datax(:,1)-1)+1];
    % y-displacements
    datay = load(sprintf('%s/n2_uym.txt',datafoldername));
    mDOFs = [mDOFs;3*(datay(:,1)-1)+2];
    Nm = length(mDOFs);

    errormat      = [];
    Jtotal_mat    = zeros(size(regcoeffmat));
    for regcoeffind = 1:length(regcoeffmat)
        outputfoldername      = sprintf('%s/frame%02d/withAC_regcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,parentfolder,frameNum,regcoeffmat(regcoeffind),E_actctr,dia_actctr,wt_nodalerror,wt_line,wt_el_pen); % output and run data

        errordata = load(sprintf('%s/errordata.txt',outputfoldername));
        errormat = [errormat;errordata(end,:)];

        solfilename      = sprintf('%s/InverseProbSolution_reg%2.4E.mat',outputfoldername,regcoeffmat(regcoeffind));
        load(solfilename);

        Jtotal_mat(regcoeffind) = wt_nodalerror * JData.J_m + JData.J_reg + wt_line * JData.J_m_line + wt_el_pen * JData.J_reg_el;
    end

    lambdamat = errormat(:,1);
    [lambdamat, lamorder] = sort(lambdamat);

    resnorm   = errormat(lamorder,2)/sqrt(Nm);
    solnorm   = errormat(lamorder,3);
    Jtotal_mat = Jtotal_mat(lamorder);

    %% plot and save L-curve with only measured data
    fhLcurve = figure(2); hold on;
    hplt = plot(resnorm,solnorm);
    %hplt.DisplayName = sprintf('%2.2f',standdev);
    hplt_markers = plot(resnorm,solnorm);
    hplt_markers.Annotation.LegendInformation.IconDisplayStyle = 'off';

    xlabel('$||\mathbf{u} - \mathbf{u}_m||/\sqrt{N_m}$ $\mu$m','Interpreter','latex');
    ylabel('$||\mathbf{f}||$','Interpreter','latex');

    fhLcurve.Color = 'w';
    fhLcurve.Units = 'centimeters';
    fhLcurve.Position = [3 3 12 10];
    fhLcurve.Renderer = 'Painters';
    box off;
    grid on;

    axh = gca;
    set(axh,'Pos',[0.2 0.2 0.7 0.7]);
    set(axh,'FontSize',14);
    axh.XLabel.FontSize = 18;
    axh.YLabel.FontSize = 18;
    axh.XScale = 'log';
    axh.YScale = 'log';
    %axh.YTick  = 10.^[-5:1:5];
    %axh.XTick  = 10.^[-5:1:5];
    axh.XTick  = [0.1, 0.2, 0.5, 1];
    axh.XLim   = [0.01, 1.2];
    axh.YLim   = [0.1, 10^4];
    axh.XMinorGrid = 'off';
    axh.YMinorGrid = 'off';

    hplt.LineWidth       = 2.0;
    hplt.Color           = linecolor;

    hplt_markers.Color           = 'none';
    hplt_markers.LineWidth       = 1.0;
    hplt_markers.Marker          = 'square';
    hplt_markers.MarkerFaceColor = linecolor;
    hplt_markers.MarkerEdgeColor = 'k';
    hplt_markers.MarkerSize      = 6;

    if addpointlabels == 1
        labelpoints(resnorm,solnorm,lambdamat,...
            'NE',0, 1,...
            'FontSize', 12);
    end

    Lcurvefilename = sprintf('%s/LCurve_Frame%02d',parentfolder,frameNum);
    print(fhLcurve, [Lcurvefilename,'.png'],'-dpng','-r300');
    print(fhLcurve, [Lcurvefilename,'.eps'],'-depsc');

end

