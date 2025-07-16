%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crosshatch Nanonet force microscopy
% To analyze frames in series using the same mesh
% Running instructions and tips:
%
% (1) Use preprocessing_fluoroimages.m to generate the filtered images for
% active contour matching.
% (2) Use preprocessing_series.m to generte the input data for the inverse
% analysis.
% (3) If the computational grid does not snap into
% the fluoroscent fiber shapes, try increasing the wtLmat value. 
% (4) The first step would be to generate the L-curve for a range of
% regularization parmeters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
addpath(genpath('./src'));

%% Simulation parameters
% runparameters_CNFM_Schematic;
% runparameters_Cell_Spreading_series;
% runparameters_Fig3C_Control_series;
% runparameters_Fig3E_series;
% runparameters_Fig3D_Control_series;
runparameters_Cell_Division_series;
totalruntime = tic;

%% Perform analysis
 for caseind = 1:Ncases
    frameNum         = framemat(caseind);
    parentfolder     = sprintf(parentfolderformat,cellIDmat(caseind));
    framefolder      = sprintf(['%s/',framefolderformat],parentfolder,frameNum);

    jobparams                   = []; % to be able to work with parfor

    % Plot and save intermediate converged step solutions
    jobparams.plotupdates       = plotupdates;

    % To copy last regcoefficient Gammamat files for restart
    useprevGamma_set            = 1;
    % To restart the minimization from existing Gammamat and u
    % restart_u.txt and restart_Gamma.txt
    jobparams.restart_gamma     = 1;

    % perturb initial guess for restarts
    jobparams.randomize_Gammainit = 0;
    jobparams.std_rand_Gammainit  = 0;
   
    imageresolutionfile         = sprintf('%s/pix2um.txt',framefolder);
    jobparams.pix2um            = load(imageresolutionfile);

    % internal use of restart from previous known solutions for speed-up
    jobparams.restart_solver    = 1;
    
    jobparams.E                  = E;             % [kPa]
    jobparams.dia                = dia;           % [um]
    jobparams.pretensionF        = pretensionF;   % [nN], 200 nN
    
    % Energy term weights
    jobparams.wt_nodalerror    = wtNmat(caseind);
    jobparams.wt_line          = wtLmat(caseind);
    jobparams.wt_el_pen        = wtEpmat(caseind);
    
    % Grid spacing in X and Y in um (for Cosserat model of the surrounding
    % region)
    jobparams.gridspacing_X    = gridspacing_X;
    jobparams.gridspacing_Y    = gridspacing_Y;
    fprintf(1,'Using\n    DelX: %2.2f\n    DelY: %2.2f\n',...
        jobparams.gridspacing_X,jobparams.gridspacing_Y);
    
    jobparams.Fmaxbound         = 500; % Upper bound on the maximum force value
    jobparams.uselinear         = 0;
    jobparams.solutionmethod    = 'fminunc'; %'gradientdescent'
    
    % Gradient descent parameters
    jobparams.Niter_actctr      = 300;
    jobparams.viscosity         = 0.5;  % viscosity (1/learning rate), increase if diverging solution
    jobparams.alpha             = 0.85; % momentum parameter
    jobparams.graddescentmethod = 'momentum';
    jobparams.plotinterm        = 0;
    
    % fmincon parameters
    jobparams.optimAlgorithm = 'quasi-newton';
    jobparams.optimStepTol   = 1E-4;
    jobparams.optimOptTol    = 1E-5;
    jobparams.optimMaxIters  = 200;
    jobparams.optimMaxFevals = 500;
    
    % Active contour model parameters
    jobparams.E_actctr         = E_actctrmat(caseind);
    jobparams.dia_actctr       = dia_actctrmat(caseind);
    jobparams.prestrain_actctr = 0.0;
    
    Idatafolder = sprintf('%s/ActiveContourImages',parentfolder);
    if ~isfolder(Idatafolder)
       mkdir(Idatafolder);
    end
    
    jobparams.frameNum              = framemat(caseind);
    jobparams.datafoldername        = framefolder;
    jobparams.filtfluoreImageFile   = sprintf('%s/FilteredFluroImage_%02d.jpg', jobparams.datafoldername, jobparams.frameNum);

    cellmaskfile = sprintf('%s/%02d_cellROI.mat',jobparams.datafoldername,jobparams.frameNum);
    cellmask = load(cellmaskfile);
    jobparams.maskcoords = cellmask.hroipoly.Position;
    
    jobrundirprefix = sprintf('%s/withAC_',jobparams.datafoldername);


    if length(regcoeffmat) > 1
        error('regcoeffmat should just be one value');
    end

    regcoeffmat_ind = 1;
    jobparams.regcoeff      = regcoeffmat(regcoeffmat_ind);
    fprintf(1,'\n\n********************\n\nRunning analysis for \nreg_param =  %2.3E\n\n',jobparams.regcoeff);
    currregparamdir = sprintf('%sregcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,jobrundirprefix,jobparams.regcoeff,jobparams.E_actctr,jobparams.dia_actctr,jobparams.wt_nodalerror,jobparams.wt_line,jobparams.wt_el_pen);

    if exist(currregparamdir,'dir') ~= 7
        fprintf(1,'Generating output folder in %s\n',currregparamdir);
        mkdir(currregparamdir);
    end
        
    if useprevGamma_set == 1
        if caseind > 1
            % Generate Gammainit and restart solution files from previous
            % step converged solution
            prevframefolder     = sprintf(['%s/',framefolderformat],parentfolder,framemat(caseind-1));
            prevjobrundirprefix = sprintf('%s/withAC_',prevframefolder);
            prevregparamdir =  sprintf('%sregcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
                ,prevjobrundirprefix,jobparams.regcoeff,jobparams.E_actctr,jobparams.dia_actctr,jobparams.wt_nodalerror,jobparams.wt_line,jobparams.wt_el_pen);

            fprintf(1,'Loading restart files from\n    %s\n',prevregparamdir);

            % Load solution
            % Unknown forces nodes
            prevunNODEs = load(sprintf('%s/n2_unNodes.txt',prevframefolder));
            prevunDOFs = [(3*(prevunNODEs-1) +1) ; (3*(prevunNODEs-1) +2)]; 
            prevGammamat = load(sprintf('%s/restart_Gamma.txt',prevregparamdir));
            
            unNODEs = load(sprintf('%s/n2_unNodes.txt',framefolder));
            unDOFs = [(3*(unNODEs-1) +1) ; (3*(unNODEs-1) +2)];
            Gammainit = zeros(length(unDOFs),1);
            for i = 1:length(unDOFs)
                matchind = find(prevunDOFs==unDOFs(i));
                if ~isempty(matchind)
                    Gammainit(i) = prevGammamat(matchind);
                end
            end

            % Initializing solver files
            copyfile(sprintf('%s/restart_u_init.txt',prevregparamdir),sprintf('%s/restart_u_init.txt',currregparamdir));
            fp = fopen(sprintf('%s/restart_Gamma_init.txt',currregparamdir),'w+');
            fprintf(fp,'%2.5E\n',Gammainit);
            fclose(fp);

            % Restart solver files
            copyfile(sprintf('%s/restart_u.txt',prevregparamdir),sprintf('%s/restart_u.txt',currregparamdir));
            fp = fopen(sprintf('%s/restart_Gamma.txt',currregparamdir),'w+');
            fprintf(fp,'%2.5E\n',Gammainit);
            fclose(fp);

            jobparams.restart_gamma = 1;
        else
            jobparams.restart_gamma = 0;
        end
    end

    if jobparams.restart_gamma==1 && ~isfile(sprintf('%s/restart_Gamma.txt',currregparamdir))
        fprintf(1,'Setting jobparams.restart_gamma = 0\nCould not find restart files:\n%s\n',sprintf('%s/restart_Gamma.txt',currregparamdir));
        jobparams.restart_gamma = 0;
    end

    jobparams.outputfoldername      = currregparamdir;

    fprintf(1,'Running analysis in %s for frame %d\nOutput folder: %s\n\n',parentfolder, jobparams.frameNum,jobparams.outputfoldername);

    regNanonetTFM_Hiperlife(jobparams);
    fprintf(1,'Total run time: %2.2f mins\n',toc(totalruntime)/60);
end

% delete(poolobj);
