%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLFM
% Refer to 'readme.md' file for running instructions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
addpath(genpath('./src'));

%% Simulation configuration file
runparameters_Fig2C;
% runparameters_Fig4A_inversevalidation;

%% Perform analysis
% Replace with 'parfor caseind = 1:Ncases' to run multiple cases parallely
totalruntime = tic;
for caseind = 1:Ncases
    frameNum         = framemat(caseind);
    parentfolder     = sprintf(parentfolderformat,cellIDmat(caseind));
    framefolder      = sprintf(['%s/',framefolderformat],parentfolder,frameNum);

    jobparams                   = []; % to be able to work with parfor

    % Plot and save intermediate converged step solutions
    jobparams.plotupdates        = plotupdates;

    % To copy last regcoefficient Gammamat files for restart
    useprevGamma_set            = 1;
    % To restart the minimization from existing Gammamat and u
    % restart_u.txt and restart_Gamma.txt
    jobparams.restart_gamma     = 1;

    % perturb initial guess for restarts
    jobparams.randomize_Gammainit = 0;
    jobparams.std_rand_Gammainit  = 0.05;
   
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
    
    % Minimization algorithm parameters
    jobparams.uselinear         = 0;
    jobparams.solutionmethod    = 'fminunc';
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
        
    for regcoeffmat_ind = 1:length(regcoeffmat)
        
        fprintf(1,'\n\n********************\n\nRunning analysis for \nreg_param =  %2.3E\n\n',regcoeffmat(regcoeffmat_ind));
        
        jobparams.regcoeff      = regcoeffmat(regcoeffmat_ind);
        currregparamdir = sprintf('%sregcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
            ,jobrundirprefix,jobparams.regcoeff,jobparams.E_actctr,jobparams.dia_actctr,jobparams.wt_nodalerror,jobparams.wt_line,jobparams.wt_el_pen);
        
        if useprevGamma_set == 1
            % First reg index does not have a previous Gamma solution
            if regcoeffmat_ind > 1
                % Copy Gammainit and restart solution files from previous
                % converged solution folder to the current folder
                prevregparamdir =  sprintf('%sregcoeff%2.3E_Eactctr%2.3E_dia%2.3E_wN%2.3E_wAC%2.3E_wER%2.3E'...
                    ,jobrundirprefix,regcoeffmat(regcoeffmat_ind-1),jobparams.E_actctr,jobparams.dia_actctr,jobparams.wt_nodalerror,jobparams.wt_line,jobparams.wt_el_pen);

                if exist(currregparamdir,'dir') ~= 7
                    fprintf(1,'Generating output folder in %s\n',currregparamdir);
                    mkdir(currregparamdir);
                end

                fprintf(1,'Loading restart files from\n    %s\n',prevregparamdir);

                % Initializing solver files
                copyfile(sprintf('%s/restart_u_init.txt',prevregparamdir),sprintf('%s/restart_u_init.txt',currregparamdir));
                copyfile(sprintf('%s/restart_Gamma_init.txt',prevregparamdir),sprintf('%s/restart_Gamma_init.txt',currregparamdir));

                % Restart solver files
                copyfile(sprintf('%s/restart_u.txt',prevregparamdir),sprintf('%s/restart_u.txt',currregparamdir));
                copyfile(sprintf('%s/restart_Gamma.txt',prevregparamdir),sprintf('%s/restart_Gamma.txt',currregparamdir));

                jobparams.restart_gamma = 1;
            end
        end

        if jobparams.restart_gamma==1 && ~isfile(sprintf('%s/restart_Gamma.txt',currregparamdir))
            fprintf(1,'Setting jobparams.restart_gamma = 0\nCould not find restart files:\n%s\n',sprintf('%s/restart_Gamma.txt',currregparamdir));
            jobparams.restart_gamma = 0;
        end
        
        jobparams.outputfoldername      = currregparamdir;
        
        fprintf(1,'Running analysis in %s for frame %d\nOutput folder: %s\n\n',parentfolder, jobparams.frameNum,jobparams.outputfoldername);
        
        % Perform inverse analysis and save solution
        regNanonetTFM_Hiperlife(jobparams);
        fprintf(1,'Total run time: %2.2f mins\n',toc(totalruntime)/60);
    end 
end

% Delete parpool if opened
if exist('poolobj','var')
    delete(poolobj);
end

%% Other runparameter files
% runparameters_Fig3C_Control_series;
% runparameters_Fig3C_Inhibition;
% runparameters_Fig3E_series;
% runparameters_Fig3D_Control_series;
% runparameters_Fig4A_Hybrid_Grid;
% runparameters_Fig4B_Cell_Division_series;
% runparameters_Fig4C_Cell_Differentiation_Day0;
% runparameters_Fig4C_Cell_Differentiation_Day3_Control;
% runparameters_Fig4C_Cell_Differentiation_Day3_Osteo;
% runparameters_Fig4C_Cell_Differentiation_Day3_Adipo;
% runparameters_Fig4C_Cell_Differentiation_Day6_Control;
% runparameters_Fig4C_Cell_Differentiation_Day6_Osteo;
% runparameters_Fig4C_Cell_Differentiation_Day6_Adipo;
% runparameters_Fig4C_Cell_Differentiation_Day9_Control;
% runparameters_Fig4C_Cell_Differentiation_Day9_Osteo;
% runparameters_Fig4C_Cell_Differentiation_Day9_Adipo;
