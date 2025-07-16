# Deep-learning-enabled live-cell fiber-force microscopy (DLFM)

This is a code to perform inverse analysis for force estimation for DLFM accompanying the paper **"Deep Learning Reveals How Cells Pull, Buckle, and Navigate Fibrous Environments"**.

## System requirements
The code is tested to run on *MATLAB R2022a* (installed on Windows 10). The following toolboxes are used: optimization, image processing.

## Installation guide
No additional installation other than Matlab (with the toolboxes mentioned above enabled) is necessary.

## Running instructions
The analysis is controlled through the configuration files stored in './src/configfiles'.

### Validation analysis
We have provided analysis performed for Fig. 4A as an example which uses configuration file 'runparameters_Fig4A_inversevalidation.m'. In each script mentioned below the configuration file should be set as follows:
```Matlab
%% Simulation parameters
runparameters_Fig4A_inversevalidation;
% runparameters_Fig2C;
```
1. Set the sigma_e value for the analysis by modifing this line in 'runparameters_Fig4A_inversevalidation.m'
```Matlab
% Set error introduced in the synthetic dataset
sigma_e           = 0.0; % 0.1, 0.2, 0.5
```
2. Plot the L-curve (Fig. 4Bi) by running 'plot_Lcurve.m' with appropriate configuration file.
3. Similarly, plot the error analysis of Fig. 4Bii by running 'plotError_validation_analysis.m' with this configuration file.
4. Run 'plotSolution_validation_analysis.m' with this configuration file to plot Fig. 4A.
5. To perform the inverse analysis, run 'runDLFM.m' with this configuration file. Note that on a typical desktop computer, it takes about 5-6 hours to sweep a range of regularization parameters. Since a nonlinear beam problem needs to be solved for each iteration, the simulations can be slower if errors and displacements are large.

### Experimental data analysis
We have provided analysis performed for Fig. 2C as an example which uses configuration file 'runparameters_Fig2C.m'. In all scripts mentiioned below use:
```Matlab
%% Simulation parameters
% runparameters_Fig4A_inversevalidation;
runparameters_Fig2C;
```
1. Run 'plotSolution.m' with the configuration file to generate Fig. 2C.
2. Plot the L-curve by running 'plot_Lcurve.m' with appropriate configuration file.
3. To perform the inverse analysis, run 'runDLFM.m' with this configuration file.

## Input and Output data
The data needed to perform the analysis is stored in './src/testcases'.

## Contact
Please contact Sohan Kale (kale@vt.edu) if there are any issues producing the expected results.
