%% -------------------- Documentation --------------------
% This script processes nuclear reactor data for the specified CASE and BURNUP.
% It performs the following tasks:

% 1. Sets options for data processing, including fitting parameters and plotting options.
CASE = "R4C38";
BURNUP = "EOC_25_percent_long";
OUTPUTBASEPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\results";
INPUTBASEPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\input";
MSCNPPPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\reference";
opts.skip_mscnpp = 50;
opts.skip_xerom = 50;
opts.jump_xerom =5;
opts.maxPeaks = 6;

opts.plotting.LineWidth=2;
opts.plotting.FontSize=16;


opts.powerplot.n_zero_xerom = 1;
opts.powerplot.n_zero_mscnpp = 1;
opts.powerplot.n_peak_xerom = 1;
opts.powerplot.n_peak_mscnpp = 2;

opts.AOplot.n_zero_xerom = 1;
opts.AOplot.n_zero_mscnpp = 1;
opts.AOplot.n_peak_xerom = 1;
opts.AOplot.n_peak_mscnpp = 2;

opts.modeplot.n_zero_xerom = 1;
opts.modeplot.n_zero_mscnpp = 1;
opts.modeplot.n_peak_xerom = 1;
opts.modeplot.n_peak_mscnpp = 2;

opts.powerfit_xerom.f0= "Default";
opts.powerfit_xerom.A0 = "Default"; 
opts.powerfit_xerom.alpha0 = -0.01;
opts.powerfit_xerom.phi0 = "Default";
opts.powerfit_xerom.c0 = "Default";
opts.powerfit_xerom.c1 = "Default";
opts.powerfit_xerom.showPlot = 1;

opts.powerfit_mscnpp.f0= "Default";
opts.powerfit_mscnpp.A0 = "Default"; 
opts.powerfit_mscnpp.alpha0 = -0.03;
opts.powerfit_mscnpp.phi0 = "Default";
opts.powerfit_mscnpp.c0 = "Default";
opts.powerfit_mscnpp.c1 = "Default";
opts.powerfit_mscnpp.showPlot = 1;

opts.AOfit_xerom.f0= "Default";
opts.AOfit_xerom.A0 = "Default"; 
opts.AOfit_xerom.alpha0 = -0.03;
opts.AOfit_xerom.phi0 = pi/2;
opts.AOfit_xerom.c0 = "Default";
opts.AOfit_xerom.c1 = "Default";
opts.AOfit_xerom.showPlot = 1;

opts.AOfit_mscnpp.f0= "Default";
opts.AOfit_mscnpp.A0 = "Default"; 
opts.AOfit_mscnpp.alpha0 = -0.03;
opts.AOfit_mscnpp.phi0 = "Default";
opts.AOfit_mscnpp.c0 = "Default";
opts.AOfit_mscnpp.c1 = "Default";
opts.AOfit_mscnpp.showPlot = 1;

opts.modefit_xerom.f0= "Default";
opts.modefit_xerom.A0 = "Default"; 
opts.modefit_xerom.alpha0 = -0.01;
opts.modefit_xerom.phi0 = pi/2;
opts.modefit_xerom.c0 = "Default";
opts.modefit_xerom.c1 = "Default";
opts.modefit_xerom.showPlot = false;

opts.modefit_mscnpp.f0= "Default";
opts.modefit_mscnpp.A0 = "Default"; 
opts.modefit_mscnpp.alpha0 = -0.01;
opts.modefit_mscnpp.phi0 = "Default";
opts.modefit_mscnpp.c0 = "Default";
opts.modefit_mscnpp.c1 = "Default";
opts.modefit_mscnpp.showPlot = false;

%% Set up paths for input and output directories
inputDir = fullfile(INPUTBASEPATH, CASE, BURNUP ,"Refinement1/");
outputDir = fullfile(OUTPUTBASEPATH, CASE, BURNUP );% 2. Sets up paths for input and output directories.
referenceDir = fullfile(MSCNPPPATH, BURNUP,"results_R4C38_Refinement1_EOC_25_percent_long.mat");
%% 3. Loads the necessary data from the specified directories.
data = loadData(inputDir,outputDir,referenceDir);
%% 4. Modifies the loaded data according to the specified options.
data = modifyData(data,opts);
%% 5. Computes delta flux, iodine, and xenon concentrations.
data = computeDeltaFlux(data);
%% 6. Calculates delta power based on the computed fluxes.
data = preComputeDpower(data);
%% 7. Computes scaling factors for power and axial offset.
data = computePowerScaling(data,opts);
data = computePowerandAO(data);
%% 8. Prepares data for power comparisons.
data = preparePowerComparison(data,opts);
%% 9. Prepares data for AO comparisons. 
data = prepareAOComparison(data,opts);
%% 10. Prepares data for mode comparisons.
data = prepareModeComparison(data,opts);
%% 10. Generates plots for mode, power, and axial offset comparisons.
modeComparisonPlot(data,outputDir,opts);
%%
powerComparisonPlot(data,outputDir,opts);
%%
AOComparisonPlot(data,outputDir,opts);
% 11. Saves the results and returns the processed data.
