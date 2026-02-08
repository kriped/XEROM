%% -------------------- Documentation --------------------
% This script processes nuclear reactor data for the specified CASE and BURNUP.
% It performs the following tasks:

% 1. Sets options for data processing, including fitting parameters and plotting options.
opts = importOpts();
%% Set up paths for input and output directories
inputDir = fullfile(opts.path.INPUTBASEPATH, opts.path.CASE, opts.path.BURNUP ,"Refinement1/");
outputDir = fullfile(opts.path.OUTPUTBASEPATH, opts.path.CASE, opts.path.BURNUP );% 2. Sets up paths for input and output directories.
referenceDir = fullfile(opts.path.MSCNPPPATH, opts.path.BURNUP,"results_R4C38_Refinement1_EOC_25_percent_long.mat");
%% 3. Loads the necessary data from the specified directories.
data = loadData(inputDir,outputDir,referenceDir);
%% 4. Modifies the loaded data according to the specified options.
data = modifyData(data,opts);
%% 5. Computes delta flux, iodine, and xenon concentrations.
data = computeDeltaFlux(data);
%% 5. Computes delta iodine, and xenon concentrations.
data = computeIodineXenon(data);
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
%% 11. Prepares data for Iodineand Xenon plotting.
data = prepareIodineXenoncomparison(data,opts);
%% 12. Generates plots for all modes
modeComparisonPlot(data,outputDir,opts);
%% 13. Generates Plots for Power comparison
powerComparisonPlot(data,outputDir,opts);
%% 14. Generates Plots for AO comparison
AOComparisonPlot(data,outputDir,opts);
