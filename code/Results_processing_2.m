%% -------------------- Documentation --------------------
% This script processes nuclear reactor data for the specified CASE and BURNUP.
% It performs the following tasks:

% 1. Sets options for data processing, including fitting parameters and plotting options.
BURNUP = "UNSTABLE";                                                                                                                                                              
CASE   = "R4C38";
opts_path = fullfile("C:\Users\amand\Documents\Kristoffer\XEROM\input",CASE,BURNUP,"Refinement1/");
addpath(opts_path);
opts = importOpts(CASE, BURNUP);
rmpath(opts_path);
close all
%% Set up paths for input and output directories
inputDir    = fullfile(opts.path.INPUTBASEPATH, opts.path.CASE, opts.path.BURNUP ,"Refinement1/");
outputDir   = fullfile(opts.path.OUTPUTBASEPATH, opts.path.CASE, opts.path.BURNUP );% 2. Sets up paths for input and output directories.
referenceDir= fullfile(opts.path.MSCNPPPATH, opts.path.CASE,opts.path.BURNUP,opts.path.CASE+"_"+opts.path.BURNUP+".mat");
%% 3. Loads the necessary data from the specified directories.
fprintf("Now running loadData.m\n")
data        = loadData(inputDir,outputDir,referenceDir,opts);
%% 4. Modifies the loaded data according to the specified options.
fprintf("Now running modifyData.m\n")
data        = modifyData(data,opts);
%% 5. Computes delta flux, iodine, and xenon concentrations.
fprintf("Now running computeDeltaFlux.m\n")
data        = computeDeltaFlux(data);
%% 5. Computes delta iodine, and xenon concentrations.
fprintf("Now running computeIodineXenon.m\n")
data        = computeIodineXenon(data);
%% 6. Calculates delta power based on the computed fluxes.
fprintf("Now running preComputeDpower\n")
data        = preComputeDpower(data);
%% 7. Computes scaling factors for power and axial offset.
fprintf("Now running computePowerScaling.m\n")
data        = computePowerScaling(data,opts);
fprintf("Now running computePowerandAO.m\n")
data        = computePowerandAO(data);
%% 8. Prepares data for power comparisons.
fprintf("Now running preparePowerComparison.m\n")
data        = preparePowerComparison(data,opts);
%% 9. Prepares data for AO comparisons. 
fprintf("Now running prepareAOComparison.m\n")
data        = prepareAOComparison(data,opts);
% %% 10. Prepares data for mode comparisons.
% fprintf("Now running prepareModeComparison.m\n")
% data        = prepareModeComparison(data,opts);
%% 11. Prepares data for Iodine and Xenon plotting.
fprintf("Now running prepareIodineXenoncomparison\n")
data        = prepareIodineXenoncomparison(data,opts);
% %% 12. Generates plots for all modes
% fprintf("Now running modeComparisonPlot.m\n")
% modeComparisonPlot(data,outputDir,opts);
%% 13. Generates Plots for Power comparison
fprintf("Now running powerComparisonPlot.m\n")
powerComparisonPlot(data,outputDir,opts);
%% 14. Generates Plots for AO comparison
fprintf("Now running AOComparisonPlot.m\n")
AOComparisonPlot(data,outputDir,opts);
%% 15. Generates Plots for Iodine comparison
fprintf("Now running IodineComparisonPlot.m\n")
IodineComparisonPlot(data,outputDir,opts)
%% 16. Generates Plots for Xenon comparison
fprintf("Now running XenonComparisonPlot.m\n")
XenonComparisonPlot(data,outputDir,opts)
