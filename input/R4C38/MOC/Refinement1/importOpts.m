function opts = importOpts(CASE, BURNUP)
%Verify CASE
% Validate CASE and BURNUP inputs
if isempty(CASE) || isempty(BURNUP)
    error('CASE and BURNUP must be provided.');
end
opts.path.CASE = "R4C38";
opts.path.BURNUP = "MOC";
if opts.path.CASE ~= CASE || opts.path.BURNUP ~= BURNUP
    error("Specified case or burnup does not match option case or option burnup! \n")
else
    fprintf("Loading options for Case %s %s \n",CASE,BURNUP)
end
opts.path.OUTPUTBASEPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\results";
opts.path.INPUTBASEPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\input";
opts.path.MSCNPPPATH = "C:\Users\amand\Documents\Kristoffer\XEROM\reference";
opts.skip_mscnpp = 1;
opts.skip_xerom = 50;
opts.jump_xerom =5;
opts.maxPeaks = 10;

opts.overwrite_data = false;

opts.plotting.LineWidth=2;
opts.plotting.FontSize=16;

opts.powerplot.n_zero_xerom = 1;
opts.powerplot.n_zero_mscnpp = 1;
opts.powerplot.n_peak_xerom = 1;
opts.powerplot.n_peak_mscnpp = 1;

opts.AOplot.n_zero_xerom = 1;
opts.AOplot.n_zero_mscnpp = 1;
opts.AOplot.n_peak_xerom = 1;
opts.AOplot.n_peak_mscnpp = 1;

opts.modeplot.n_zero_xerom = 1;
opts.modeplot.n_zero_mscnpp = 1;
opts.modeplot.n_peak_xerom = 1;
opts.modeplot.n_peak_mscnpp = 1;

opts.xenonplot.n_zero_xerom = 1;
opts.xenonplot.n_zero_mscnpp = 1;

opts.iodineplot.n_zero_xerom = 1;
opts.iodineplot.n_zero_mscnpp = 1;

opts.xenonplot.n_peak_xerom = 1;
opts.xenonplot.n_peak_mscnpp = 1;

opts.iodineplot.n_peak_xerom = 1;
opts.iodineplot.n_peak_mscnpp = 1;

opts.powerfit_xerom.skip = 50;
opts.powerfit_xerom.f0= "Default";
opts.powerfit_xerom.fU= "Default";
opts.powerfit_xerom.fL= "Default";
opts.powerfit_xerom.A0 = "Default"; 
opts.powerfit_xerom.alpha0 = "Default";
opts.powerfit_xerom.alphaU = "Default";
opts.powerfit_xerom.alphaL = "Default";
opts.powerfit_xerom.phi0 = "Default";
opts.powerfit_xerom.c0 = "Default";
opts.powerfit_xerom.c1 = "Default";
opts.powerfit_xerom.showPlot = false;
opts.powerfit_xerom.name = "Power XEROM fit";

opts.powerfit_mscnpp.skip = 10;
opts.powerfit_mscnpp.f0= 0.05;
opts.powerfit_mscnpp.fU= 0.06;
opts.powerfit_mscnpp.fL= 0.035;
opts.powerfit_mscnpp.A0 = "Default"; 
opts.powerfit_mscnpp.alpha0 = -0.05;
opts.powerfit_mscnpp.alphaU = -0.07;
opts.powerfit_mscnpp.alphaL = -0.085;
opts.powerfit_mscnpp.phi0 = "Default";
opts.powerfit_mscnpp.c0 = "Default";
opts.powerfit_mscnpp.c1 = "Default";
opts.powerfit_mscnpp.showPlot = true;
opts.powerfit_mscnpp.name = "Power MScNPP fit";

opts.AOfit_xerom.skip = 50; %not used
opts.AOfit_xerom.f0= "Default";
opts.AOfit_xerom.fU= "Default";
opts.AOfit_xerom.fL= "Default";
opts.AOfit_xerom.A0 = "Default"; 
opts.AOfit_xerom.alpha0 = "Default";
opts.AOfit_xerom.alphaU = "Default";
opts.AOfit_xerom.alphaL = "Default";
opts.AOfit_xerom.phi0 = "Default";
opts.AOfit_xerom.c0 = "Default";
opts.AOfit_xerom.c1 = "Default";
opts.AOfit_xerom.showPlot = true;
opts.AOfit_xerom.name = "AO XEROM fit";

opts.AOfit_mscnpp.skip = 10; %not used
opts.AOfit_mscnpp.f0= "Default";
opts.AOfit_mscnpp.fU= "Default";
opts.AOfit_mscnpp.fL= "Default";
opts.AOfit_mscnpp.A0 = "Default"; 
opts.AOfit_mscnpp.alpha0 = "Default";
opts.AOfit_mscnpp.alphaU = "Default";
opts.AOfit_mscnpp.alphaL = "Default";
opts.AOfit_mscnpp.phi0 = "Default";
opts.AOfit_mscnpp.c0 = "Default";
opts.AOfit_mscnpp.c1 = "Default";
opts.AOfit_mscnpp.showPlot = true;
opts.AOfit_mscnpp.name = "AO MScNPP fit";

opts.modefit_xerom.skip = 50; %not used
opts.modefit_xerom.f0= "Default";
opts.modefit_xerom.fU= "Default";
opts.modefit_xerom.fL= "Default";
opts.modefit_xerom.A0 = "Default"; 
opts.modefit_xerom.alpha0 = "Default";
opts.modefit_xerom.alphaU = "Default";
opts.modefit_xerom.alphaL = "Default";
opts.modefit_xerom.phi0 = "Default";
opts.modefit_xerom.c0 = "Default";
opts.modefit_xerom.c1 = "Default";
opts.modefit_xerom.showPlot = true;
opts.modefit_xerom.name = "Mode XEROM fit";

opts.modefit_mscnpp.skip = 1; %not used
opts.modefit_mscnpp.f0= "Default";
opts.modefit_mscnpp.fU= "Default";
opts.modefit_mscnpp.fL= "Default";
opts.modefit_mscnpp.A0 = "Default"; 
opts.modefit_mscnpp.alpha0 = "Default";
opts.modefit_mscnpp.alphaU = "Default";
opts.modefit_mscnpp.alphaL = "Default";
opts.modefit_mscnpp.phi0 = "Default";
opts.modefit_mscnpp.c0 = "Default";
opts.modefit_mscnpp.c1 = "Default";
opts.modefit_mscnpp.showPlot = true;
opts.modefit_mscnpp.name = "Mode MScNPP fit";

opts.xenonfit_xerom.skip = 1; %not used
opts.xenonfit_xerom.f0= "Default";
opts.xenonfit_xerom.fU= "Default";
opts.xenonfit_xerom.fL= "Default";
opts.xenonfit_xerom.A0 = "Default"; 
opts.xenonfit_xerom.alpha0 = "Default";
opts.xenonfit_xerom.alphaU = "Default";
opts.xenonfit_xerom.alphaL = "Default";
opts.xenonfit_xerom.phi0 = "Default";
opts.xenonfit_xerom.c0 = "Default";
opts.xenonfit_xerom.c1 = "Default";
opts.xenonfit_xerom.showPlot = true;
opts.xenonfit_xerom.name = "Xenon XEROM fit";

opts.iodinefit_xerom.skip = 40; %not used
opts.iodinefit_xerom.f0= "Default";
opts.iodinefit_xerom.fU= "Default";
opts.iodinefit_xerom.fL= "Default";
opts.iodinefit_xerom.A0 = "Default"; 
opts.iodinefit_xerom.alpha0 = "Default";
opts.iodinefit_xerom.alphaU = "Default";
opts.iodinefit_xerom.alphaL = "Default";
opts.iodinefit_xerom.phi0 = "Default";
opts.iodinefit_xerom.c0 = "Default";
opts.iodinefit_xerom.c1 = "Default";
opts.iodinefit_xerom.showPlot = true;
opts.iodinefit_xerom.name = "Iodine XEROM fit";

opts.xenonfit_mscnpp.skip = 15; %not used
opts.xenonfit_mscnpp.f0= "Default";
opts.xenonfit_mscnpp.fU= "Default";
opts.xenonfit_mscnpp.fL= "Default";
opts.xenonfit_mscnpp.A0 = "Default"; 
opts.xenonfit_mscnpp.alpha0 = "Default";
opts.xenonfit_mscnpp.alphaU = "Default";
opts.xenonfit_mscnpp.alphaL = "Default";
opts.xenonfit_mscnpp.phi0 = "Default";
opts.xenonfit_mscnpp.c0 = "Default";
opts.xenonfit_mscnpp.c1 = "Default";
opts.xenonfit_mscnpp.showPlot = true;
opts.xenonfit_mscnpp.name = "Xenon MScNPP fit";

opts.iodinefit_mscnpp.skip = 15; %not used
opts.iodinefit_mscnpp.f0= "Default";
opts.iodinefit_mscnpp.fU= "Default";
opts.iodinefit_mscnpp.fL= "Default";
opts.iodinefit_mscnpp.A0 = "Default"; 
opts.iodinefit_mscnpp.alpha0 = "Default";
opts.iodinefit_mscnpp.alphaU = "Default";
opts.iodinefit_mscnpp.alphaL = "Default";
opts.iodinefit_mscnpp.phi0 = "Default";
opts.iodinefit_mscnpp.c0 = "Default";
opts.iodinefit_mscnpp.c1 = "Default";
opts.iodinefit_mscnpp.showPlot = true;
opts.iodinefit_mscnpp.name = "Iodine MScNPP fit";

end