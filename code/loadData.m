function [data] = loadData(inputDir, resultsDir, refPath)
    data = struct();
    filePath = fullfile(inputDir,"GEOM_data.mat");
    s = load(filePath, "DX","DY","DZ"); 
    [data.DX,data.DY,data.DZ] = deal(s.DX,s.DY,s.DZ);
    
    filePath = fullfile(resultsDir,"RESULTS_HET.mat");
    s = load(filePath,"state_values_2G","time_2G");
    data.xerom.state_values_2G = s.state_values_2G;
    data.xerom.time_2G = s.time_2G;

    filePath = fullfile(inputDir,"RESULTS.mat");
    s = load(filePath,"MOD1","MOD2","keff");
    [data.xerom.MOD1,data.xerom.MOD2,data.keff] = deal(s.MOD1,s.MOD2,s.keff);

    filePath = fullfile(inputDir,"XS_data.mat");
    s = load(filePath,"NUFIS1","NUFIS2");
    [data.NUFIS1,data.NUFIS2] = deal(s.NUFIS1,s.NUFIS2);

    filePath = fullfile(inputDir,"POWER_data.mat");
    s = load(filePath,"KAPPA1","KAPPA2","REFERENCE_POWER");
    [data.KAPPA1,data.KAPPA2,data.mscnpp.power] = deal(s.KAPPA1,s.KAPPA2,s.REFERENCE_POWER);
    
    s = load("../input/CONSTANTS_data.mat","NU"); data.NU = s.NU;
    filePath = fullfile(resultsDir, "PARAMETERS_data.mat");
    s = load(filePath,"M","PHID_PHI","PHID_F_PHI");
    [data.M, data.xerom.PHID_PHI,data.xerom.PHID_F_PHI] = deal(s.M,s.PHID_PHI,s.PHID_F_PHI);
    
    s = load(refPath,"FLX1","FLX2","Xe135","I135","inputData","t");
    [data.mscnpp.FLX1,data.mscnpp.FLX2,data.mscnpp.Xe135,data.mscnpp.I135,data.mscnpp.InputData,data.mscnpp.t]...
        = deal(s.FLX1,s.FLX2,s.Xe135,s.I135,s.inputData,s.t);
end
