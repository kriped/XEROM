function XEROM_initialiser(CASE,BURNUP)

% NAME DIRECTORIES
SIMULATE_INPUT_DATA_PATH =       sprintf("../input/%s/%s/SIMULATE_data/Input_data.mat",CASE,BURNUP);
XEROM_INPUT_REFINEMENT0_dir =    sprintf("../input/%s/%s/Refinement0/",CASE,BURNUP);
XEROM_INPUT_REFINEMENT1_dir =    sprintf("../input/%s/%s/Refinement1/",CASE,BURNUP);
XEROM_FEEDBACK_REFINEMENT0_dir = sprintf("../input/%s/%s/feedback/Refinement0/",CASE,BURNUP);
XEROM_FEEDBACK_REFINEMENT1_dir = sprintf("../input/%s/%s/feedback/Refinement1/",CASE,BURNUP);
%SIMULATE_REFERENCE_D1_PATH =     sprintf("../input/%s/%s/SIMULATE_data/D1/Reference_data.mat",CASE,BURNUP);
%SIMULATE_REFERENCE_D2_PATH =     sprintf("../input/%s/%s/SIMULATE_data/D2/Reference_data.mat",CASE,BURNUP);
CONSTANTS_PATH =                 sprintf("../input/CONSTANTS_data.mat");
%% OPTIONS TO OVERWRITE DATA
OVERWRITE_EIGENFUNCTIONS = false;
OVERWRITE_FEEDBACK = true;
OVERWRITE_STREAMING = false;
OVERWRITE_PRECURSORS = false;
OVERWRITE_RESULTS = false;
%% Initialise and load directories
if exist(SIMULATE_INPUT_DATA_PATH,"file")
    load(SIMULATE_INPUT_DATA_PATH,"*Case_2*","Case_4*")
else
    error("SIMULATE input files not found")
end
if ~exist(XEROM_FEEDBACK_REFINEMENT0_dir,"dir")
    mkdir(XEROM_FEEDBACK_REFINEMENT0_dir)
end
if ~exist(XEROM_FEEDBACK_REFINEMENT1_dir,"dir")
    mkdir(XEROM_FEEDBACK_REFINEMENT1_dir)
end
if ~exist(XEROM_INPUT_REFINEMENT0_dir,"dir")
    mkdir(XEROM_INPUT_REFINEMENT0_dir)
end
if ~exist(XEROM_INPUT_REFINEMENT1_dir,"dir")
    mkdir(XEROM_INPUT_REFINEMENT1_dir)
end
if exist(CONSTANTS_PATH,"file")
    load(CONSTANTS_PATH,"NU")
else
    error("Constants not found")
end

% DEFINE NODE DIMENSIONS
Assembly_pitch = 21.5036; %cm
Node_height = 15.2400; % cm
DX = Assembly_pitch/2;
DY = Assembly_pitch/2;
DZ = Node_height;
%% EXTRACT UNREFINED CROSS SECTIONS
D1 = Case_2_D1;
D2 = Case_2_D2;
NUFIS1 = Case_2_NF1;
NUFIS2 = Case_2_NF2;
ABS1 = Case_2_SA1;
ABS2 = Case_2_SA2;
REM = Case_2_SS12;
KF1 = Case_2_KF1;
KF2 = Case_2_KF2;
KAPPA1 = NU*KF1./NUFIS1;
KAPPA2 = NU*KF2./NUFIS2;
REFERENCE_POWER = Case_2_Power*1e6; % W

% REPLACE NAN WITH ZEROS
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KAPPA1(isnan(KAPPA1)) = 0;
KAPPA2(isnan(KAPPA2)) = 0;

% SAVE UNREFINED DATA
save(XEROM_INPUT_REFINEMENT0_dir+"GEOM_data","DZ","DY","DX")
save(XEROM_INPUT_REFINEMENT0_dir+"XS_data","ABS1","ABS2","NUFIS1","NUFIS2","REM","D1","D2")
save(XEROM_INPUT_REFINEMENT0_dir+"POWER_data","KAPPA1","KAPPA2","REFERENCE_POWER")

%REFINE CROSS SECTIONS
ABS1 = repelem(ABS1,2,2,2);
ABS2= repelem(ABS2,2,2,2);
NUFIS1= repelem(NUFIS1,2,2,2);
NUFIS2= repelem(NUFIS2,2,2,2);
REM = repelem(REM,2,2,2);
D1= repelem(D1,2,2,2);
D2= repelem(D2,2,2,2);
KAPPA1 = repelem(KAPPA1,2,2,2);
KAPPA2 = repelem(KAPPA2,2,2,2);
DX = DX/2;
DY = DY/2;
DZ = DZ/2;

% SAVE REFINED DATA
save(XEROM_INPUT_REFINEMENT1_dir+"GEOM_data","DZ","DY","DX")
save(XEROM_INPUT_REFINEMENT1_dir+"XS_data","ABS1","ABS2","NUFIS1","NUFIS2","REM","D1","D2")
save(XEROM_INPUT_REFINEMENT1_dir+"POWER_data","KAPPA1","KAPPA2","REFERENCE_POWER")

clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","REFERENCE_POWER")

%% Create feedback data
% DEFINE NODE DIMENSIONS
Assembly_pitch = 21.5036; %cm
Node_height = 15.2400; % cm
DX = Assembly_pitch/2;
DY = Assembly_pitch/2;
DZ = Node_height;
% EXTRACT UNREFINED CROSS SECTIONS
D1 = Case_4_D1;
D2 = Case_4_D2;
NUFIS1 = Case_4_NF1;
NUFIS2 = Case_4_NF2;
ABS1 = Case_4_SA1;
ABS2 = Case_4_SA2;
REM = Case_4_SS12;
KF1 = Case_4_KF1;
KF2 = Case_4_KF2;
KAPPA1 = NU*KF1./NUFIS1;
KAPPA2 = NU*KF2./NUFIS2;

%EXTRACT POWER DATA
REFERENCE_POWER = Case_4_Power*1e6; % W

% REPLACE ZEROS WITH NAN
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KAPPA1(isnan(KAPPA1)) = 0;
KAPPA2(isnan(KAPPA2)) = 0;

%SAVE UNREFINED CROSS SECTIONS
save(XEROM_FEEDBACK_REFINEMENT0_dir+"GEOM_data","DZ","DY","DX")
save(XEROM_FEEDBACK_REFINEMENT0_dir+"XS_data","ABS1","ABS2","NUFIS1","NUFIS2","REM","D1","D2")
save(XEROM_FEEDBACK_REFINEMENT0_dir+"POWER_data","KAPPA1","KAPPA2","REFERENCE_POWER")

%REFINE CROSS SECTIONS
ABS1 = repelem(ABS1,2,2,2);
ABS2= repelem(ABS2,2,2,2);
NUFIS1= repelem(NUFIS1,2,2,2);
NUFIS2= repelem(NUFIS2,2,2,2);
REM = repelem(REM,2,2,2);
D1= repelem(D1,2,2,2);
D2= repelem(D2,2,2,2);
KAPPA1 = repelem(KAPPA1,2,2,2);
KAPPA2 = repelem(KAPPA2,2,2,2);
DX=DX/2;
DY=DY/2;
DZ=DZ/2;

%SAVE REFINED CROSS SECTIONS
save(XEROM_FEEDBACK_REFINEMENT1_dir+"GEOM_data","DZ","DY","DX")
save(XEROM_FEEDBACK_REFINEMENT1_dir+"XS_data","ABS1","ABS2","NUFIS1","NUFIS2","REM","D1","D2")
save(XEROM_FEEDBACK_REFINEMENT1_dir+"POWER_data","KAPPA1","KAPPA2","REFERENCE_POWER")

%clearvars -except Assembly_pitch Node_height SIMULATE_REFERENCE_D1_PATH SIMULATE_REFERENCE_D2_PATH
%% Calculate Eigenmodes
% CALCULATE input EIGENMODES FOR UNREFINED DATA
if ~exist(XEROM_INPUT_REFINEMENT0_dir+"RESULTS.mat","file") || OVERWRITE_EIGENFUNCTIONS
    fprintf("Calculating Eigenfunctions for %s %s unrefined data\n", CASE, BURNUP)
    Calculate_eigenfunction(XEROM_INPUT_REFINEMENT0_dir) 
end

if ~exist(XEROM_FEEDBACK_REFINEMENT0_dir+"RESULTS.mat","file") || OVERWRITE_EIGENFUNCTIONS
    fprintf("Calculating Eigenfunctions for %s %s unrefined data\n", CASE, BURNUP)
    Calculate_eigenfunction(XEROM_FEEDBACK_REFINEMENT0_dir) 
end
% SPLIT THE REFINED RESULTS.MAT 
files = ["STREAMING_MATRICES_data.mat","HELPER_data.mat"];

% Check if any file in the list does not exist
missing_file = false;

for i = 1:length(files)
    if ~exist(fullfile(XEROM_INPUT_REFINEMENT1_dir, files(i)), "file")
        missing_file = true;
        break; % Stop checking as we found a missing file
    end
end

% Run function if any file is missing
if missing_file || OVERWRITE_STREAMING
    Calculate_Streaming_matrices(XEROM_INPUT_REFINEMENT1_dir);
end

% Check if any file in the list does not exist
missing_file = false;

for i = 1:length(files)
    if ~exist(fullfile(XEROM_FEEDBACK_REFINEMENT1_dir, files(i)), "file")
        missing_file = true;
        break; % Stop checking as we found a missing file
    end
end

% Run function if any file is missing
if missing_file || OVERWRITE_STREAMING
    Calculate_Streaming_matrices(XEROM_FEEDBACK_REFINEMENT1_dir);
end

if ~exist(XEROM_INPUT_REFINEMENT1_dir+"FUM_data.mat","file")|| OVERWRITE_RESULTS
    Result_splitting(XEROM_INPUT_REFINEMENT1_dir)
end

if ~exist(XEROM_FEEDBACK_REFINEMENT1_dir+"FUM_data.mat","file")|| OVERWRITE_RESULTS
    Result_splitting(XEROM_FEEDBACK_REFINEMENT1_dir)
end
%% Calculate feedback terms

% CALCULATE FEEDBACK TERMS FOR UNREFINED DATA
if ~exist(XEROM_INPUT_REFINEMENT0_dir+"FEEDBACK_data.mat",'file') || OVERWRITE_FEEDBACK
    fprintf("Calculating Feedback terms for for %s %s unrefined data\n", CASE, BURNUP)
    Calculate_feedback_coefficients(XEROM_INPUT_REFINEMENT0_dir,XEROM_FEEDBACK_REFINEMENT0_dir)
end
% CALCULATE FEEDBACK TERMS FOR REFINED DATA
if ~exist(XEROM_INPUT_REFINEMENT1_dir+"FEEDBACK_data.mat",'file') || OVERWRITE_FEEDBACK
    fprintf("Calculating Feedback terms for for %s %s refined data\n", CASE, BURNUP)
    Calculate_feedback_coefficients(XEROM_INPUT_REFINEMENT1_dir,XEROM_FEEDBACK_REFINEMENT1_dir)
end

%% Create feedback data
if ~exist(XEROM_INPUT_REFINEMENT0_dir+"PRECURSORS_data.mat","file") ||~exist(XEROM_INPUT_REFINEMENT1_dir+"PRECURSORS_data.mat","file") || OVERWRITE_PRECURSORS
    fprintf("Please input the relevant precursor data for %s %s \n", CASE, BURNUP)
    [lambda,beta] = take_precursor_data(CASE,BURNUP); 
    save(XEROM_INPUT_REFINEMENT0_dir+"PRECURSORS_data.mat","beta","lambda")
    fprintf("Precursor data saved to unrefined dir \n")
    save(XEROM_INPUT_REFINEMENT1_dir+"PRECURSORS_data.mat","beta","lambda")
    fprintf("Precursor data saved to refined dir \n")
end
%% Create reference data
REF = 0;

if REF
    for ref_case = ["D1","D2"]
        load(sprintf("../%s/SIMULATE_data/%s/Reference_data.mat",input_dir,ref_case))

        DX_ref = Assembly_pitch/2;
        DY_ref = Assembly_pitch/2;
        DZ_ref = Node_height;
        DV_ref = DX_ref*DY_ref*DZ_ref;
        FLX1_EQ_ref = Case_1_FLX_1;
        FLX1_EQ_ref(isnan(FLX1_EQ_ref))=0; % Remove Nan values
        FLX2_EQ_ref = Case_1_FLX_2;
        FLX2_EQ_ref(isnan(FLX2_EQ_ref))=0; % Remove Nan values
        X0_ref = Case_1_XENON;
        X0_ref(isnan(X0_ref)) = 0;
        I0_ref = Case_1_IODINE;
        I0_ref(isnan(I0_ref)) = 0;
        Power_ref = Case_1_Power*1e6; % Convert power from MW to W
        Power_rel_ref = Case_1_Rel_Power;
        KF1_ref = Case_1_KF1;
        KF1_ref(isnan(KF1_ref)) = 0;
        KF2_ref = Case_1_KF2;
        KF2_ref(isnan(KF2_ref)) = 0;
        FLX1_ref = Case_5_output_FLX_1;
        FLX1_ref(isnan(FLX1_ref)) = 0;
        FLX2_ref = Case_5_output_FLX_2;
        FLX2_ref(isnan(FLX2_ref)) = 0;
        X_ref = Case_5_output_XENON*1e24;
        X_ref(isnan(X_ref)) = 0;
        I_ref = Case_5_output_IODINE*1e24;
        I_ref(isnan(I_ref)) = 0;

        if ~exist(sprintf("../%s/%s",reference_dir,ref_case),"dir")
            mkdir(sprintf("../%s/%s",reference_dir,ref_case))
        end

        save(sprintf("../%s/%s/reference_time_series.mat",reference_dir,ref_case), "I_ref","X_ref","FLX1_ref","FLX2_ref","KF1_ref","KF2_ref","FLX1_EQ_ref","FLX2_EQ_ref","X0_ref","I0_ref","Power_rel_ref","Power_ref","DV_ref","DX_ref","DY_ref","DZ_ref")
    end
end

fprintf("All data for %s - %s has been initialised! \n",CASE,BURNUP)
end