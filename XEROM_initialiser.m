clear variables; clc; close all;
load("input/intrusive1_Tkorr_caseAconv1.mat")


if exist('output','dir')~=7
    mkdir('output')
end
if exist('feedback','dir')~=7
    mkdir('feedback')
end
if exist('input','dir')~=7
    mkdir('input')
end

Assembly_pitch = 21.5036; %cm
Node_height = 15.2400; % cm
DX = Assembly_pitch/2;
DY = Assembly_pitch/2;
DZ = Node_height;

%%
D1 = Case_2_D1;
D1(:,:,1) = Case_1_D1(:,:,1);
D2 = Case_2_D2;
NUFIS1 = Case_2_NF1;
NUFIS2 = Case_2_NF2;
ABS1 = Case_2_SA1;
ABS2 = Case_2_SA2;
REM = Case_2_SS12;
% Replace NaN values with zeroes
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KF1 = Case_2_KF1;
KF1(isnan(KF1)) = 0;
KF2 = Case_2_KF2;
KF2(isnan(KF2)) = 0;
STA_FLX1 = Case_2_FLX1;
STA_FLX1(isnan(STA_FLX1)) = 0;
STA_FLX2 = Case_2_FLX2;
STA_FLX2(isnan(STA_FLX2)) = 0;
reactor_power = 3287.6699E6 * 99.8685/100; % W
sigmaX = 2.7000e-18; %cm^2
lambdaI = 2.87e-5; % s^-1
lambdaX = 2.09e-5; %s^-1
gammaI = 0.062; %branching ration iodine
gammaX = 0.002; %branching ratio xenon

REFINE = 1;

if REFINE
    ABS1 = repelem(ABS1,2,2,2);
    ABS2= repelem(ABS2,2,2,2);
    NUFIS1= repelem(NUFIS1,2,2,2);
    NUFIS2= repelem(NUFIS2,2,2,2);
    REM = repelem(REM,2,2,2);
    D1= repelem(D1,2,2,2);
    D2= repelem(D2,2,2,2);
    STA_FLX1 = repelem(STA_FLX1,2,2,2);
    STA_FLX2 = repelem(STA_FLX2,2,2,2);
    KF1 = repelem(KF1,2,2,2);
    KF2 = repelem(KF2,2,2,2);
    DX = DX/2;
    DY = DY/2;
    DZ = DZ/2;
end



save("input/XS_data_new.mat","REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2")
save("input/additional_data_new","STA_FLX1","STA_FLX2","reactor_power","gammaX","gammaI","lambdaX","lambdaI","sigmaX")
clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")

%% Create feedback data

%%
D1 = Case_4_D1;
D1(:,:,1) = Case_3_D1(:,:,1);
D2 = Case_4_D2;
NUFIS1 = Case_4_NF1;
NUFIS2 = Case_4_NF2;
ABS1 = Case_4_SA1;
ABS2 = Case_4_SA2;
REM = Case_4_SS12;
KF1 = Case_4_KF1;
KF2 = Case_4_KF2;
STA_FLX1 = Case_4_FLX1;
STA_FLX2 = Case_4_FLX2;
% Replace NaN values with zeroes
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KF1(isnan(KF1)) = 0;
KF2(isnan(KF2)) = 0;
STA_FLX1(isnan(STA_FLX1)) = 0;
STA_FLX2(isnan(STA_FLX2)) = 0;

if REFINE
    ABS1 = repelem(ABS1,2,2,2);
    ABS2= repelem(ABS2,2,2,2);
    NUFIS1= repelem(NUFIS1,2,2,2);
    NUFIS2= repelem(NUFIS2,2,2,2);
    REM = repelem(REM,2,2,2);
    D1= repelem(D1,2,2,2);
    D2= repelem(D2,2,2,2);
    STA_FLX1 = repelem(STA_FLX1,2,2,2);
    STA_FLX2 = repelem(STA_FLX2,2,2,2);
    KF1 = repelem(KF1,2,2,2);
    KF2 = repelem(KF2,2,2,2);
end

save("feedback/XS_data_new.mat","REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2")
save("feedback/additional_data_new","STA_FLX1","STA_FLX2","reactor_power","gammaX","gammaI","lambdaX","lambdaI","sigmaX")
clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")


save("input/GEOM_data.mat","DX","DY","DZ")
save("feedback/GEOM_data.mat","DX","DY","DZ")
clear("DZ","DY","DX")

clear variables