clear variables; clc; close all;
load("input/test_data.mat")


if exist('output','dir')~=7
    mkdir('output')
end
if exist('feedback','dir')~=7
    mkdir('feedback')
end
if exist('input','dir')~=7
    mkdir('input')
end


%%
D1 = Case_1_D1;
D2 = Case_1_D2;
NUFIS1 = Case_1_NF1;
NUFIS2 = Case_1_NF2;
ABS1 = Case_1_SA1;
ABS2 = Case_1_SA2;
REM = Case_1_SS12;
% Replace NaN values with zeroes
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KF1 = Case_1_KF1;
KF1(isnan(KF1)) = 0;
KF2 = Case_2_KF2;
KF2(isnan(KF2)) = 0;
STA_FLX1 = Case_1_FLX1;
STA_FLX1(isnan(STA_FLX1)) = 0;
STA_FLX2 = Case_1_FLX2;
STA_FLX2(isnan(STA_FLX2)) = 0;

save("input/XS_data.mat","REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")
clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")

%% Create feedback data

%%
D1 = Case_3_D1;
D2 = Case_3_D2;
NUFIS1 = Case_3_NF1;
NUFIS2 = Case_3_NF2;
ABS1 = Case_3_SA1;
ABS2 = Case_3_SA2;
REM = Case_3_SS12;
% Replace NaN values with zeroes
D1(isnan(D1)) = 0;
D2(isnan(D2)) = 0;
NUFIS1(isnan(NUFIS1)) = 0;
NUFIS2(isnan(NUFIS2)) = 0;
ABS1(isnan(ABS1)) = 0;
ABS2(isnan(ABS2)) = 0;
REM(isnan(REM)) = 0;
KF1 = Case_3_KF1;
KF1(isnan(KF1)) = 0;
KF2 = Case_2_KF2;
KF2(isnan(KF2)) = 0;
STA_FLX1 = Case_3_FLX1;
STA_FLX1(isnan(STA_FLX1)) = 0;
STA_FLX2 = Case_3_FLX2;
STA_FLX2(isnan(STA_FLX2)) = 0;

save("feedback/XS_data.mat","REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")
clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1","KF1","KF2","STA_FLX1","STA_FLX2")

Assembly_pitch = 21.5036; %cm
Node_height = 15.2400; % cm
DX = Assembly_pitch/2;
DY = Assembly_pitch/2;
DZ = Node_height;
save("input/GEOM_data.mat","DX","DY","DZ")
save("feedback/GEOM_data.mat","DX","DY","DZ")
clear("DZ","DY","DX")

clear variables