clear variables; clc; close all;
load("input/test_data.mat")


if exist('output','dir')~=7
    mkdir('output')
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
KF2 = Case_2_KF2;
STA_FLX1 = Case_1_FLX1;
STA_FLX2 = Case_1_FLX2;

save("input/XS_data.mat","REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1")
clear("REM","ABS2","ABS1","NUFIS2","NUFIS1","D2","D1")
Assembly_pitch = 21.5036; %cm
Node_height = 15.2400; % cm
DX = Assembly_pitch/2;
DY = Assembly_pitch/2;
DZ = Node_height;
save("input/GEOM_data.mat","DX","DY","DZ")
clear("DZ","DY","DX")

input_XS.D1 =Case_1_D1;
input_XS.D2 = Case_1_D2;
input_XS.NF1 = Case_1_NF1;
input_XS.NF2 = Case_1_NF2;
input_XS.SA1 = Case_1_SA1;
input_XS.SA2 = Case_1_SA2;
input_XS.REM = Case_1_SS12;
input_XS.STA_FLX1 = Case_1_FLX1;
input_XS.STA_fLX2 = Case_1_FLX2;

feedback_XS.D1 =Case_3_D1;
feedback_XS.D2 = Case_3_D2;
feedback_XS.NF1 = Case_3_NF1;
feedback_XS.NF2 = Case_3_NF2;
feedback_XS.SA1 = Case_3_SA1;
feedback_XS.SA2 = Case_3_SA2;
feedback_XS.REM = Case_3_SS12;
feedback_XS.STA_FLX1 = Case_3_FLX1;
feedback_XS.STA_fLX2 = Case_3_FLX2;

save("input\XEROM_data.mat","KF1","KF2","input_XS","feedback_XS")
clear variables