function Calculate_feedback_coefficients(input_dir, feedback_dir)
%% load data for heterogeneous feedback coefficients

%Prepare input data
fprintf("Preparing input data for feedback term calculation \n")
load("../input/CONSTANTS_data.mat","NU");
power_input = load(input_dir+"POWER_data.mat");
FUM_input = load(input_dir+"FUM_data.mat","FLX1","FLX2","keff");
XS_input = load(input_dir+"XS_data");
HELPER_input = load(input_dir+"HELPER_data.mat");
STEAMING_input = load(input_dir+"STREAMING_MATRICES_data.mat");
GEOM_input = load(input_dir+"GEOM_data.mat");

input = struct();
input = concatenate_structs(input,power_input);
input = concatenate_structs(input,FUM_input);
input = concatenate_structs(input,XS_input);
input = concatenate_structs(input,HELPER_input);
input = concatenate_structs(input,STEAMING_input);
input = concatenate_structs(input,GEOM_input);
input.KFIS1 = input.KAPPA1.*input.NUFIS1/NU;
input.KFIS2 = input.KAPPA2.*input.NUFIS2/NU;
clear *_input input.KAPPA1 input.KAPPA2
fprintf("Input data collected \n")


%Prepare feedback data
fprintf("Preparing feedback data for feedback term calculation \n")
power_feedback = load(feedback_dir+"POWER_data.mat");
FUM_feedback = load(feedback_dir+"FUM_data.mat","FLX1","FLX2","keff");
XS_feedback = load(feedback_dir+"XS_data");
HELPER_feedback = load(feedback_dir+"HELPER_data.mat");
STEAMING_feedback = load(feedback_dir+"STREAMING_MATRICES_data.mat");
GEOM_feedback = load(feedback_dir+"GEOM_data.mat");

feedback = struct();
feedback = concatenate_structs(feedback,power_feedback);
feedback = concatenate_structs(feedback,FUM_feedback);
feedback = concatenate_structs(feedback,XS_feedback);
feedback = concatenate_structs(feedback,HELPER_feedback);
feedback = concatenate_structs(feedback,STEAMING_feedback);
feedback = concatenate_structs(feedback,GEOM_feedback);
feedback.KFIS1 = feedback.KAPPA1.*feedback.NUFIS1/NU;
feedback.KFIS2 = feedback.KAPPA2.*feedback.NUFIS2/NU;
clear *_feedback feedback.KAPPA1 feedback.KAPPA2 NU
fprintf("Feedback data collected \n")


DV = input.DX*input.DY*input.DZ;

sizex = size(input.ABS1,1);
sizey = size(input.ABS1,2);
sizez = size(input.ABS1,3);

SHIFT = input.SHIFT;
SHIFT_XYZ = input.SHIFT_XYZ;
SHIFT_XY = input.SHIFT_XY;
I_MAX = input.I_MAX;
J_MAX = input.J_MAX;
K_MAX = input.K_MAX;
TYP = input.TYP;
CONV = input.CONV;
input.MOD_EQ = abs([input.FLX1(:,:,:,1);input.FLX2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
input.KFISINT =  DV*sum(G2_inner_product([input.KFIS1,input.KFIS2],input.MOD_EQ,"vector","vector"),"all");
input.PS = input.REFERENCE_POWER*input.keff/input.KFISINT;

feedback.MOD_EQ = abs([feedback.FLX1(:,:,:,1);feedback.FLX2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
feedback.KFISINT = DV*sum(G2_inner_product([feedback.KFIS1,feedback.KFIS2],feedback.MOD_EQ,"vector","vector"),"all");
feedback.PS = feedback.REFERENCE_POWER*feedback.keff/feedback.KFISINT;

input.MOD_EQ_scaled= input.PS*input.MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the REFERENCE_POWER
input.MOD_EQ_1_scaled= input.MOD_EQ_scaled(1:sizex,:,:);
input.MOD_EQ_2_scaled= input.MOD_EQ_scaled(sizex+1:end,:,:);

feedback.MOD_EQ_scaled= feedback.PS*feedback.MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the REFERENCE_POWER
feedback.MOD_EQ_1_scaled= feedback.MOD_EQ_scaled(1:sizex,:,:);
feedback.MOD_EQ_2_scaled= feedback.MOD_EQ_scaled(sizex+1:end,:,:);


input.STA_FLX_col=zeros(1,SHIFT_XYZ*2);
feedback.STA_FLX_col=zeros(1,SHIFT_XYZ*2);
% Calculate streaming terms
for I=1:I_MAX
    for J=1:J_MAX
        for K=1:K_MAX
            if (TYP(I,J,K)~=0)
                % Group 1
                input.STA_FLX_col(CONV(I,J,K)) = input.MOD_EQ_1_scaled(I,J,K);
                feedback.STA_FLX_col(CONV(I,J,K)) = feedback.MOD_EQ_1_scaled(I,J,K);
                % Group 2
                input.STA_FLX_col(CONV(I,J,K)+SHIFT_XYZ) = input.MOD_EQ_2_scaled(I,J,K);
                feedback.STA_FLX_col(CONV(I,J,K)+SHIFT_XYZ) = feedback.MOD_EQ_2_scaled(I,J,K);
            end
        end
    end
end


input.NapJ = (input.AX/input.DX+input.AY/input.DY+input.AZ/input.DZ+input.BX/input.DX+input.BY/input.DY+input.BZ/input.DZ+input.CX/input.DX+input.CY/input.DY+input.CZ/input.DZ) * input.STA_FLX_col';
feedback.NapJ = (feedback.AX/feedback.DX+feedback.AY/feedback.DY+feedback.AZ/feedback.DZ+feedback.BX/feedback.DX+feedback.BY/feedback.DY+feedback.BZ/feedback.DZ+feedback.CX/feedback.DX+feedback.CY/feedback.DY+feedback.CZ/feedback.DZ) * feedback.STA_FLX_col';


for I = 1:I_MAX
    for J = 1:J_MAX
        for K = 1:K_MAX
            input.NJ1(I,J,K) = input.NapJ(CONV(I,J,K)+SHIFT_XY);
            feedback.NJ1(I,J,K) = feedback.NapJ(CONV(I,J,K)+SHIFT_XY);
            input.NJ2(I,J,K) = input.NapJ(CONV(I,J,K)+SHIFT_XYZ);
            feedback.NJ2(I,J,K) = feedback.NapJ(CONV(I,J,K)+SHIFT_XYZ);
        end
    end
end

%Calculate input reaction rates 
input.NUFIS1PHI1 = input.NUFIS1.*input.MOD_EQ_1_scaled; 
input.NUFIS2PHI2 = input.NUFIS2.*input.MOD_EQ_2_scaled;
input.ABS1PHI1 = input.ABS1.*input.MOD_EQ_1_scaled;
input.ABS2PHI2 = input.ABS2.*input.MOD_EQ_2_scaled;
input.REMPHI1 = input.REM.*input.MOD_EQ_1_scaled;

%Calculate feedback reaction rates 
feedback.NUFIS1PHI1 = feedback.NUFIS1.*feedback.MOD_EQ_1_scaled;
feedback.NUFIS2PHI2 = feedback.NUFIS2.*feedback.MOD_EQ_2_scaled;
feedback.ABS1PHI1 = feedback.ABS1.*feedback.MOD_EQ_1_scaled;
feedback.ABS2PHI2 = feedback.ABS2.*feedback.MOD_EQ_2_scaled;
feedback.REMPHI1 = feedback.REM.*feedback.MOD_EQ_1_scaled;

%Calculate reaction rate differences
DNapJ = feedback.NapJ - input.NapJ;
DNUFIS1PHI1 = feedback.NUFIS1PHI1 - input.NUFIS1PHI1;
DNUFIS2PHI2 = feedback.NUFIS2PHI2 - input.NUFIS2PHI2;  
DABS1PHI1 = feedback.ABS1PHI1 - input.ABS1PHI1;
DABS2PHI2 = feedback.ABS2PHI2 - input.ABS2PHI2;
DREMPHI1 = feedback.REMPHI1 - input.REMPHI1;

%Calculate flux differences
DFLX1 = feedback.MOD_EQ_1_scaled - input.MOD_EQ_1_scaled;
DFLX2 = feedback.MOD_EQ_2_scaled - input.MOD_EQ_2_scaled;

%convert back to 3D array
DNJ1 = zeros(I_MAX,J_MAX,K_MAX);
DNJ2 = zeros(I_MAX,J_MAX,K_MAX);
for I = 1:I_MAX
    for J = 1:J_MAX
        for K = 1:K_MAX
            DNJ1(I,J,K) = DNapJ(CONV(I,J,K)+SHIFT_XY);
            DNJ2(I,J,K) = DNapJ(CONV(I,J,K)+SHIFT_XYZ);
        end
    end
end
%Divide change in reaction rates by changes in the flux
DABS1PHI1DPHI1 = DABS1PHI1 ./DFLX1;
DNUFIS1PHI1DPHI1 = DNUFIS1PHI1 ./ DFLX1;
DNapJ1DPHI1 = DNJ1 ./ DFLX1;
DREMPHI1DPHI1 = DREMPHI1 ./ DFLX1;
DABS2PHI2DPHI2 = DABS2PHI2./ DFLX2;
DNapJDPHI2 = DNJ2 ./ DFLX2;
DNUFIS2PHI2DPHI2 = DNUFIS2PHI2 ./ DFLX2;

%Remove infinities arising from zeros in the flux difference variable
DABS1PHI1DPHI1(isinf(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isinf(DNUFIS1PHI1DPHI1)) = 0;
DNapJ1DPHI1(isinf(DNapJ1DPHI1)) = 0;
DREMPHI1DPHI1(isinf(DREMPHI1DPHI1)) = 0;
DABS2PHI2DPHI2(isinf(DABS2PHI2DPHI2)) = 0;
DNapJDPHI2(isinf(DNapJDPHI2)) = 0;
DNUFIS2PHI2DPHI2(isinf(DNUFIS2PHI2DPHI2)) = 0;

%Remove any NaN values
DABS1PHI1DPHI1(isnan(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isnan(DNUFIS1PHI1DPHI1)) = 0;
DNapJ1DPHI1(isnan(DNapJ1DPHI1)) = 0;
DREMPHI1DPHI1(isnan(DREMPHI1DPHI1)) = 0;
DABS2PHI2DPHI2(isnan(DABS2PHI2DPHI2)) = 0;
DNapJDPHI2(isnan(DNapJDPHI2)) = 0;
DNUFIS2PHI2DPHI2(isnan(DNUFIS2PHI2DPHI2)) = 0;

%Calculate the matrix elements of the heterogeneous feedback matrix
K11 = DNUFIS1PHI1DPHI1-DREMPHI1DPHI1-DNapJ1DPHI1-DABS1PHI1DPHI1;
K12 = DNUFIS2PHI2DPHI2;
K21 = DREMPHI1DPHI1;
K22 = -DNapJDPHI2-DABS2PHI2DPHI2;

%Save the final matrix elements
save(input_dir+"FEEDBACK_data.mat","K11","K12","K21","K22")
end
function struct_predator = concatenate_structs(struct_predator,struct_prey)
    for fn = fieldnames(struct_prey)'
        struct_predator.(fn{1}) = struct_prey.(fn{1});
    end
end