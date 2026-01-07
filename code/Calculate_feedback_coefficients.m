function Calculate_feedback_coefficients(input_dir, feedback_dir)
%% load data for heterogeneous feedback coefficients
close all
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
DFLX1 = (feedback.MOD_EQ_1_scaled - input.MOD_EQ_1_scaled);
DFLX2 = (feedback.MOD_EQ_2_scaled - input.MOD_EQ_2_scaled);

%DFLX1=threshold_filter_3d(DFLX1,1e+10);
%DFLX2=threshold_filter_3d(DFLX2,1e+10);
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
 %DNJ1=block_sum_3d(DNJ1);
 %DNJ2=block_sum_3d(DNJ2);

%% Investigate delta terms

DABS1PHI1 = average_assemblies(DABS1PHI1);
DABS2PHI2 = average_assemblies(DABS2PHI2);
DNUFIS1PHI1 = average_assemblies(DNUFIS1PHI1);
DNUFIS2PHI2 = average_assemblies(DNUFIS2PHI2);
DNJ1 = average_assemblies(DNJ1);
DNJ2 = average_assemblies(DNJ2);
DREMPHI1 = average_assemblies(DREMPHI1);
DFLX1 = average_assemblies(DFLX1);
DFLX2 = average_assemblies(DFLX2);
% 
% %% Calculate errors
% 
% DABS1PHI1_RRMSE = RRMSE(DABS1PHI1_avg,DABS1PHI1)*100;
% DABS2PHI2_RRMSE = RRMSE(DABS2PHI2_avg,DABS2PHI2)*100;
% DNUFIS1PHI1_RRMSE = RRMSE(DNUFIS1PHI1_avg,DNUFIS1PHI1)*100;
% DNUFIS2PHI2_RRMSE = RRMSE(DNUFIS2PHI2_avg,DNUFIS2PHI2)*100;
% DNJ1_RRMSE = RRMSE(DNJ1_avg,DNJ1)*100;
% DNJ2_RRMSE = RRMSE(DNJ2_avg,DNJ2)*100;
% DREMPHI1_RRMSE = RRMSE(DREMPHI1_avg,DREMPHI1)*100;
% DFLX1_RRMSE = RRMSE(DFLX1_avg,DFLX1)*100;
% DFLX2_RRMSE = RRMSE(DFLX2_avg,DFLX2)*100;
% 
% sprintf("Delta ABS1 RR average RRMSE = %.1f%%",DABS1PHI1_RRMSE)
% sprintf("Delta ABS2 RR average RRMSE = %.1f%%",DABS2PHI2_RRMSE)
% sprintf("Delta NUFIS1 RR average RRMSE = %.1f%%",DNUFIS1PHI1_RRMSE)
% sprintf("Delta NUFIS2 RR average RRMSE = %.1f%%",DNUFIS2PHI2_RRMSE)
% sprintf("Delta FAST CURRENT average RRMSE = %.1f%%",DNJ1_RRMSE)
% sprintf("Delta THERMAL CURRENT average RRMSE = %.1f%%",DNJ2_RRMSE)
% sprintf("Delta REM RR average RRMSE = %.1f%%",DREMPHI1_RRMSE)
% sprintf("Delta FLX1 average RRMSE = %.1f%%",DFLX1_RRMSE)
% sprintf("Delta FLX2 average RRMSE = %.1f%%",DFLX2_RRMSE)
% 
% RRMSE_table_headers = ["Delta ABS1 RR"; "Delta ABS2 RR"; "Delta NUFIS1 RR"; "Delta NUFIS2 RR"; "Delta FAST CURRENT"; "Delta THERMAL CURRENT"; "Delta REM RR"; "Delta FLX1"; "Delta FLX2"];
% RRMSE_values = [DABS1PHI1_RRMSE;DABS2PHI2_RRMSE;DNUFIS1PHI1_RRMSE;DNUFIS2PHI2_RRMSE;DNJ1_RRMSE;DNJ2_RRMSE;DREMPHI1_RRMSE;DFLX1_RRMSE;DFLX2_RRMSE];
% RRMSE_table = table(RRMSE_table_headers,RRMSE_values, 'VariableNames', {'Variable', 'RRMSE(%)'});
% if size(DABS1PHI1,1) == 34
%     writetable(RRMSE_table, '../input/R3C40/EOC/Refinement0/RRMSE_table.csv');
% else
%     writetable(RRMSE_table, '../input/R3C40/EOC/Refinement1/RRMSE_table.csv');
% end
% 
% L2norm = @(averaged,original) norm(averaged - original,'fro') / norm(original,'fro');
% 
% DABS1PHI1_L2norm = L2norm(DABS1PHI1_avg,DABS1PHI1)*100;
% DABS2PHI2_L2norm = L2norm(DABS2PHI2_avg,DABS2PHI2)*100;
% DNUFIS1PHI1_L2norm = L2norm(DNUFIS1PHI1_avg,DNUFIS1PHI1)*100;
% DNUFIS2PHI2_L2norm = L2norm(DNUFIS2PHI2_avg,DNUFIS2PHI2)*100;
% DNJ1_L2norm = L2norm(DNJ1_avg,DNJ1)*100;
% DNJ2_L2norm = L2norm(DNJ2_avg,DNJ2)*100;
% DREMPHI1_L2norm = L2norm(DREMPHI1_avg,DREMPHI1)*100;
% DFLX1_L2norm = L2norm(DFLX1_avg,DFLX1)*100;
% DFLX2_L2norm = L2norm(DFLX2_avg,DFLX2)*100;
% 
% mre = @(averaged,original) abs(averaged - original) ./ abs(original);
% 
% DABS1PHI1_mre = mre(DABS1PHI1_avg,DABS1PHI1)*100;
% DABS2PHI2_mre = mre(DABS2PHI2_avg,DABS2PHI2)*100;
% DNUFIS1PHI1_mre = mre(DNUFIS1PHI1_avg,DNUFIS1PHI1)*100;
% DNUFIS2PHI2_mre = mre(DNUFIS2PHI2_avg,DNUFIS2PHI2)*100;
% DNJ1_mre = mre(DNJ1_avg,DNJ1)*100;
% DNJ2_mre = mre(DNJ2_avg,DNJ2)*100;
% DREMPHI1_mre = mre(DREMPHI1_avg,DREMPHI1)*100;
% DFLX1_mre = mre(DFLX1_avg,DFLX1)*100;
% DFLX2_mre = mre(DFLX2_avg,DFLX2)*100;

%Divide change in reaction rates by changes in the flux
DABS1PHI1DPHI1 = DABS1PHI1 ./ DFLX1;
DNUFIS1PHI1DPHI1 = DNUFIS1PHI1 ./ DFLX1;
DNJ1DPHI1 = DNJ1 ./ DFLX1;
DREMPHI1DPHI1 = DREMPHI1 ./ DFLX1;
DABS2PHI2DPHI2 = DABS2PHI2./ DFLX2;
DNJ2DPHI2 = DNJ2 ./ DFLX2;
DNUFIS2PHI2DPHI2 = DNUFIS2PHI2 ./ DFLX2;


%Remove infinities arising from zeros in the flux difference variable
DABS1PHI1DPHI1(isinf(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isinf(DNUFIS1PHI1DPHI1)) = 0;
DNJ1DPHI1(isinf(DNJ1DPHI1)) = 0;
DREMPHI1DPHI1(isinf(DREMPHI1DPHI1)) = 0;
DABS2PHI2DPHI2(isinf(DABS2PHI2DPHI2)) = 0;
DNJ2DPHI2(isinf(DNJ2DPHI2)) = 0;
DNUFIS2PHI2DPHI2(isinf(DNUFIS2PHI2DPHI2)) = 0;

%Remove any NaN values
DABS1PHI1DPHI1(isnan(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isnan(DNUFIS1PHI1DPHI1)) = 0;
DNJ1DPHI1(isnan(DNJ1DPHI1)) = 0;
DREMPHI1DPHI1(isnan(DREMPHI1DPHI1)) = 0;
DABS2PHI2DPHI2(isnan(DABS2PHI2DPHI2)) = 0;
DNJ2DPHI2(isnan(DNJ2DPHI2)) = 0;
DNUFIS2PHI2DPHI2(isnan(DNUFIS2PHI2DPHI2)) = 0;

%Calculate the matrix elements of the heterogeneous feedback matrix
K11 = DNUFIS1PHI1DPHI1-DREMPHI1DPHI1-DNJ1DPHI1-DABS1PHI1DPHI1;
K12 = DNUFIS2PHI2DPHI2;
K21 = DREMPHI1DPHI1;
K22 = -DNJ2DPHI2-DABS2PHI2DPHI2;

%K11(:,:,end) = K11(:,:,end-1);
%K12(:,:,end) = K12(:,:,end-1);
%K21(:,:,end) = K21(:,:,end-1);
%K22(:,:,end) = K22(:,:,end-1);
% K11 = threshold_filter_3d(K11,0.05,'a');
% K12 = threshold_filter_3d(K12,0.38,"a");
% K21 = threshold_filter_3d(K21,0.08,"a");
% K22 = threshold_filter_3d(K22,5,"a");
% K22 = threshold_filter_3d(K22,2.5,"a");
% K22 = threshold_filter_3d(K22,1,"a");
% K22 = threshold_filter_3d(K22,0.5,"a");
% K22 = threshold_filter_3d(K22,0.25,"a");
% K22 = threshold_filter_3d(K22,0.2,"a");
% K22 = threshold_filter_3d(K22,0.2,"a");
%View4D(abs(K22))



%Save the final matrix elements
save(input_dir+"FEEDBACK_data.mat","K11","K12","K21","K22")
end
function struct_predator = concatenate_structs(struct_predator,struct_prey)
    for fn = fieldnames(struct_prey)'
        struct_predator.(fn{1}) = struct_prey.(fn{1});
    end
end

function B = block_sum_3d(A)
    if sum(size(A)) == sum([68,68,52]) 
        divider = 4;
    elseif sum(mod(size(A),2)) == 0
        divider = 2;
    else 
        error("all matrix sizes must be even numbers")
    end
    % Get the original dimensions
    [m, n, p] = size(A);
    shift = divider-1;
    % Create output matrix (half the size in each dimension)
    B = zeros(m/divider, n/divider, p/divider);
    
    % Sum blocks of 2x2x2 elements
    for i = 1:m/divider
        for j = 1:n/divider
            for k = 1:p/divider
                % Sum the 2x2x2 block starting at (2i-1, 2j-1, 2k-1)
                B(i,j,k) = sum(sum(sum(A(divider*i-shift:divider*i, divider*j-shift:divider*j, divider*k-shift:divider*k))));
            end
        end
    end
end

function [B] = threshold_filter_3d(A, threshold,AB)
    % Simple threshold filter for 3D Matrices
    % Inputs:
    %   A - 3D array
    %   threshold - Minimum absolute value
    % Output:
    %   B - 3D array of filtered values where points above the threshold
    %   are substituted by mean values of valid neighboring points.
    
    % Create a copy of the input array
    B = A;
    
    switch AB
        case "a"
            % Find points where flux difference is above threshold
            [peak_idx, peak_values]=find_peaks_3d(abs(A),'MinPeakHeight',threshold,'NeighborhoodSize',[3,3,3]);
            small_diff_mask = zeros(size(A));
            linear_indices = sub2ind(size(small_diff_mask),peak_idx(:,1),peak_idx(:,2),peak_idx(:,3));
            small_diff_mask(linear_indices) = 1;
            display_results = @(name, indices, values) ...
            disp_top_peaks(name, indices, values, min(10, size(indices, 1)));
            display_results('VariableName', peak_idx, peak_values);
        case "b"
            % Find points where flux difference is below threshold
            small_diff_mask = abs(A) < threshold;
        otherwise
            error("State whether problematic values are found above (a) or below (b) the threshold")
    end
    
    % Process only if we have points below threshold
    if ~any(small_diff_mask(:))
        if AB == "a"
            fprintf("No problematic points were found above the threshold value")
        else
            fprintf("No problematic points were found below the threshold value")
        end
        return
    end

    % Loop through each problematic point
    for i = 1:length(peak_idx)
        x = peak_idx(i,1);
        y = peak_idx(i,2);
        z = peak_idx(i,3);

        % Define neighborhood bounds
        x_min = max(1, x-1);
        x_max = min(size(A, 1), x+1);
        y_min = max(1, y-1);
        y_max = min(size(A, 2), y+1);
        z_min = max(1, z-1);
        z_max = min(size(A, 3), z+1);

        % Get local neighborhood
        neighborhood = A(x_min:x_max, y_min:y_max, z_min:z_max);
        neighborhood_mask = ~small_diff_mask(x_min:x_max, y_min:y_max, z_min:z_max);

        % Extract valid neighboring values (above threshold)
        valid_values = neighborhood(neighborhood_mask);

        % Replace with interpolated value if valid neighbors exist
        if ~isempty(valid_values)
            B(x, y, z) = mean(valid_values(:),"omitmissing");
        else
            B(x, y, z) = threshold;
        end
    end
end

% Helper function to display peak results
function disp_top_peaks(name, indices, values, num)
    disp(['------- ', name, ' -------']);
    if isempty(indices)
        disp('No peaks found.');
        return;
    end
    
    for i = 1:num
        fprintf('Peak at [%d, %d, %d] with value %.4f\n', ...
            indices(i,1), indices(i,2), indices(i,3), values(i));
    end
    fprintf('Total peaks found: %d\n\n', size(indices, 1));
end

