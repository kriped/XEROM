function XEROM(CASE,BURNUP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          MAIN PROGRAM FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For explanatory notes about the use of this tool, see the file
% USERS_GUIDE.PDF in the directory "manuals".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT_REFINEMENT0_dir = sprintf("../input/%s/%s/Refinement0/",CASE,BURNUP);
INPUT_REFINEMENT1_dir = sprintf("../input/%s/%s/Refinement1/",CASE,BURNUP);
RESULTS_REFINEMENT0_dir = sprintf("../results/%s/%s/Refinement0/",CASE,BURNUP);
RESULTS_REFINEMENT1_dir = sprintf("../results/%s/%s/Refinement1/",CASE,BURNUP);

%Check data files are present
input_files = ["XS_data.mat","RESULTS.mat","POWER_data.mat","GEOM_data.mat","FEEDBACK_data.mat"];
for file = input_files
    if ~exist(sprintf(INPUT_REFINEMENT1_dir+file),"file")
        error("%s/%s was not found",INPUT_REFINEMENT1_dir,file)
    end
end

%% Load Xerom_data
load(INPUT_REFINEMENT1_dir+"GEOM_data.mat");
load(INPUT_REFINEMENT1_dir+"RESULTS.mat","MOD1","MOD2","MOD1_adj","MOD2_adj","lambda");
load(INPUT_REFINEMENT1_dir+"POWER_data.mat","KAPPA1","KAPPA2","REFERENCE_POWER");
load(INPUT_REFINEMENT1_dir+"XS_data.mat");
load(INPUT_REFINEMENT1_dir+"FEEDBACK_data.mat")
load("../input/CONSTANTS_data.mat");

%% Create vectors and matrices
M=size(MOD1,4); % number of modes (Change in SETTINGS.m and rerun calculate_eigenfunction)
sizex = size(MOD1,1);
sizey = size(MOD1,2);
sizez = size(MOD1,3);
INCLUDED_MODES = 1:M; %all modes included
%INCLUDED_MODES = []; %Select which modes to include;
MOD1 = MOD1(:,:,:,INCLUDED_MODES);
MOD2 = MOD2(:,:,:,INCLUDED_MODES);
MOD1_adj = MOD1_adj(:,:,:,INCLUDED_MODES);
MOD2_adj = MOD2_adj(:,:,:,INCLUDED_MODES);
power = REFERENCE_POWER; 
keff = lambda(1);
ZERO = zeros(size(MOD1_adj));

%FB = 1.3e-18;


V_inv = [1/v1, 0 ; 0, 1/v2];
K_VALUE = lambda; % K values / eigenvalues of the modes
DV = DX*DY*DZ; % Discreet volume element

SIGF1 = NUFIS1/NU; % Fast fission cross section
SIGF2 = NUFIS2/NU; % Thermal fission cross section
KFIS1 = KAPPA1.*SIGF1;
KFIS2 = KAPPA2.*SIGF2;
ZERO = zeros(size(NUFIS1)); % zero element matching the size of the reactor
ONE = ones(size(NUFIS1)); % unit element matching the size of the reactor
F = 1/keff.*[NUFIS1, NUFIS2;ZERO,ZERO]; % Fission matrix
MOD = [MOD1;MOD2]; % vector of solutions to the forward problem
MOD_adj = [MOD1_adj,MOD2_adj]; % vector of solutions to the adjoint problem
MOD_EQ = abs([MOD1(:,:,:,1);MOD2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
KFISINT =  DV*1/keff*sum(G2_inner_product([KFIS1,KFIS2],MOD_EQ,"vector","vector"),"all");
PS = power*keff/KFISINT;
MOD_EQ_scaled= PS*MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the power
MOD_EQ_1_scaled= MOD_EQ_scaled(1:sizey,:,:);
MOD_EQ_2_scaled= MOD_EQ_scaled(sizey+1:end,:,:);
MOD_UPPER = [MOD_EQ_2_scaled, ZERO ; ZERO, ZERO]; % Costom matrix used in the equations 
MOD_LOWER = [ZERO, ZERO; MOD_EQ_2_scaled, ZERO]; % Costom matrix used in the equations
%Make matrices of delta values only for termal feedbacks

%Original feedback matrix
DeltaXS = [K11,K12;K21,K22];

intial_AO = (DV*1/keff*sum(KFIS1(:,:,sizez/2+1:sizez).*MOD_EQ_1_scaled(:,:,sizez/2+1:sizez)+KFIS2(:,:,sizez/2+1:sizez).*MOD_EQ_2_scaled(:,:,sizez/2+1:sizez),'all')-DV*1/keff*sum(KFIS1(:,:,1:sizez/2).*MOD_EQ_1_scaled(:,:,1:sizez/2)+KFIS2(:,:,1:sizez/2).*MOD_EQ_2_scaled(:,:,1:sizez/2),'all'))/(1/keff*DV*sum(KFIS1(:,:,:).*MOD_EQ_1_scaled(:,:,:)+KFIS2(:,:,:).*MOD_EQ_2_scaled(:,:,:),'all'));
sprintf("Equilibrium axial offset = %.2f",intial_AO)
GAMMAI = 1/keff.*[gammaI*SIGF1,gammaI*SIGF2 ; ZERO, ZERO ]; % matrix containing the production of Iodine from fission
GAMMAX = 1/keff.*[gammaX*SIGF1,gammaX*SIGF2 ; ZERO, ZERO ]; % matrix containing the production of xenon from fission
VECXT = [ONE,ZERO]; % transformation vector

I0 = 1/(lambdaI)*G2_inner_product(VECXT, G2_inner_product(GAMMAI,MOD_EQ_scaled,"matrix","vector"),"vector","vector"); % I_0(r) = \hat(X)^T \cdot 1/(keff*lambda_I)Gamma_I\times Phi_0
X0 = ((gammaX+gammaI)*(SIGF1.*MOD_EQ_1_scaled+SIGF2.*MOD_EQ_2_scaled))./((lambdaX)+sigmaX*MOD_EQ_2_scaled); 
Xe_UPPER = [ZERO,X0;ZERO,ZERO]; % Custom matrix used in the equations

%% test properties
% close all
Volume = DV * sum(MOD1(:,:,:,1)~=0,'all');
active_volume = DV * sum(NUFIS1~=0,'all');
% SIGF_PHI = SIGF1.*MOD_EQ_scaled(1:sizex,:,:)+SIGF2.*MOD_EQ_scaled(sizex+1:end,:,:);
% SIGF_PHI_average = DV * sum(SIGF_PHI,'all')/active_volume;
% SIGF1_average = DV * sum(SIGF1,'all')/active_volume;
% SIGF2_average = DV * sum(SIGF2,'all')/active_volume;
% MOD_EQ1_average = DV * sum(MOD_EQ_scaled(1:sizex,:,:),'all')/Volume;
% MOD_EQ2_average = DV * sum(MOD_EQ_scaled(sizex+1:end,:,:),'all')/Volume;
% X0_Average = DV * sum(X0,'all')/active_volume;
% I0_Average = DV * sum(I0,'all')/active_volume;
% MOD_EQ2_max = max(MOD_EQ_scaled(sizex+1:end,:,:),[],'all');
% X0_max = max(X0,[],'all');
% I0_max = max(I0,[],'all');
% DSA1_average = DV*sum(DABS1PHI1DPHI1_CS,'all')/Volume
% DSF1_average = DV*sum(DNUFIS1PHI1DPHI1_CS,'all')/active_volume
% DNaJ1_average = DV*sum( DNapJ1DPHI1_CS,'all')/Volume
% DREM_average = DV*sum(DREMPHI1DPHI1_CS,'all')/Volume
% DSA2_average = DV*sum(DABS2PHI2DPHI2_CS,'all')/Volume
% DSF2_average = DV*sum(DNUFIS2PHI2DPHI2_CS,'all')/active_volume
% DNaJ2_average = DV*sum(DNapJDPHI2_CS,'all')/Volume
% delta_NUFIS1 = feedback.NUFIS1-input.NUFIS1;
% delta_ABS1 = feedback.ABS1 - input.ABS1;
% delta_NUFIS2 = feedback.NUFIS2-input.NUFIS2;
% delta_ABS2 = feedback.ABS2 - input.ABS2;
% delta_REM = feedback.REM - input.REM;
% Deltarho = DV * sum( (delta_NUFIS1-delta_ABS1-delta_REM).*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1) + delta_NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1)+delta_REM.*MOD2_adj(:,:,:,1).*MOD1(:,:,:,1) - delta_ABS2.*MOD2_adj(:,:,:,1).*MOD2(:,:,:,1),'all') ...
%     / (DV * sum(input.NUFIS1.*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1)+input.NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1),'all'));
% Deltarho_reduced = DV * sum(  delta_NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1) - delta_ABS2.*MOD2_adj(:,:,:,1).*MOD2(:,:,:,1),'all') ...
%     / (DV * sum(input.NUFIS1.*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1)+input.NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1),'all'));
% keff2_SIM = 1.00143;
% reac1_SIM = 0;
% reac2_SIM = (keff2_SIM-1)/keff2_SIM;
% x = DX*(1:sizex);
% y = DY*(1:sizey);
% z = DZ*(1:sizez);
% [X,Y] = meshgrid(x,y);
% mat_DSA(:,:) = DABS2PHI2DPHI2(:,:,13);
% vec_DSA(:) = DABS2PHI2DPHI2(17,17,:);
% figure()
% surf(X,Y,mat_DSA);
% view(2)
% colorbar
% zlabel('^{\Delta\Sigma_{a,2}\phi_{2}}/_{\Delta\phi_{1}}',"FontSize",20)
% xlabel("X (cm)")
% ylabel("Y (cm)")
% title("Absorption")
% figure()
% plot(vec_DSA);
% xlabel("Height (cm)")
% ylabel('^{\Delta\Sigma_{a,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
% title("Absorption")
% mat_DSF(:,:) = DNUFIS2PHI2DPHI2(:,:,13);
% vec_DSF(:) = DNUFIS2PHI2DPHI2(17,17,:);
% figure()
% surf(mat_DSF)
% view(2)
% colorbar
% zlabel('^{\Delta\Sigma_{f,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
% xlabel("X (cm)")
% ylabel("Y (cm)")
% title("Fission")
% figure()
% plot(vec_DSF);
% xlabel("Height (cm)")
% ylabel('^{\Delta\Sigma_{f,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
% title("Fission")
% mat_DNaJ(:,:) = DNapJDPHI2(:,:,13);
% vec_DNaJ(:) = DNapJDPHI2(17,17,:);
% figure()
% surf(mat_DNaJ)
% view(2)
% colorbar
% zlabel('^{\Delta(\nabla J_{2}\phi_{2})}/_{\Delta\phi_{2}}',"FontSize",20)
% xlabel("X (cm)")
% ylabel("Y (cm)")
% title("Leakage")
% figure()
% plot(vec_DNaJ);
% xlabel("Height (cm)")
% ylabel('^{\Delta(\nabla J_{2}\phi_{2})}/_{\Delta\phi_{2}}',"FontSize",20)
% title("Leakage")
%% clear variables
%clear NUFIS1 NUFIS2 n_iter n_restart ABS1 ABS2 FLX1 FLX2 EIG_MET D1 D2 MOD1_adj MOD2_adj REM RES_FLX RES_MOD RES_MOD_adj ...
 %   XS KN DX DY DZ conv_ERAM conv_POW lambda lambda_adj gammaI gammaX v1 v2
%% initialise parameters
PHID_F_PHI = zeros(1,M);
PHID_V_PHI = zeros(1,M);
PHID_PHI = zeros(1,M);
PHID_GAMMAI_PHI = zeros(M);
PHID_GAMMAX_PHI = zeros(M);
PHID_FB_PHI = zeros(M);
PHID_PHIUPPER_PHI=zeros(M);
PHID_PHILOWER_PHI=zeros(M);
PHID_X0_PHI=zeros(M);

temp_GAMMAX_PHI = zeros(sizex*2,sizey,sizez,M);
temp_GAMMAI_PHI = zeros(sizex*2,sizey,sizez,M);
temp_F_PHI = zeros(sizex*2,sizey,sizez,M);
temp_V_PHI = zeros(sizex*2,sizey,sizez,M);
%temp_PHI_eq_mat_PHI=zeros(sizex*2,sizey,sizez,M);
temp_FB_PHI= zeros(sizex*2,sizey,sizez,M);
% temp_FB_PHI_CS= zeros(sizex*2,sizey,sizez,M);
temp_PHIUPPER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_PHILOWER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_X0_PHI = zeros(sizex*2,sizey,sizez,M);
%temp_CR_PHI = zeros(sizex*2,sizey,sizez,M);
%% Calculate kets
for n = 1:M
    temp_F_PHI(:,:,:,n)=G2_inner_product(F(:,:,:),MOD(:,:,:,n),"matrix","vector"); % | F* Phi_n>
    %temp_FB_PHI(:,:,:,n) = G2_inner_product(MOD_eq_MAT,MOD(:,:,:,n),"matrix","vector"); %|Phi_0_mat * Phi_n>
    temp_FB_PHI(:,:,:,n) = G2_inner_product(DeltaXS,MOD(:,:,:,n),'matrix','vector'); % | FB *Phi_n>
    %temp_FB_PHI(:,:,:,n) =
    %[DeltaXS_CS(1,1)*MOD1(:,:,:,n)+DeltaXS_CS(1,2)*MOD2(:,:,:,n);DeltaXS_CS(2,1)*MOD1(:,:,:,n)+DeltaXS_CS(2,2)*MOD2(:,:,:,n)];
    %% | FB *Phi_n> % using Homogeneous feedback matrix
    temp_GAMMAX_PHI(:,:,:,n) = G2_inner_product(GAMMAX,MOD(:,:,:,n),"matrix","vector"); %  |1/k_0 *Gamma_X * Phi_n>
    temp_GAMMAI_PHI(:,:,:,n) = G2_inner_product(GAMMAI,MOD(:,:,:,n),"matrix","vector"); % | 1/k_0 Gamma_I * Phi_n>
    temp_PHIUPPER_PHI(:,:,:,n) = G2_inner_product(MOD_UPPER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \bar{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_PHILOWER_PHI(:,:,:,n) = G2_inner_product(MOD_LOWER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \tilde{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_X0_PHI(:,:,:,n) = G2_inner_product(Xe_UPPER,MOD(:,:,:,n),"matrix","vector"); % |\bar{X} * X0 * Phi_n >
end

%% Calculate bra-kets 

for m = 1:M
    PHID_F_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_F_PHI(:,:,:,m),"vector","vector"),'all'); %<Phi^dagger_m| F Phi_m>
    temp_V_PHI(:,:,:,m) = G2_inner_product(V_inv,MOD(:,:,:,m),"scalar_matrix","vector"); % |v^-1 Phi_m>
    PHID_V_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_V_PHI(:,:,:,m),"vector","vector"),"all"); %<Phi^dagger_m| v^-1 Phi_m>
    PHID_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),MOD(:,:,:,m),"vector","vector"),"all"); % <Phi^dagger_m|Phi_m>
    for n = 1:M
        PHID_GAMMAX_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAX_PHI(:,:,:,n),"vector","vector"),"all"); % <Phi^dagger_m | Gamma_X * Phi_n >
        PHID_GAMMAI_PHI(m,n) =  DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAI_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | Gamma_I * Phi_n > 
        PHID_FB_PHI(m,n) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_FB_PHI(:,:,:,n),"vector","vector"),"all");
        PHID_PHIUPPER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHIUPPER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_PHILOWER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHILOWER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \tilde{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_X0_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_X0_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} X_0 Phi_n >
    end
end


LAMBDA = PHID_V_PHI./ PHID_F_PHI; % <Phi^dagger_m |v^-1 Phi_n>/ <Phi^dagger_m |F Phi_m> 

%% clear temp variables
%clearvars temp*
%% Test feedback terms
% make dimensionless values
MOD1_max(:) = max(MOD1,[],[1,2,3]);
MOD2_max(:) = max(MOD2,[],[1,2,3]);
scaled_flux1_max = max(MOD_EQ_1_scaled,[],'all');
scaled_flux2_max = max(MOD_EQ_2_scaled,[],'all');
%% test magnitudes of terms
FB = zeros(M);
P1 = zeros(M);
X = zeros(M);
P2 = zeros(M);
for m= 1:M
    for n = 1:M
        FB(m,n) = 1./LAMBDA(m).*PHID_FB_PHI(m,n)./PHID_F_PHI(m);
        P1(m,n) = 1./LAMBDA(m).*PHID_PHILOWER_PHI(m,n)/PHID_F_PHI(m)*PHID_PHI(n)/PHID_F_PHI(n);
        X(m,n) = PHID_X0_PHI(m,n)/PHID_PHI(m);
        P2(m,n) = PHID_PHIUPPER_PHI(m,n)/PHID_PHI(m)*PHID_PHI(n)/PHID_F_PHI(n);
    end
end
FB_norm = FB/FB(1,1);
P1_norm = P1/P1(1,1);
X_norm = X/X(1,1);
P2_norm = P2/P2(1,1);

fprintf("Printing eig separation term \n")
1./LAMBDA'.*(1./keff-1./K_VALUE) 
fprintf("printing feedback term \n")
1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI
fprintf("printing diagonal of feedback term \n")
diag(1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI)
C11 = diag(1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI) + 1./LAMBDA'.*(1./keff-1./K_VALUE)
fprintf("printing flux xenon term \n")
C13 = diag(sigmaX.*P1);   
fprintf("Printing Iodine creation term \n")
C21 = gammaI/NU*PHID_F_PHI./PHID_PHI
fprintf("Printing first xenon absorption term \n")
C22 = -lambdaI;
C31 = diag(gammaX/NU *PHID_F_PHI./PHID_PHI-sigmaX*X)
fprintf("Printing second xenon absorption term \n")
C22 = lambdaI;
C33 = diag(-sigmaX.*P2)
save(RESULTS_REFINEMENT1_dir+"PARAMETERS_data.mat","PHID_GAMMAI_PHI","PHID_GAMMAX_PHI","PHID_PHI","PHID_PHI","PHID_F_PHI","PHID_FB_PHI","PHID_X0_PHI","PHID_V_PHI","PHID_PHIUPPER_PHI","PHID_PHILOWER_PHI","LAMBDA","M","keff","K_VALUE")
save(RESULTS_REFINEMENT1_dir+"COUPLING_data.mat","FB_norm","P1_norm","X_norm","P2_norm","C11","C13","C21","C22","C31","C33","temp_X0_PHI","temp_PHILOWER_PHI","temp_PHIUPPER_PHI","PHID_PHI","PHID_F_PHI","MOD","MOD_adj","LAMBDA")
%%
% matrix = 1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI/(1./LAMBDA(1).*PHID_FB_PHI(1,1)./PHID_F_PHI(1));
% % Create a heatmap
% figure;
% h = heatmap(matrix);
% % Apply the formatted labels to the heatmap
% h.CellLabelFormat = '%0.2f'; % Format for the heatmap cells
% h.ColorLimits = [min(matrix(:)), max(matrix(:))]; % Adjust color scaling to matrix range
% labels = {"Fundamental Mode", "First Axial Harmonic Mode", ...
%           "First Radial Harmonic Mode", "Second Radial Harmonic Mode", ...
%           "Second Axial Harmonic Mode"};
% % Add custom axis labels
% h.XDisplayLabels = labels; % Set column (X-axis) labels
% h.YDisplayLabels = labels; % Set row (Y-axis) labels
% % Add axis labels
% %h.XLabel = 'Column Index';
% %h.YLabel = 'Row Index';
% 
% % Customize color scheme
% colormap('parula'); % Change to any preferred colormap ('parula', 'hot', etc.)
% colorbar; % Add a colorbar for reference
% 
% % Set title
% %title('Heatmap of Matrix Values');
% 
% % Format the numbers to display 2 decimals
% % Generate custom labels for each cell
% dataLabels = arrayfun(@(x) sprintf('%.2f', x), matrix, 'UniformOutput', false);

%%
[time_2G,state_values_2G] = ode_Nsolve(RESULTS_REFINEMENT1_dir);
%clearvars -except time_2G state_values_2G C11 results_dir;
%investigation = "Feedback_coefficient_variations_10_24";
if ~exist(RESULTS_REFINEMENT1_dir,"dir")
    mkdir(RESULTS_REFINEMENT1_dir)
end
C11
save(RESULTS_REFINEMENT1_dir+"RESULTS_HET.mat","state_values_2G","time_2G")

end