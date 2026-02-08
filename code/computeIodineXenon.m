function data = computeIodineXenon(data)
MOD1 = data.xerom.MOD1;
MOD2 = data.xerom.MOD2;
iodine_modes = data.xerom.iodineModes;
xenon_modes = data.xerom.xenonModes;
nt = length(data.xerom.reduced_time_hours);
M = data.M;
PHID_PHI = data.xerom.PHID_PHI;
PHID_F_PHI = data.xerom.PHID_F_PHI;
FIS1 = data.FIS1;
FIS2 = data.FIS2;
top_indices = data.top_indices;
bottom_indices = data.bottom_indices;

% Get dimensions
[NX, NY, NZ, ~] = size(FIS1);

% Calculate Xenon and Iodine concentrations
A1 = reshape(MOD1, [], M); % (NX*NY*NZ) x M
A2 = reshape(MOD2, [], M); % (NX*NY*NZ) x M
% Expand FIS to match mode dimensions
FIS1_expanded = repmat(FIS1, [1, 1, 1, M]); % NX x NY x NZ -> NX x NY x NZ x M
FIS2_expanded = repmat(FIS2, [1, 1, 1, M]); % NX x NY x NZ -> NX x NY x NZ x M
FIS1_reshaped = reshape(FIS1_expanded, [], M); % (NX*NY*NZ) x M
FIS2_reshaped = reshape(FIS2_expanded, [], M); % (NX*NY*NZ) x M

% Calculate Q_n scaling factors
QN = PHID_PHI ./ PHID_F_PHI; % 1 x M

% Calculate spatial part: Q_n * F(r) * Phi_n(r)
Spatial_part = A1 .* FIS1_reshaped + A2 .* FIS2_reshaped; % (NX*NY*NZ) x M
Scaled_spatial_part = Spatial_part .* QN; % (NX*NY*NZ) x M

% Initialize output arrays
delta_iodine_concentration = zeros(NX, NY, NZ, nt);
delta_xenon_concentration = zeros(NX, NY, NZ, nt);

%Stationary xenon and iodine
Xe_eq = data.xerom.X0;
I_eq = data.xerom.I0;
% Loop over time steps to compute concentrations
for t = 1:nt
    % Sum over modes: Σ Q_n * F(r) * I_n(t) * Phi_n(r)
    iodine_sum = Scaled_spatial_part * iodine_modes(t, :)'; % (NX*NY*NZ) x 1
    xenon_sum = Scaled_spatial_part * xenon_modes(t, :)'; % (NX*NY*NZ) x 1
    
    % Reshape back to 3D spatial grid
    delta_iodine_concentration(:, :, :, t) = reshape(iodine_sum, NX, NY, NZ);
    delta_xenon_concentration(:, :, :, t) = reshape(xenon_sum, NX, NY, NZ);
end
clearvars iodine_modes xenon_modes Scaled_spatial_part Scaled_spatial_part FIS1_reshaped FIS2_reshaped iodine_sum xenon_sum

iodine_concentration = delta_iodine_concentration+I_eq;
xenon_concentration = delta_xenon_concentration + Xe_eq;
iodine_mean(:) = mean(iodine_concentration,[1,2,3]);
iodine_sum(:) = sum(iodine_concentration,[1,2,3]);
iodine_top(:) = sum(iodine_concentration(:,:,top_indices,:),[1,2,3]);
iodine_bottom(:) = sum(iodine_concentration(:,:,bottom_indices,:),[1,2,3]);
iodine_top_norm(:) = iodine_top./iodine_sum;
iodine_bottom_norm(:) = iodine_bottom./iodine_sum;

xenon_mean(:) = mean(xenon_concentration,[1,2,3]);
xenon_sum(:) = sum(xenon_concentration,[1,2,3]);
xenon_top(:) = sum(xenon_concentration(:,:,top_indices,:),[1,2,3]);
xenon_bottom(:) = sum(xenon_concentration(:,:,bottom_indices,:),[1,2,3]);
xenon_top_norm(:) = xenon_top./xenon_sum;
xenon_bottom_norm(:) = xenon_bottom./xenon_sum;
%Calculate mscnpp mean and top and bottom
mscnpp_xenon = data.mscnpp.Xe135;
mscnpp_iodine = data.mscnpp.I135;

mscnpp_xenon_mean(:) = mean(mscnpp_xenon,[1,2,3]);
mscnpp_iodine_mean(:) = mean(mscnpp_iodine,[1,2,3]);
% Calculate mscnpp top and bottom normalized values
mscnpp_xenon_top(:) = sum(mscnpp_xenon(:,:,top_indices,:),[1,2,3]);
mscnpp_xenon_bottom(:) = sum(mscnpp_xenon(:,:,bottom_indices,:),[1,2,3]);
mscnpp_iodine_top(:) = sum(mscnpp_iodine(:,:,top_indices,:),[1,2,3]);
mscnpp_iodine_bottom(:) = sum(mscnpp_iodine(:,:,bottom_indices,:),[1,2,3]);
mscnpp_xenon_top_norm(:) = mscnpp_xenon_top ./ mscnpp_xenon_mean;
mscnpp_xenon_bottom_norm(:) = mscnpp_xenon_bottom ./ mscnpp_xenon_mean;
mscnpp_iodine_top_norm(:) = mscnpp_iodine_top ./ mscnpp_iodine_mean;
mscnpp_iodine_bottom_norm(:) = mscnpp_iodine_bottom ./ mscnpp_iodine_mean;

% Store mscnpp results in data structure
data.mscnpp.xenon_top_norm = mscnpp_xenon_top_norm;
data.mscnpp.xenon_bottom_norm = mscnpp_xenon_bottom_norm;
data.mscnpp.iodine_top_norm = mscnpp_iodine_top_norm;
data.mscnpp.iodine_bottom_norm = mscnpp_iodine_bottom_norm;

% Store results in data structure
data.xerom.iodine_concentration = iodine_concentration;
data.xerom.xenon_concentration = xenon_concentration;
data.xerom.mean_iodine_concentration = iodine_mean;
data.xerom.mean_xenon_concentration = xenon_mean;

data.xerom.iodine_top_norm = iodine_top_norm;
data.xerom.iodine_bottom_norm = iodine_bottom_norm;
data.xerom.iodine_top = iodine_top;
data.xerom.iodine_bottom = iodine_bottom;

data.xerom.xenon_top_norm = xenon_top_norm;
data.xerom.xenon_bottom_norm = xenon_bottom_norm;
data.xerom.xenon_top = xenon_top;
data.xerom.xenon_bottom = xenon_bottom;

data.mscnpp.xenon_mean = mscnpp_xenon_mean;
data.mscnpp.xenon_top = mscnpp_xenon_top;
data.mscnpp.xenon_bottom = mscnpp_xenon_bottom;
data.mscnpp.xenon_top_norm = mscnpp_xenon_top_norm;
data.mscnpp.xenon_bottom_norm = mscnpp_xenon_bottom_norm;
data.mscnpp.iodine_mean = mscnpp_iodine_mean;
data.mscnpp.iodine_top = mscnpp_iodine_top;
data.mscnpp.iodine_bottom = mscnpp_iodine_bottom;
data.mscnpp.iodine_top_norm = mscnpp_iodine_top_norm;
data.mscnpp.iodine_bottom_norm = mscnpp_iodine_bottom_norm;

end