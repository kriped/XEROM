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

% Get dimensions
[NX, NY, NZ, ~] = size(FIS1);

% Calculate Xenon and Iodine concentrations
A1 = reshape(MOD1, [], M); % (NX*NY*NZ) x M
A2 = reshape(MOD2, [], M); % (NX*NY*NZ) x M
% Expand FIS to match mode dimensions
FIS1_expanded = repmat(FIS1, [1, 1, 1, M]); % NX x NY x NZ -> NX x NY x NZ x M
FIS2_expanded = repmat(FIS2, [1, 1, 1, M]); % NX x NY x NZ -> NX x NY x NZ x M
FIS1_reshaped = reshape(FIS1_expanded, [], M); % (NX*NY*NZ) x M
FIS2_reshaped = reshape(FIS1_expanded, [], M); % (NX*NY*NZ) x M

% Calculate Q_n scaling factors
QN = PHID_PHI ./ PHID_F_PHI; % 1 x M

% Calculate spatial part: Q_n * F(r) * Phi_n(r)
Spatial_part = A1 .* FIS1_reshaped + A2 .* FIS2_reshaped; % (NX*NY*NZ) x M
Scaled_spatial_part = Spatial_part .* QN; % (NX*NY*NZ) x M

% Initialize output arrays
iodine_concentration = zeros(NX, NY, NZ, nt);
xenon_concentration = zeros(NX, NY, NZ, nt);

% Loop over time steps to compute concentrations
for t = 1:nt
    % Sum over modes: Σ Q_n * F(r) * I_n(t) * Phi_n(r)
    iodine_sum = Scaled_spatial_part * iodine_modes(t, :)'; % (NX*NY*NZ) x 1
    xenon_sum = Scaled_spatial_part * xenon_modes(t, :)'; % (NX*NY*NZ) x 1
    
    % Reshape back to 3D spatial grid
    iodine_concentration(:, :, :, t) = reshape(iodine_sum, NX, NY, NZ);
    xenon_concentration(:, :, :, t) = reshape(xenon_sum, NX, NY, NZ);
end

% Store results in data structure
data.xerom.iodine_concentration = iodine_concentration;
data.xerom.xenon_concentration = xenon_concentration;
end