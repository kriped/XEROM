function data = computePowerandAO(data)
% computePowerandAO - Computes power and area of operation (AO) for XEROM and MSCNPP.
%
% Syntax: data = computePowerandAO(data, opts)
%
% Inputs:
%   data - A structure containing the necessary data for calculations, including:
%       - scaling: scaling factors for conversion
%       - xerom: structure with fields DFLX1 and DFLX2 for XEROM calculations
%       - mscnpp: structure with fields DFLX1 and DFLX2 for MSCNPP calculations
%       - FLX1_EQ_scaled: scaled equilibrium flux for FLX1
%       - FLX2_EQ_scaled: scaled equilibrium flux for FLX2
%       - power: original power values
%       - KFIS1: first set of coefficients for flux calculations
%       - KFIS2: second set of coefficients for flux calculations
%       - DV: scaling factor for power calculations
%   opts - Options structure for additional parameters (not used in this function).
%
% Outputs:
%   data - Updated structure with calculated power and AO for both XEROM and MSCNPP.
%       - data.xerom.power: calculated power for XEROM
%       - data.xerom.power_top: power for the top half of the XEROM
%       - data.xerom.power_bottom: power for the bottom half of the XEROM
%       - data.xerom.AO: area of operation for XEROM
%       - data.mscnpp.power_signal: calculated power for MSCNPP
%       - data.mscnpp.power_top: power for the top half of the MSCNPP
%       - data.mscnpp.power_bottom: power for the bottom half of the MSCNPP
%       - data.mscnpp.AO: area of operation for MSCNPP
%load in data
power_scaling = data.scaling.xerom_to_mscnpp;
DFLX1_xerom =  power_scaling * data.xerom.DFLX1;
DFLX2_xerom =  power_scaling * data.xerom.DFLX2;
DFLX1_mscnpp =  data.mscnpp.DFLX1;
DFLX2_mscnpp =  data.mscnpp.DFLX2;
FLX1_eq = data.FLX1_EQ_scaled;
FLX2_eq = data.FLX2_EQ_scaled;
power_eq = data.mscnpp.power;
KFIS1 = data.KFIS1;
KFIS2 = data.KFIS2;
DV = data.DV;

    %% Calculate power   
    
    sizez =  size(DFLX1_xerom,3);
    top_indices = (sizez/2+1):sizez;
    bottom_indices = 1:(sizez/2);
    n_timesteps_xerom = size(DFLX1_xerom,4);
    n_timesteps_mscnpp = size(DFLX1_mscnpp,4); 
    KFIS1_expanded = repmat(KFIS1,[1,1,1,n_timesteps_xerom]);
    KFIS2_expanded = repmat(KFIS2,[1,1,1,n_timesteps_xerom]);
    %% Calculate Power for XEROM
    FLX1_xerom = DFLX1_xerom+FLX1_eq;
    FLX2_xerom = DFLX2_xerom+FLX2_eq;
    data.xerom.power(:) = DV*sum(KFIS1_expanded.*FLX1_xerom+KFIS2_expanded.*FLX2_xerom,[1,2,3]);
    power_eq_test(:) = DV*sum(KFIS1.*FLX1_eq+KFIS2.*FLX2_eq,'all');
    fprintf("Original power = %0.4e, Calculated power = %0.4e \n",power_eq,power_eq_test)
    
    %% Calculate AO for XEROM
    data.xerom.power_top(:) = DV *sum(KFIS1_expanded(:,:,top_indices,:).*FLX1_xerom(:,:,top_indices,:)+KFIS2_expanded(:,:,top_indices,:).*FLX2_xerom(:,:,top_indices,:),[1,2,3]);
    data.xerom.power_bottom(:) = DV *sum(KFIS1_expanded(:,:,bottom_indices,:).*FLX1_xerom(:,:,bottom_indices,:)+KFIS2_expanded(:,:,bottom_indices,:).*FLX2_xerom(:,:,bottom_indices,:),[1,2,3]);
    data.xerom.AO(:) = (data.xerom.power_top-data.xerom.power_bottom)./(data.xerom.power_top+data.xerom.power_bottom);
    %% Reexpand KFIS
    clear KFIS1_expanded KFIS2_expanded
    KFIS1_expanded = repmat(KFIS1,[1,1,1,n_timesteps_mscnpp]);
    KFIS2_expanded = repmat(KFIS2,[1,1,1,n_timesteps_mscnpp]);
    %% Calculate Power for MSCNPP
    FLX1_mscnpp = FLX1_eq + DFLX1_mscnpp;
    FLX2_mscnpp = FLX2_eq + DFLX2_mscnpp;
    data.mscnpp.power_signal(:) = DV * sum(KFIS1_expanded.*FLX1_mscnpp+KFIS2_expanded.*FLX2_mscnpp,[1,2,3]);
    %% Calculate AO for MSCNPP
    data.mscnpp.power_top(:) = DV * sum(KFIS1_expanded(:,:,top_indices,:).*FLX1_mscnpp(:,:,top_indices,:)+KFIS2_expanded(:,:,top_indices,:).*FLX2_mscnpp(:,:,top_indices,:),[1,2,3]);
    data.mscnpp.power_bottom(:) = DV * sum(KFIS1_expanded(:,:,bottom_indices,:).*FLX1_mscnpp(:,:,bottom_indices,:)+KFIS2_expanded(:,:,bottom_indices,:).*FLX2_mscnpp(:,:,bottom_indices,:),[1,2,3]);
    data.mscnpp.AO(:) = (data.mscnpp.power_top - data.mscnpp.power_bottom) ./ (data.mscnpp.power_top + data.mscnpp.power_bottom);
end