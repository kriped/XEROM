function data = computeDeltaFlux(data)

    sizex = data.sizex;
    sizey = data.sizey;
    sizez = data.sizez;
    nt = length(data.xerom.reduced_time_hours);
    MOD = data.xerom.MOD;
    M = data.M;
    neutronModes = data.xerom.neutronModes;
    %Test nt is same as length of neutronModes
    if ~size(neutronModes,1) == nt
        error("nt is not equal to number of timesteps in neutronModes")
    end

    %% Calculate fluxes and FP concentrations

    A = reshape(MOD, [], M);
    R = A * neutronModes.'; 
    recreated_flux = reshape(R, [], sizey, sizez, nt);
    DFLX1(:,:,:,:) = recreated_flux(1:sizex,:,:,:);
    DFLX2(:,:,:,:) = recreated_flux(sizex+1:end,:,:,:);
    data.xerom.DFLX1 = DFLX1;
    data.xerom.DFLX2 = DFLX2;
    if ~size(data.xerom.DFLX1) == [sizex, sizey,sizez,nt]
        error("Size of FLX1 is not correct!")
    end
    if ~size(data.xerom.DFLX2) == [sizex, sizey,sizez,nt]
        error("Size of FLX2 is not correct!")
    end
    
end