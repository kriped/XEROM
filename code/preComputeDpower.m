function data = preComputeDpower(data)
%   load in data
    DFLX1_xerom = data.xerom.DFLX1;
    DFLX2_xerom = data.xerom.DFLX2;
    DFLX1_mscnpp = data.mscnpp.DFLX1;
    DFLX2_mscnpp = data.mscnpp.DFLX2;
    
    KFIS1 = data.KFIS1;
    KFIS2 = data.KFIS2;
    DV = data.DV;
    if size(DFLX1_xerom,1) ~= size(DFLX1_xerom,2)
        warning("X and Y dimensions of the fast flux are not equal! Check if the format is correct! \n")
    end
    if size(DFLX1_mscnpp,1) ~= size(DFLX1_mscnpp,2)
        warning("X and Y dimensions of the fast flux are not equal! Check if the format is correct! \n")
    end
    [~,max_dim]=max(size(DFLX1_xerom)); 
    if max_dim ~= 4
        warning("Time dimension is not the largest dimension of the XEROM dataset. Check if the format is correct! \n");
    end
    [~,max_dim]=max(size(DFLX1_mscnpp)); 
    if max_dim ~= 4
        warning("Time dimension is not the largest dimension of the MScNPP dataset. Check if the format is correct! \n");
    end
    n_timesteps_xerom = size(DFLX1_xerom,4);
    n_timesteps_mscnpp = size(DFLX1_mscnpp,4);
    KFIS1_expanded = repmat(KFIS1,[1,1,1,n_timesteps_xerom]);
    KFIS2_expanded = repmat(KFIS2,[1,1,1,n_timesteps_xerom]);
    %% Calculate delta Power for XEROM
    data.xerom.Dpower(:) = DV*sum(KFIS1_expanded.*DFLX1_xerom+KFIS2_expanded.*DFLX2_xerom,[1,2,3]);
    %% Reexpand KFIS
    clear KFIS1_expanded KFIS2_expanded
    KFIS1_expanded = repmat(KFIS1,[1,1,1,n_timesteps_mscnpp]);
    KFIS2_expanded = repmat(KFIS2,[1,1,1,n_timesteps_mscnpp]);
    %% Calculate delta Power for MSCNPP
    data.mscnpp.Dpower(:) = DV * sum(KFIS1_expanded.*DFLX1_mscnpp+KFIS2_expanded.*DFLX2_mscnpp,[1,2,3]);
end