function data = modifyData(data, opts)

    
    skip_xerom = opts.skip_xerom; skip_mscnpp = opts.skip_mscnpp; jump = opts.jump_xerom;
    
    data.DV = data.DX*data.DY*data.DZ;    
    data.sizex = size(data.xerom.MOD1,1);
    data.sizey = size(data.xerom.MOD1,2);
    data.sizez = size(data.xerom.MOD1,3);

    data.xerom.ntimesteps = size(data.xerom.time_2G,1);
    reduced_time_indices = 1:jump:data.xerom.ntimesteps;
    reduced_time_hours = data.xerom.time_2G(reduced_time_indices)/3600;
    reduced_time_hours = reduced_time_hours - reduced_time_hours(1);
    data.xerom.reduced_time_hours = reduced_time_hours;
    data.mscnpp.t_hours = data.mscnpp.t/3600;

    KFIS1 = (1/data.NU) * data.NUFIS1 .* data.KAPPA1;
    KFIS2 = (1/data.NU) * data.NUFIS2 .* data.KAPPA2;
    FIS1 = (1/data.NU) * data.NUFIS1 ;
    FIS2 = (1/data.NU) * data.NUFIS2;
    data.KFIS1 = KFIS1; data.KFIS2 = KFIS2; data.FIS1 = FIS1; data.FIS2 = FIS2;
    

    MOD = cat(1,data.xerom.MOD1,data.xerom.MOD2);
    MOD_EQ = abs(cat(1,data.xerom.MOD1(:,:,:,1), data.xerom.MOD2(:,:,:,1)));
    data.xerom.MOD = MOD; data.xerom.MOD_EQ = MOD_EQ;

    KFISINT = data.DV * sum(G2_inner_product([KFIS1,KFIS2],MOD_EQ,"vector","vector"),"all");
    PS = data.mscnpp.power / KFISINT;
    MOD_EQ_scaled = MOD_EQ * PS;
    data.FLX1_EQ_scaled = MOD_EQ_scaled(1:data.sizex,:,:);
    data.FLX2_EQ_scaled = MOD_EQ_scaled(data.sizex+1:2*data.sizex,:,:);

    data.mscnpp.DFLX1 = data.mscnpp.FLX1 - data.FLX1_EQ_scaled;
    data.mscnpp.DFLX2 = data.mscnpp.FLX2 - data.FLX2_EQ_scaled;
    data.mscnpp.DFLX = cat(1,data.mscnpp.DFLX1,data.mscnpp.DFLX2);
    data.mscnpp.a = computeModalAmplitudes(data.mscnpp.DFLX, MOD);

    neutron_mode_indices = ((1:data.M)-1)*3 + 1; % neutron entries in state vector
    iodine_mode_indices = ((1:data.M)-1)*3 + 2; % iodine entries in state vector
    xenon_mode_indices = ((1:data.M)-1)*3 + 3; % xenon entries in state vector
    %data.xerom.neutron_mode_indices = neutron_mode_indices;  data.xerom.iodine_mode_indices = iodine_mode_indices; data.xerom.xenon_mode_indices = xenon_mode_indices; 
    data.xerom.neutronModes = data.xerom.state_values_2G(reduced_time_indices, neutron_mode_indices);
    data.xerom.iodineModes = data.xerom.state_values_2G(reduced_time_indices, iodine_mode_indices);
    data.xerom.xenonModes = data.xerom.state_values_2G(reduced_time_indices, xenon_mode_indices);
end
