function data = modifyData(data, opts)

    
    jump = opts.jump_xerom;
    
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
    
    data.top_indices = (data.sizez/2+1):data.sizez;
    data.bottom_indices = 1:(data.sizez/2);
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
    data.mscnpp.mean_xenon = mean(data.mscnpp.Xe135,[1, 2, 3]);
    tot_xenon(:) = sum(data.mscnpp.Xe135,[1, 2, 3]);
    upper_xenon(:) = sum(data.mscnpp.Xe135(:,:,data.top_indices,:),[1, 2, 3]);
    lower_xenon(:) = sum(data.mscnpp.Xe135(:,:,data.bottom_indices,:),[1, 2, 3]);
    data.mscnpp.upper_xenon_norm = upper_xenon./tot_xenon;
    data.mscnpp.lower_xenon_norm = lower_xenon./tot_xenon;
    data.mscnpp.mean_iodine = mean(data.mscnpp.I135,[1, 2, 3]);
    tot_iodine = sum(data.mscnpp.I135,[1, 2, 3]);
    upper_iodine(:) = sum(data.mscnpp.I135(:,:,data.top_indices,:),[1, 2, 3]);
    lower_iodine(:) = sum(data.mscnpp.I135(:,:,data.bottom_indices,:),[1, 2, 3]);
    data.mscnpp.upper_iodine_norm = upper_iodine./tot_iodine;
    data.mscnpp.lower_iodine_norm = lower_iodine./tot_iodine;

    neutron_mode_indices = ((1:data.M)-1)*3 + 1; % neutron entries in state vector
    iodine_mode_indices = ((1:data.M)-1)*3 + 2; % iodine entries in state vector
    xenon_mode_indices = ((1:data.M)-1)*3 + 3; % xenon entries in state vector
    %data.xerom.neutron_mode_indices = neutron_mode_indices;  data.xerom.iodine_mode_indices = iodine_mode_indices; data.xerom.xenon_mode_indices = xenon_mode_indices; 
    data.xerom.neutronModes = data.xerom.state_values_2G(reduced_time_indices, neutron_mode_indices);
    data.xerom.iodineModes = data.xerom.state_values_2G(reduced_time_indices, iodine_mode_indices);
    data.xerom.xenonModes = data.xerom.state_values_2G(reduced_time_indices, xenon_mode_indices);
end
