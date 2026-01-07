function computeModalSignals(MOD,a,InputFolder,OutputFolder)
    load(InputFolder+"XS_data.mat","NUFIS1","NUFIS2");
    load(OutputFolder+"RESULTS_HET.mat","state_values_2G","time_2G");
    load(OutputFolder+"PARAMETERS_data.mat","PHID_PHI","PHID_F_PHI","M");

    %ZEROS = zeros(size(NUFIS1));
    %skip = 100;
    %jump = 10;
    sizex = size(NUFIS1,1);
    sizey = size(NUFIS1,2);
    sizez = size(NUFIS1,3);
    %MOD = [MOD1;MOD2]; % vector of solutions to the forward problem
    %F = [NUFIS1, NUFIS2;ZEROS,ZEROS];
    %A = PHID_PHI./PHID_F_PHI;
    %% Create time and state variables
    nt = size(a,1);
    %reduced_indeces = skip:jump:n_timesteps_full;
    %t_reduced = time_2G(reduced_indeces);
    %state_values_flux = state_values_2G(reduced_indeces,1:3:end);
    %state_values_iodine = state_values_2G(reduced_indeces,2:3:end);
    %state_values_xenon = state_values_2G(reduced_indeces,3:3:end);
    %n_timesteps_reduced = size(t_reduced,1);


    %% Calculate fluxes and FP concentrations
    recreated_flux = zeros(sizex*2,sizey,sizez,nt); %Initialise Flux state vector
    %I = zeros(sizex,sizey,sizez,n_timesteps_reduced); % Initialise I-135 state vector
    %X = zeros(sizex,sizey,sizez,n_timesteps_reduced); % Initialise X-135 state vector
    state_flux_temp = zeros(sizex*2,sizey,sizez,M);
    %state_iodine_temp = zeros(sizex,sizey,sizez,M);
    %state_xenon_temp = zeros(sizex,sizey,sizez,M);
    %IX_vector = zeros(sizex,sizey,sizez,M);
    for t = 1:nt
        for m = 1:M
            %IX_vector(:,:,:,m) = A(m).*(F(1,1)*MOD(1:sizex,:,:,m)+F(1,2)*MOD(sizex+1:end,:,:,m));
            state_flux_temp(:,:,:,m)=a(t,m)*MOD(:,:,:,m);
            %state_iodine_temp(:,:,:,m) = state_values_iodine(t,m).*IX_vector(:,:,:,m);
            %state_xenon_temp(:,:,:,m) = state_values_xenon(t,m).*IX_vector(:,:,:,m);
        end
        recreated_flux(:,:,:,t) = sum(state_flux_temp,4);
        %I(:,:,:,t)= sum(state_iodine_temp,4);
        %X(:,:,:,t) = sum(state_xenon_temp,4);
    end
    FLX1(:,:,:,:) = recreated_flux(1:sizex,:,:,:);
    FLX2(:,:,:,:) = recreated_flux(sizex+1:end,:,:,:);
    save(OutputFolder+"Signals.mat", "FLX1","FLX2")
end