function Results_processing(CASE,BURNUP)
INPUT_REFINEMENT0_dir = sprintf("../input/%s/%s/Refinement0/",CASE,BURNUP);
INPUT_REFINEMENT1_dir = sprintf("../input/%s/%s/Refinement1/",CASE,BURNUP);
RESULTS_REFINEMENT0_dir = sprintf("../results/%s/%s/Refinement0/",CASE,BURNUP);
RESULTS_REFINEMENT1_dir = sprintf("../results/%s/%s/Refinement1/",CASE,BURNUP);

load(INPUT_REFINEMENT1_dir+"GEOM_data.mat")
load(INOUT_REFINEMENT)
load(RESULTS_REFINEMENT1_dir+"Results_HET");
load(RESULTS_REFINEMENT1_dir"PARAMETERS_data.mat")

load("output\dataFile.mat", "MOD", "KFIS1", "KFIS2", "DV","DX","DY", "DZ","PHID_PHI", "F", "DV","MOD1","MOD2","h","PHID_F_PHI","PHID_FB_PHI","sizex","sizey","sizez","M");
load("reference\reference_time_series.mat");
skip = 60;
jump = 5;
x_HET = HET_data.DX:HET_data.DX:HET_data.DX*HET_data.sizex;
y_HET = HET_data.DY:HET_data.DY:HET_data.DY*HET_data.sizey;
z_HET = HET_data.DZ:HET_data.DZ:HET_data.DZ*HET_data.sizez;
M = HET_data.M;
sizex = HET_data.sizex;
sizey = HET_data.sizey;
sizez = HET_data.sizez;
sizex_ref = size(FLX1_ref,1);
sizey_ref = size(FLX1_ref,2);
sizez_ref = size(FLX1_ref,3);
%% Create time variables
number_of_timesteps_HET = size(HET_results.state_values_2G,1);
number_of_timesteps_ref = size(FLX1_ref,4);
time_ref = linspace(0,75,number_of_timesteps_ref);
number_of_equations = size(HET_results.state_values_2G,2);
reduced_time_HET = HET_results.time_2G(skip:jump:number_of_timesteps_HET);
reduced_time_hours_HET = reduced_time_HET/3600;
reduced_number_of_timesteps_HET = size(reduced_time_HET,1);
%state_values_flux_reshaped = state_values_all_reshaped(:,:,:,:,1:3:end);
%state_values_iodine_reshaped = state_values_all_reshaped(:,:,:,:,2:3:end);
%state_values_xenon_reshaped = state_values_all_reshaped(:,:,:,:,3:3:end);
%% Reduce size of time series
%number_of_timesteps_HOM = size(HOM_results.state_values_2G,1);

A_HET = HET_data.PHID_PHI./HET_data.PHID_F_PHI;

%state_values_all_reshaped_HOM = reshape(HOM_results.state_values_2G,[number_of_timesteps_HOM,1,1,1,number_of_equations]);
%state_values_all_reshaped_HET = reshape(HET_results.state_values_2G,[number_of_timesteps_HET,1,1,1,number_of_equations]);
state_values_flux_HET = HET_results.state_values_2G(skip:jump:number_of_timesteps_HET,1:3:end);
state_values_iodine_HET = HET_results.state_values_2G(skip:jump:number_of_timesteps_HET,2:3:end);
state_values_xenon_HET = HET_results.state_values_2G(skip:jump:number_of_timesteps_HET,3:3:end);

%% Calculate the reference power and delta axial offset
% Calculate power from flux 
KFIS1_expanded_ref = reshape(repmat(KF1_ref,[1,1,1,number_of_timesteps_ref]),[sizex_ref,sizey_ref,sizez_ref,number_of_timesteps_ref]);
KFIS2_expanded_ref = reshape(repmat(KF2_ref,[1,1,1,number_of_timesteps_ref]),[sizex_ref,sizey_ref,sizez_ref,number_of_timesteps_ref]);

FLX1_EQ_ref_expanded = reshape(repmat(FLX1_EQ_ref,[1,1,1,number_of_timesteps_ref]),[sizex_ref,sizey_ref,sizez_ref,number_of_timesteps_ref]);
FLX2_EQ_ref_expanded = reshape(repmat(FLX2_EQ_ref,[1,1,1,number_of_timesteps_ref]),[sizex_ref,sizey_ref,sizez_ref,number_of_timesteps_ref]);

Delta_FLX1_ref = FLX1_ref - FLX1_EQ_ref_expanded;
Delta_FLX2_ref = FLX2_ref - FLX2_EQ_ref_expanded;
Power_EQ = DV_ref*sum(KF1_ref.*FLX1_EQ_ref+KF2_ref.*FLX2_EQ_ref,"all");
power_ref_test(:) = DV_ref * sum(KFIS1_expanded_ref.*FLX1_ref+KFIS2_expanded_ref.*FLX2_ref,[1,2,3]);
delta_Power_ref(:) = DV_ref * sum(KFIS1_expanded_ref.*Delta_FLX1_ref+KFIS2_expanded_ref.*Delta_FLX2_ref,[1,2,3]);
delta_Power_ref_max = max(delta_Power_ref);
delta_Power_ref_rescaled = delta_Power_ref./delta_Power_ref_max;

% Calculate the axial offset
Power_top_ref(:) = DV_ref*sum(KFIS1_expanded_ref(:,:,sizez_ref/2+1:sizez_ref,:).*Delta_FLX1_ref(:,:,sizez_ref/2+1:sizez_ref,:)+KFIS2_expanded_ref(:,:,sizez_ref/2+1:sizez_ref,:).*Delta_FLX2_ref(:,:,sizez_ref/2+1:sizez_ref,:),[1,2,3]);
Power_bottom_ref(:) = DV_ref*sum(KFIS1_expanded_ref(:,:,1:sizez_ref/2,:).*Delta_FLX1_ref(:,:,1:sizez_ref/2,:)+KFIS2_expanded_ref(:,:,1:sizez_ref/2,:).*Delta_FLX2_ref(:,:,1:sizez_ref/2,:),[1,2,3]);
Power_top_ref_rescaled = Power_top_ref./Power_EQ;
Power_bottom_ref_rescaled = Power_bottom_ref./Power_EQ;
delta_AO_ref = Power_top_ref_rescaled-Power_bottom_ref_rescaled; 
%%
% figure
% plot(time_ref,delta_Power_ref_rescaled)
% hold on
% Calculate AO from flux
% [peaks, loc_peak] = findpeaks(delta_Power_ref_rescaled,time_ref,"MinPeakDistance",15);
% peaks = peaks(peaks>=0);
% loc_peak = loc_peak(peaks>=0);
% Plot the identified peaks
% plot(loc_peak, peaks, 'ro', 'MarkerFaceColor', 'r');
% peak_intervals = diff(loc_peak); % Time intervals between consecutive peaks
% avg_period = mean(peak_intervals); % Average period
%     ax = gca;
%     x_limits = ax.XLim;
%     y_limits = ax.YLim;
% 
%     annotation_pos_x = x_limits(1) + 0.7 * (x_limits(2) - x_limits(1));
%     annotation_pos_y = y_limits(1) + 0.1 * (y_limits(2) - y_limits(1));
%     annotationText = ['Period: ', num2str(avg_period, '%.1f'), ' h'];
%     Add annotation at a specific position (adjust coordinates for desired location)
%     text(annotation_pos_x, ...
%          annotation_pos_y, ...
%          annotationText, 'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');
Signal_processing(time_ref,delta_Power_ref_rescaled,"Power deviation")
hold off
figure()
plot(time_ref,Power_top_ref_rescaled)
hold on
plot(time_ref,Power_bottom_ref_rescaled)
hold off
legend("Top", "Bottom")
title("SIMULATE 5 power top vs bottom")
Signal_processing(time_ref,delta_AO_ref,"Axial offset deviation")

%% Calculate fluxes and FP concentrations
%Flux_time_dep = sum(bsxfun(@times,state_values_flux_reshaped,MOD),5);
recreated_flux_HET = zeros(length(skip:jump:number_of_timesteps_HET),HET_data.sizex*2,HET_data.sizey,HET_data.sizez);
delta_iodine_HET = zeros(length(skip:jump:number_of_timesteps_HET),HET_data.sizex,HET_data.sizey,HET_data.sizez);
delta_xenon_HET = zeros(length(skip:jump:number_of_timesteps_HET),HET_data.sizex,HET_data.sizey,HET_data.sizez);
state_flux_temp_HET = zeros(HET_data.sizex*2,HET_data.sizey,HET_data.sizez,HET_data.M);
state_iodine_temp_HET = zeros(HET_data.sizex,HET_data.sizey,HET_data.sizez,HET_data.M);
state_xenon_temp_HET = zeros(HET_data.sizex,HET_data.sizey,HET_data.sizez,HET_data.M);
IX_vector_HET = zeros(HET_data.sizex,HET_data.sizey,HET_data.sizez,HET_data.M);
for t = 1:reduced_number_of_timesteps_HET
    for m = 1:HET_data.M
        IX_vector_HET(:,:,:,m) = A_HET(m).*(HET_data.F(1,1)*HET_data.MOD(1:HET_data.sizex,:,:,m)+HET_data.F(1,2)*HET_data.MOD(HET_data.sizex+1:end,:,:,m));
        state_flux_temp_HET(:,:,:,m)=state_values_flux_HET(t,m)*HET_data.MOD(:,:,:,m);
        state_iodine_temp_HET(:,:,:,m) = state_values_iodine_HET(t,m).*IX_vector_HET(:,:,:,m);
        state_xenon_temp_HET(:,:,:,m) = state_values_xenon_HET(t,m).*IX_vector_HET(:,:,:,m);
    end
    recreated_flux_HET(t,:,:,:) = sum(state_flux_temp_HET,4);
    delta_iodine_HET(t,:,:,:)=sum(state_iodine_temp_HET,4);
    delta_xenon_HET(t,:,:,:) = sum(state_xenon_temp_HET,4);
end
clear state_flux_temp* state_xenon_temp* state_iodine_temp* IX_vector* time_int state_values_flux_HET state_values_iodine_HET state_values_xenon_HET 
Flux1_HET(:,:,:,:) = recreated_flux_HET(:,1:HET_data.sizex,:,:);
Flux2_HET(:,:,:,:) = recreated_flux_HET(:,HET_data.sizex+1:end,:,:);
save(sprintf("%s/Results.mat",folder_name), "delta_iodine_HET","Flux2_HET","delta_xenon_HET","Flux1_HET");

%% Calculate the power and AO
% Calculate the power
KFIS1_expanded_HET = reshape(repmat(HET_data.KFIS1,[reduced_number_of_timesteps_HET,1,1,1]),[reduced_number_of_timesteps_HET,HET_data.sizex,HET_data.sizey,HET_data.sizez]);
KFIS2_expanded_HET = reshape(repmat(HET_data.KFIS2,[reduced_number_of_timesteps_HET,1,1,1]),[reduced_number_of_timesteps_HET,HET_data.sizex,HET_data.sizey,HET_data.sizez]);

delta_power_HET(:) = HET_data.DV*sum(KFIS1_expanded_HET.*Flux1_HET+KFIS2_expanded_HET.*Flux2_HET,[2,3,4]);
%calculate the maximum values of the parameters

delta_power_max_HET = max(delta_power_HET);
%delta_iodine_max_HET = max(delta_iodine_HET,[],'all');
%delta_xenon_max_HET = max(delta_xenon_HET,[],'all');

%rescale the parameters by the maximum value

delta_power_rescaled_HET = delta_power_HET/delta_power_max_HET;
%delta_iodine_rescaled_HET = delta_iodine_HET/delta_iodine_max_HET;
%delta_xenon_rescaled_HET = delta_xenon_HET/delta_xenon_max_HET;

% delta_power_right_HET = zeros(reduced_number_of_timesteps_HET,1);
% delta_power_left_HET = zeros(reduced_number_of_timesteps_HET,1);
% for t = 1:length(reduced_time_hours_HET)
%     for k = 1:sizez 
%         flipped_slice_HET = fliplr(squeeze(Temp_power_matrix_HET(t,:,:,k)));
%         delta_power_right_HET(t) = delta_power_right_HET(t)+ HET_data.DV*sum(tril(flipped_slice_HET),'all');    
%         delta_power_left_HET(t) = delta_power_left_HET(t)+ HET_data.DV*sum(triu(flipped_slice_HET),'all');
%     end
% end
%Calculate the radial power values homogenous model
delta_power_top_HET(:) = HET_data.DV *sum(KFIS1_expanded_HET(:,:,:,HET_data.sizez/2+1:end).*Flux1_HET(:,:,:,HET_data.sizez/2+1:end)+KFIS2_expanded_HET(:,:,:,HET_data.sizez/2+1:end).*Flux2_HET(:,:,:,HET_data.sizez/2+1:end),[2,3,4]);
delta_power_top_rescaled_HET = delta_power_top_HET/delta_power_max_HET;
delta_power_bottom_HET(:) = HET_data.DV *sum(KFIS1_expanded_HET(:,:,:,1:HET_data.sizez/2).*Flux1_HET(:,:,:,1:HET_data.sizez/2)+KFIS2_expanded_HET(:,:,:,1:HET_data.sizez/2).*Flux2_HET(:,:,:,1:HET_data.sizez/2),[2,3,4]);
delta_power_bottom_rescaled_HET = delta_power_bottom_HET/delta_power_max_HET;
delta_AO_HET(:) = (delta_power_top_HET-delta_power_bottom_HET);
vecs1_HET(:,:) = Flux1_HET(:,34,34,:);
vecs2_HET(:,:) = Flux2_HET(:,34,34,:);
ymax_fast_HET = max(vecs1_HET,[],'all');
ymax_thermal_HET = max(vecs2_HET,[],'all');
clear Temp_power_matrix_HET KFIS2_expanded_HET KFIS1_expanded_HET

%% Plot modes
time_int = figure('Position', get(0, 'Screensize'));
plot(HET_results.time_2G(skip:end)/3600,HET_results.state_values_2G(skip:end,(1:8-1)*3+1),"LineWidth",2,"LineStyle","-")

grid on
xlim([0 150])    
%ylim([-4E10 6E10])
ylabel("Amplitude [cm^{-2}s^{-1}]",'Fontsize', 22)
xlabel("Time [h]",'Fontsize', 22)
%legend("Fundamental mode","First axial harmonic", "First radial harmonic","Second radial harmonic", "Second axial harmonic","Fontsize",22,"location","best")
%legend("First axial harmonic", ...
 %       "First radial harmonic", "Second radial harmonic","mode 5","mode 6","Second axial harmonic","Mode 8", "mode 9", "Mode 10");

legend("First axial harmonic", ...
        "First radial harmonic", "Second radial harmonic","mode 5","mode 6","Mode 8", "mode 9", "Mode 10");
ax2 = gca;
ax2.FontSize = 22;

%saveas(time_int,"Results/time_int.png");
%saveas(time_int,"Results/time_int.fig");
saveas(time_int,sprintf("%s/time_int_HET.svg",folder_name));
saveas(time_int,sprintf("%s/time_int_HET.fig",folder_name));

%% Print all the modes sequentially

% Set up the figure size
time_int = figure('Position', get(0, 'Screensize'));
%max_values_HOM = max(abs(HOM_results.state_values_2G(skip:end,:)),[],1);
max_values_HET = max(abs(HET_results.state_values_2G(skip:end,:)),[],1);
% Define line styles and colors for each mode
%lineStyles = ['-', '-', '-', '-', '-'];  % Line styles
%colors = lines(5);  % Get 5 distinct colors
%legendEntries = ["Fundamental mode", "First axial harmonic", ...
%        "First radial harmonic", "Second radial harmonic", "Second axial harmonic"];
%legendEntries=["First axial harmonic"];
%legendEntries = ["First axial harmonic", ...
%        "First radial harmonic", "Second radial harmonic","mode 5","mode 6", "Second axial harmonic","Mode 8", "mode 9", "Mode 10"];
legendEntries = ["First axial harmonic", ...
        "First radial harmonic", "Second radial harmonic","mode 5","mode 6","Mode 8", "mode 9", "Mode 10"];

% Loop through modes and generate all plots
for i = 1:M
    % Plot the current mode with full opacity
    %plot(HOM_results.time_2G(skip:end)/3600, HOM_results.state_values_2G(skip:end, (i-1)*3+1)/max_values_HOM((i-1)*3+1), ...
    %    'LineWidth', 2, 'LineStyle', lineStyles{i});  % Full opacity for current mode
    plot(HET_results.time_2G(skip:end)/3600, HET_results.state_values_2G(skip:end, (i-1)*3+1)/max_values_HET((i-1)*3+1), ...
        'LineWidth', 2);  % Full opacity for current mode    
        % Customize the plot
    grid on;
    xlim([0 70]);
    %ylim([-4E10 6E10]);
    ylabel('Amplitude [cm^{-2}s^{-1}]', 'FontSize', 14);
    xlabel('Time [h]', 'FontSize', 14);

    % Create a legend for modes
    
    legend(["2G HET "+legendEntries(i)], 'FontSize', 14, 'Location', 'best');

    % Set font size for axes
    ax2 = gca;
    ax2.FontSize = 14;

    % Save the figure with a suitable file name
    % saveas(time_int, sprintf('Results/mode_plot_%d.png', i));
    %saveas(time_int,sprintf('Results_comparison_10_24/mode_plot_%d_HOM.png', i'))
    %saveas(time_int,sprintf("%s/mode_plot_%d_HOM.png",investigation,folder_name2,i))
    saveas(time_int,sprintf("%s/mode_plot_%d_HET.svg",folder_name,i))
  
    % Clear the figure for the next iteration
    hold off;
end


%% Plot power and AO
delta_power_figure = figure();
plot(reduced_time_hours_HET,delta_power_rescaled_HET,"LineWidth",2);

ylim([-1,1])
%yticks(-0.5:0.1:1)
xlim([0,150])
grid on
xlabel("Time (h)","FontSize",12)
ylabel("Power deviation (AU)","FontSize",12)
legend("Heterogeneous model")
saveas(delta_power_figure,sprintf("%s/delta_power.svg",folder_name))

%saveas(delta_power_figure,"../Mode_variation_10_24/diag_05/delta_power_HOM.png")
delta_AO_figure = figure();
target_time = 40;
[~, target_index] = min(abs(reduced_time_hours_HET - target_time));
%plot(reduced_time_hours_HET(target_index:end)-target_time,delta_AO_HET(target_index:end)*100,"LineWidth",2);

grid on
xlim([0,150-target_time])
ylim([-40,40])
yticks(-30:5:40)
xlabel("Time (h)","FontSize",12)
ylabel("Axial offset deviation (%)","FontSize",12)
% legend("Heterogeneous model")
Signal_processing(reduced_time_hours_HET(target_index:end)-target_time,delta_AO_HET(target_index:end)*100,"Axial Offset Deviation HET")
saveas(delta_AO_figure,sprintf("%s/delta_AO.svg",folder_name))
%saveas(delta_AO_figure,"Results_comparison_10_24/delta_AO.png")

power_TB_plot_HET = figure();
plot(reduced_time_hours_HET,delta_power_top_rescaled_HET,"LineWidth",2)
hold on 
plot(reduced_time_hours_HET,delta_power_bottom_rescaled_HET,"LineWidth",2)
% ylim([-0.4,0.7])
% yticks([-0.4:0.1:0.7])
xlim([0,150])
grid on
xlabel("Time (h)","FontSize",12)
ylabel("Power deviation (AU)","FontSize",12)
legend("Top","Bottom")
%saveas(power_top_plot_HET,"Results_comparison_10_24/delta_power_TB_HET.png")
saveas(power_TB_plot_HET,sprintf("%s/delta_power_TB_HET.svg",folder_name))

% power_LR_plot_HET = figure();
% plot(reduced_time_hours_HET,delta_power_right_rescaled_HET,"LineWidth",2)
% hold on 
% plot(reduced_time_hours_HET,delta_power_left_rescaled_HET,"LineWidth",2)
% % ylim([-0.4,0.7])
% %yticks([-0.4:0.1:0.7])
% xlim([0,70])
% grid on
% xlabel("Time (h)","FontSize",12)
% ylabel("Power deviation (AU)","FontSize",12)
% legend("Right","Left")
% %saveas(power_top_plot_HET,"Results_comparison_10_24/delta_power_TB_HET.png")
% saveas(power_LR_plot_HET,sprintf("%s/%s/delta_power_LR_HET.svg",investigation,folder_name))

% delta_RO_plot = figure();
% plot(reduced_time_hours_HOM,delta_RO_HOM*100,"LineWidth",2)
% hold on
% plot(reduced_time_hours_HET,delta_RO_HET*100,"LineWidth",2)
% xlim([0,70])
% grid on
% xlabel("Time (h)","FontSize",12)
% ylabel("Radial offset deviation (%)","FontSize",12)
% legend("Homogeneous model","Heteronegous model")
% saveas(delta_RO_plot,sprintf("%s/%s/delta_RO.svg",investigation,folder_name))
%%
% f_fast = figure("Visible","off");
% f_thermal = figure("Visible","off");
% filename_fast = 'Results/fast_xenon_osicllation_HOM.gif';
% filename_thermal = 'Results/thermal_xenon_oscillation_HOM.gif';
% for frame = 60:20:length(vecs1)
%         %draw fast neutron flux
%         f_fast;
%         plot(z,vecs1(frame,:),"LineWidth",2)
%         ylim([0,ymax_fast])
%         xlabel("Height (cm)","FontSize",14)
%         ylabel("Fast neutron flux (cm^{-2}s^{-1})","FontSize",14)
%         % Add time stamp as annotation
%         current_time = time_2G_hours(frame);  % Get the current time
%         annotation_text = sprintf('Time: %.1f hours', current_time);  % Format the time string
%         text(1, ymax_fast * 0.9, annotation_text, 'FontSize', 12, 'FontWeight', 'bold')  % Add text at the top
% 
%         drawnow
% 
% 
%         % Capture the current figure as an image
%         frameImage = getframe(gcf);
%         im = frame2im(frameImage);  % Convert to image data
%         [imind, cm] = rgb2ind(im, 256);  % Index image for GIF format
% 
%         % Write to the GIF file
%         if frame == 1
%             imwrite(imind, cm, filename_fast, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
%         else
%             imwrite(imind, cm, filename_fast, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
%         end
%         % Draw thermal neutron flux
%         f_thermal;
%         plot(z,vecs2(frame,:),"LineWidth",2)
%         ylim([0,ymax_thermal])
%         xlabel("Height (cm)",FontSize=14)
%         ylabel("Thermal neutron flux (cm^{-2}s^{-1})","FontSize",14)
%         % Add time stamp as annotation
%         current_time = time_2G_hours(frame);  % Get the current time
%         annotation_text = sprintf('Time: %.1f hours', current_time);  % Format the time string
%         text(1, ymax_thermal * 0.9, annotation_text, 'FontSize', 12, 'FontWeight', 'bold')  % Add text at the top
% 
%         drawnow
% 
%         % Capture the current figure as an image
%         frameImage = getframe(gcf);
%         im = frame2im(frameImage);  % Convert to image data
%         [imind, cm] = rgb2ind(im, 256);  % Index image for GIF format
% 
%         % Write to the GIF file
%         if frame == 1
%             imwrite(imind, cm, filename_thermal, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
%         else
%             imwrite(imind, cm, filename_thermal, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
%         end
% end
% 
% %%
% delta_t = zeros(ceil(length(time_2G_hours)/50),1);
% for t = 2:50:length(time_2G_hours)
%     delta_t(t) = time_2G_hours(t)-time_2G_hours(t-1);
% end
% mean(delta_t)
%%
% delta_power_figure = figure();
% plot(reduced_time_hours_HET,delta_power_rescaled_HET,"LineWidth",2);
% ylim([-1,1])
% %yticks(-0.5:0.1:1)
% xlim([0,70])
% grid on
% xlabel("Time (h)","FontSize",12)
% ylabel("Power deviation (AU)","FontSize",12)
% saveas(delta_power_figure,sprintf("%s/delta_power_HET.svg",investigation,folder_name))
% 
% %saveas(delta_power_figure,"../Mode_variation_10_24/diag_05/delta_power_HOM.png")
% delta_AO_figure = figure();
% plot(reduced_time_hours_HET,delta_AO_HET*100,"LineWidth",2);
% grid on
% xlim([0,70])
% ylim([-30,40])
% yticks(-30:5:40)
% xlabel("Time (h)","FontSize",12)
% ylabel("Axial offset deviation (%)","FontSize",12)
% saveas(delta_AO_figure,sprintf("%s/delta_AO_HET.svg",investigation,folder_name))

end