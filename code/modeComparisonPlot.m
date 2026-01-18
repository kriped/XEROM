function modeComparisonPlot(data,resultsDir,opts)

M = data.M;
xerom_modes = data.modeplotting.xerom_modes;
mscnpp_modes = data.modeplotting.mscnpp_modes;
t_xerom = data.modeplotting.xerom_t;
t_mscnpp = data.modeplotting.mscnpp_t;
xerom_peaks = data.modeplotting.xerom_peaks;
% xerom_peak_t = data.modeplotting.xerom_peak_t;
mscnpp_peaks = data.modeplotting.mscnpp_peaks;
% mscnpp_peak_t = data.modeplotting.mscnpp_peak_t;
xerom_peaks_fundamental = xerom_peaks(:,1);
xerom_peaks_first = xerom_peaks(:,2);
xerom_peaks_second = xerom_peaks(:,7);
mscnpp_peaks_fundamental = mscnpp_peaks(:,1);
mscnpp_peaks_first = mscnpp_peaks(:,2);
mscnpp_peaks_second = mscnpp_peaks(:,7);
for mode = 1:M
    mode_figure = figure('Position', get(0, 'Screensize'));
    legendEntry = [sprintf("Mode %i XEROM", mode-1), sprintf("Mode %i MScNPP", mode-1)];
    %t_plot_xerom = reduced_time_hours(1:size(neutron_modes_normalised(:,mode), 1));
%    mode_time = t_a_shifted(1:find(t_a_shifted(:,mode), 1, 'last'), mode);
    
    mode_amplitude_xerom = xerom_modes(1:find(xerom_modes(:,mode),1,"last"),mode);
    mode_amplitude_mscnpp = mscnpp_modes(1:find(mscnpp_modes(:,mode),1,"last"),mode);
    if mode == 1
        mode1_xerom = mode_amplitude_xerom;
        mode1_mscnpp = mode_amplitude_mscnpp;
    elseif  mode == 2
        mode2_xerom = mode_amplitude_xerom;
        mode2_mscnpp = mode_amplitude_mscnpp;
    elseif mode == 7
        mode3_xerom = mode_amplitude_xerom;
        mode3_mscnpp = mode_amplitude_mscnpp;
    end
    L_xerom = length(mode_amplitude_xerom);
    t_xerom_i = t_xerom(1:L_xerom,mode);
    L_mscnpp = length(mode_amplitude_mscnpp);
    t_mscnpp_i = t_mscnpp(1:L_mscnpp,mode);
    %plot(mode_time, mode_amplitude, "LineWidth", 2)
    plot(t_xerom_i, mode_amplitude_xerom, "LineWidth", opts.plotting.LineWidth)
    hold on
    %plot(mode_time_ref, mode_amplitude_ref, "LineWidth", 2);
    plot(t_mscnpp_i, mode_amplitude_mscnpp, "LineWidth", opts.plotting.LineWidth);
    tlim = min(t_mscnpp_i(end),t_xerom_i(end));
    xlim([0,tlim]);
    ylim([-1.2,1.2])
    ylabel('Amplitude [AU]', 'FontSize', opts.plotting.FontSize);
    xlabel('Time [h]', 'FontSize', opts.plotting.FontSize);
    yline(0)
    grid on
    legend(legendEntry)
    ax2 = gca; ax2.FontSize = opts.plotting.FontSize;
    saveas(mode_figure, resultsDir + sprintf("mode_comp_%i.fig", mode-1))
    saveas(mode_figure, resultsDir + sprintf("mode_comp_%i.pdf", mode-1))
    saveas(mode_figure, resultsDir + sprintf("mode_comp_%i.png", mode-1))
    close(mode_figure)
end
% ------------------------------------
% Plot Relative amplitude between mode 0 and 1
% ------------------------------------
% Calculate relative amplitudes for the first two modes
LX1 = min( find(xerom_modes(:,1), 1, 'last'), ...
          find(xerom_modes(:,2), 1, 'last') );
LX2 = min( find(xerom_modes(:,1), 1, 'last'), ...
          find(xerom_modes(:,7), 1, 'last') );
LR1 = min( find(mscnpp_modes(:,1), 1, 'last'), ...
          find(mscnpp_modes(:,2), 1, 'last') );
LR2 = min( find(mscnpp_modes(:,1), 1, 'last'), ...
          find(mscnpp_modes(:,7), 1, 'last') );
X1 = [mode1_xerom(1:LX1,1), mode2_xerom(1:LX1)];
R1 = [mode1_mscnpp(1:LR1,1), mode2_mscnpp(1:LR1)];
X2 = [mode1_xerom(1:LX2,1), mode3_xerom(1:LX2)];
R2 = [mode1_mscnpp(1:LR2,1), mode3_mscnpp(1:LR2)];
tX1 = t_xerom_i(1:LX1);
tR1 = t_mscnpp_i(1:LR1);
tX2 = t_xerom_i(1:LX2);
tR2 = t_mscnpp_i(1:LR2);

% Call the function to plot relative amplitudes  
plotRelativeAmplitudes(X1, R1, tX1, tR1);
plotRelativeAmplitudes(X2, R2, tX2, tR2);
plotRelativePeaks(xerom_peaks_fundamental, xerom_peaks_first, ...
                  mscnpp_peaks_fundamental, mscnpp_peaks_first);
plotRelativePeaks(xerom_peaks_fundamental,xerom_peaks_second,...
    mscnpp_peaks_fundamental,mscnpp_peaks_second);

end

% =====================================================
%            Helper Functions
% =====================================================
function plotRelativeAmplitudes(X,R,tX,tR)
    if ~isempty(X) && ~isempty(R)
        figure;
        plot(tX, X(:,2)./X(:,1),'DisplayName','XEROM rel');
        hold on;
        plot(tR, R(:,2)./R(:,1),'DisplayName','MSCNPP rel');
        ylim([0 2]);
        xlabel('Time (h)');
        ylabel('Relative Amplitude');
        title('Relative Mode Amplitudes');
        legend show; grid on; hold off;
    end
end

function plotRelativePeaks(xfund, xfirst, rfund, rfirst)
    if isempty(xfund) || isempty(xfirst) || isempty(rfund) || isempty(rfirst)
        return;
    end

    % Remove trailing zeros
    xfund = xfund(find(xfund, 1, 'last'):end);
    xfirst = xfirst(find(xfirst, 1, 'last'):end);
    rfund = rfund(find(rfund, 1, 'last'):end);
    rfirst = rfirst(find(rfirst, 1, 'last'):end);

    m = min(length(xfirst),length(xfund));
    n = min(length(rfirst),length(rfund));
    L = min(m,n);
    xr = xfirst(1:L) ./ xfund(1:L);
   
    rr = rfirst(1:L) ./ rfund(1:L);
 
    
    % Calculate relative error
    relative_error = abs(xr - rr) ./ max(xr, rr);
    
    figure;
    plot(xr,'DisplayName','XEROM peak ratio','Marker','o'); hold on;
    plot(rr,'DisplayName','MSCNPP peak ratio','Marker','o');
    ylabel('Relative amplitude');
    % Create a second y-axis for the relative error
    yyaxis right;
    plot(relative_error, 'DisplayName', 'Relative Error', 'LineStyle', '--',"Marker",'o', 'Color', 'k');
    ylabel('Relative Error');
    xlim([1,L])
    ax = gca; % Get current axes
    ax.YColor = 'k'; % Set the color of the right y-axis to black
    
    xlabel('Peak index');
    
    title('Peak ratios: fundamental vs first mode');
    legend show; 
    legend('Location', 'best'); % Improve the position of the legend
    grid on;
end