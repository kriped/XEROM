function [period, exponent] = Signal_processing(time,signal,PlotText)

% Plot the original signal
figure;
plot(time, signal,"LineWidth",2);
hold on;
xlabel('Time (h)');
ylabel(PlotText);

xlim([0,time(end)])

% Find the peaks of the signal
[peaks, loc_peak] = findpeaks(signal,time,"MinPeakDistance",15);
peaks = peaks(peaks>=0);
loc_peak = loc_peak(peaks>=0);
% Plot the identified peaks
plot(loc_peak, peaks, 'ro', 'MarkerFaceColor', 'r');

% Fit an exponential decay function to the peaks
fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [max(signal)*1.3, 0.02]); % Initial guesses
fitType = fittype('A*exp(-b*x)', 'independent', 'x', 'coefficients', {'A', 'b'});
[fitResult, gof] = fit(loc_peak(:), peaks(:), fitType, fitOptions);

% Display fitting results
disp('Fitting Results:');
disp(fitResult);
disp('Goodness of Fit:');
disp(gof);

% Plot the fitted decay curve
%t_fit = linspace(min(loc_peak), max(loc_peak), 500); % Fine time grid for fitting curve
plot(time, fitResult.A * exp(-fitResult.b * time), 'k--', 'LineWidth', 2);

% Calculate the period of oscillation
peak_intervals = diff(loc_peak); % Time intervals between consecutive peaks
avg_period = mean(peak_intervals); % Average period
disp(['Average Period of Oscillation: ', num2str(avg_period), ' h']);
disp(['Decay Constant: ', num2str(-fitResult.b), ' h'])
period = num2str(avg_period);
exponent = fitResult.b;

% Get the current axis limits to adjust text position accordingly
    ax = gca;
    x_limits = ax.XLim;
    y_limits = ax.YLim;
    
    annotation_pos_x = x_limits(1) + 0.1 * (x_limits(2) - x_limits(1));
    annotation_pos_y = y_limits(1) + 0.9 * (y_limits(2) - y_limits(1));
% Annotate the exponent and period on the plot
annotationText = ['Exponent: ', num2str(-fitResult.b, '%.3f'), ' h^{-1}', newline, ...
    'Period: ', num2str(avg_period, '%.1f'), ' h'];

% Add annotation at a specific position (adjust coordinates for desired location)
    text(annotation_pos_x, ...
         annotation_pos_y, ...
         annotationText, 'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Finalize plot
legend(PlotText, 'Peaks', 'Fitted Curve');
grid on
hold off;
end