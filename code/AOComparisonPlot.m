function AOComparisonPlot(data, resultsDir,opts)

% Extract time and AO data from the struct
t_xerom = data.AOplotting.xerom_t;
AO_xerom = data.AOplotting.AO_xerom;
t_mscnpp = data.AOplotting.mscnpp_t;
AO_mscnpp = data.AOplotting.AO_mscnpp;

%tlim = min(t_mscnpp(end),t_xerom(end));
% Create a figure for the plot
figure; clf; hold off

% Plot AO vs time
plot(t_xerom, AO_xerom, 'LineWidth', opts.plotting.LineWidth);
hold on
plot(t_mscnpp,AO_mscnpp, 'LineWidth', opts.plotting.LineWidth)
xlabel('Time','FontSize',opts.plotting.FontSize);
ylabel('AO',opts.plotting.FontSize);
tlim = min(t_mscnpp(end),t_xerom(end));
xlim([0,tlim]);
ylim([-1.2,1.2])
legend("XEROM","MScNPP",fontsize=opts.plotting.FontSize)
grid on;

% Save the plot to the specified results directory
saveas(gcf, fullfile(resultsDir, 'AO_comparison.png'));
saveas(gcf, fullfile(resultsDir, 'AO_comparison.pdf'));
saveas(gcf, fullfile(resultsDir, 'AO_comparison.fig'));
end