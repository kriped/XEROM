function IodineComparisonPlot(data, resultsDir,opts)

% Extract time and Iodine data from the struct
t_xerom = data.iodineplotting.t_xerom;
I_xerom = data.iodineplotting.xerom;
t_mscnpp = data.iodineplotting.t_mscnpp;
I_mscnpp = data.iodineplotting.mscnpp;

%tlim = min(t_mscnpp(end),t_xerom(end));
% Create a figure for the plot
figure; clf; hold off

% Plot Iodine concentration vs time
t_end = min(t_mscnpp(end),t_xerom(end));
%Find indexs of L
L_xerom = find(t_xerom <= t_end, 1, 'last');
L_mscnpp = find(t_mscnpp <= t_end,1,'last');
plot(t_xerom(1:L_xerom), I_xerom(1:L_xerom), 'LineWidth', opts.plotting.LineWidth)
hold on
plot(t_mscnpp(1:L_mscnpp),I_mscnpp(1:L_mscnpp), 'LineWidth', opts.plotting.LineWidth)
xlabel('Time','FontSize',opts.plotting.FontSize);
ylabel('Iodine Concentration (AU)','FontSize',opts.plotting.FontSize);
tlim = min(t_mscnpp(end),t_xerom(end));
xlim([0,tlim])
ylim([-1.2,1.2])
pbaspect([4,3,1])
legend("XEROM","MScNPP",fontsize=opts.plotting.FontSize)
grid on;

% Save the plot to the specified results directory
saveas(gcf, fullfile(resultsDir, 'iodine_comparison.png'));
saveas(gcf, fullfile(resultsDir, 'iodine_comparison.pdf'));
saveas(gcf, fullfile(resultsDir, 'iodine_comparison.fig'));
end