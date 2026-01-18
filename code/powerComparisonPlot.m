function powerComparisonPlot(data,resultsDir,opts)


t_xerom = data.powerplotting.t_xerom;
t_mscnpp = data.powerplotting.t_mscnpp;
power_xerom = data.powerplotting.xerom_power;
power_mscnpp = data.powerplotting.mscnpp_power;

power_figure = figure();
% max_t = min(t_xerom(end), t_mscnpp(end));
% [~, xerom_power_idx] = min(abs(t_hours_shifted_power - max_t));
% [~, ref_power_idx] = min(abs(Reference.t_power - max_t));
plot(t_xerom, power_xerom, "LineWidth", opts.plotting.LineWidth)
hold on
plot(t_mscnpp, power_mscnpp, "LineWidth", opts.plotting.LineWidth)
tlim = min(t_mscnpp(end),t_xerom(end));
xlim([0,tlim]);
ylim([-1.2,1.2])
grid on
xlabel("Time (h)", "FontSize", opts.plotting.FontSize)
ylabel("Power (AU)", "FontSize", opts.plotting.FontSize)
ax2 = gca; ax2.FontSize = opts.plotting.FontSize;
legend("XEROM", "MscNPP")
saveas(power_figure, sprintf("%s/power_comparison.png", resultsDir))
saveas(power_figure, sprintf("%s/power_comparison.pdf", resultsDir))
saveas(power_figure, sprintf("%s/power_comparison.fig", resultsDir))
