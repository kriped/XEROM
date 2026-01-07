function powerComparisonPlot(data,resultsDir,opts)
% p = inputParser;
% p.addParameter("skip",1)
% 
% p.parse(varargin{:});
% skip = p.Results.skip;

t_xerom = data.powerplotting.t_xerom;
t_mscnpp = data.powerplotting.t_mscnpp;
power_xerom = data.powerplotting.xerom_power;
power_mscnpp = data.powerplotting.mscnpp_power;

power_figure = figure();
% max_t = min(t_xerom(end), t_mscnpp(end));
% [~, xerom_power_idx] = min(abs(t_hours_shifted_power - max_t));
% [~, ref_power_idx] = min(abs(Reference.t_power - max_t));
plot(t_xerom(opts.plotting.skip:end), power_xerom(opts.plotting.skip:end), "LineWidth", 2)
hold on
plot(t_mscnpp(opts.plotting.skip:end), power_mscnpp(opts.plotting.skip:end), "LineWidth", 2)
tlim = min(t_mscnpp(end),t_xerom(end));
xlim([0,tlim]);
ylim([-1.2,1.2])
grid on
xlabel("Time (h)", "FontSize", 16)
ylabel("Power (AU)", "FontSize", 16)
ax2 = gca; ax2.FontSize = 16;
legend("XEROM", "MscNPP")
saveas(power_figure, sprintf("%s/power_comparison.png", resultsDir))
saveas(power_figure, sprintf("%s/power_comparison.pdf", resultsDir))
saveas(power_figure, sprintf("%s/power_comparison.fig", resultsDir))
