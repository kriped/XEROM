function [time,state_values]=ivp(M)
close all; clc; clear vars;
ti=0; tf = 30*3600;
tspan = [ti,tf];
IC = zeros(1,M*3)+0.01;
% ExMode = 2;
% exmodeidx = (ExMode-1)*3+1;
% IC(exmodeidx) = 0.01;
% CR = zeros(M);
[time, state_values] = ode15s(ode,tspan,IC,[]);

% ti = 2*3600+1; tf = 5*3600;
% tspan = [ti,tf];
% IC = state_values_b(end,:);
% % CR = 3.0e+04 *[
% %   1.3416   -0.0000   -0.0000    1.1672
% %     0.0000   -0.4781   -0.0000   -0.0000
% %    -0.0000    0.0000    0.4781    0.0000
% %    -1.1543   -0.0000   -0.0000   -1.0377];
% 
% [time_a, state_values_a] = ode15s(@ode,tspan,IC,[],CR);
% 
% time_tot = [time_b;time_a]/3600;
% 
% state_values_tot = [state_values_b;state_values_a];

%% plot
figure()

plot(time,state_values(:,1:3:end))
%ylim([-0.1,0.1])
ylabel("Amplitude (AU)")
xlabel("Time (h)")
legend("Eq mode", "rad mode 1", "rad mode 2", "ax mode")
end