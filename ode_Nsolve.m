function [time,state_values] = ode_Nsolve()


%%
load data/tempFile

%tc = 2*3600; % time (in hours) of control rod insertion
%CR = zeros(M).*(t<=tc) + PHID_CR_PHI.*(t>tc); %control rod step function

f = FunctionGen(M);
ti=0; tf = 150*3600;
tspan = [ti,tf];
%IC = zeros(1,M*3);
MinValue = -5.9e+09;
MaxValue = 5.9e+09;
rng(1235482);
IC = MinValue + (MaxValue - MinValue) * rand(1,M*3); % set all values to random values between +-0.01% of the initial perturbation 1pcm of the equilibrium flux
exmode = 4;
IC((exmode-1)*3+1) = 5.9e+11;
opts=odeset("MaxStep",10);
f_handle = eval(['@(t,s)[' f ']']);

[time, state_values] = ode15s(f_handle,tspan,IC,opts);

% figure(1)
% 
% plot(time/3600,state_values(:,(4-1)*3+1))
% %xlim([0,10])
% %title('First four flux mode amplitudes of HETROM')
% %set(gca,'yscale','log')
% ylim([-1e-7,1e-7])
% ylabel("P_{m}(t) Amplitudes (cm^{-2}*s^{-1})",'Fontsize', 14)
% xlabel("Time (h)",'Fontsize', 14)
% legend("Eq mode","ax mode", "rad mode 1", "rad mode 2", 'Location','best')
% grid on
% 
% 
% figure(2)
% yyaxis left
% plot(time/3600,state_values(:,([4 9]-1)*3+1))
% ylim([-3e-9,1e-9])
% ylabel("Amplitude [cm^{-2}s^{-1}]",'Fontsize', 14)
% yyaxis right
% plot(time/3600,state_values(:,(1-1)*3+1))
% ylim([-4E-9 , 8E-9])
% yl = ylabel("Amplitude [cm^{-2}s^{-1}]",'Fontsize', 14,"Rotation",-90);
% yl.Position(1) = 78; 
% legend([ "First axial harmonic", "Third axial harmonic" "Fundamental mode (exited)"],"Location","southwest")
% xlabel("Time [h]",'Fontsize', 14)
% xlim([0 70]);
% grid on
% hold off
% figure(2)
% plot(time/3600,state_values(:,([1, 4, 2, 3]-1)*3+1))
% xlim([0,1])
% ylim([-5E-8 , 5E-7])
% %title('First four flux mode amplitudes of HETROM')
% ylabel("P_{m}(t) Amplitudes (cm^{-2}*s^{-1})",'Fontsize', 14)
% xlabel("Time (h)",'Fontsize', 14)
% legend("Eq mode","ax mode", "rad mode 1", "rad mode 2", 'Location','best')
%grid on
end