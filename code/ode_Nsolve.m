function [time,state_values] = ode_Nsolve(input_dir)


%%
load("../input/CONSTANTS_data.mat")
load(input_dir+"PARAMETERS_data.mat")

%tc = 2*3600; % time (in hours) of control rod insertion
%CR = zeros(M).*(t<=tc) + PHID_CR_PHI.*(t>tc); %control rod step function

f = FunctionGen(M);
ti=0; tf = 200*3600;
tspan = [ti,tf];

% IC = zeros(1,M*3);
MinValue = -5.9e+09;
MaxValue = 5.9e+09;
rng(1235482);
% Include 10 modes
%Include all modes. Assign random numbers all modes
IC_values = MinValue + (MaxValue - MinValue) * rand(1,M*3); % set all values to random values between +-0.01% of the initial perturbation 1pcm of the equilibrium flux
%Include 9 modes
%IC_values = 1E11 * [5.6211    5.3118   -2.4537    2.0210    1.6172    5.7508    5.4815    2.8594    0.8657   -5.3484   -2.8468 3.7124    0.2117    3.6790   -3.7046    1.8795    4.0070   -4.3037    4.2582    1.5652    4.5605   -2.2567    3.9460    2.0000    4.5438  -3.2849   -0.0480];
%include 8 modes
%IC_values = 1E11 * [5.6211    5.3118   -2.4537    2.0210    1.6172    5.7508    5.4815    2.8594    0.8657   -5.3484   -2.8468 3.7124    0.2117    3.6790   -3.7046 4.2582    1.5652    4.5605   -2.2567    3.9460    2.0000    4.5438  -3.2849   -0.0480];
%Include 5 modes
%IC_values = 1E11 * [4.4029   -2.9283    5.5014    5.6211    5.3118   -2.4537    2.0210    1.6172    5.7508    5.4815    2.8594    0.8657   -5.3484   -2.8468    3.7124];
IC = IC_values;
exmode = 2;
IC((exmode-1)*3+1) = 3e+14;
%Include 1 mode
%IC = [3e+14  5.3118e+11   -2.4537+11];
opts=odeset("MaxStep",180);
f_handle = eval(['@(t,s)[' f ']']);

[time, state_values] = ode15s(f_handle,tspan,IC,opts);
end