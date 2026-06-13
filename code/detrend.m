function results = detrend(t, y,opts)
    
% DETREND DAMPED SINUSOID
% This function detrends a given signal using a damped sinusoidal model.
%
% Inputs:
%   t      - Time vector (1D array)
%   y      - Signal vector (1D array) corresponding to time vector t
%   opts   - Options structure containing:
%            .f0      - Frequency guess (default: "Default")
%            .A0      - Amplitude guess (default: "Default")
%            .alpha0  - Decay rate guess (default: "Default")
%            .phi0    - Phase guess (default: "Default")
%            .c0      - Intercept guess (default: "Default")
%            .c1      - Slope guess (default: "Default")
%            .showPlot - Boolean to display plots (default: false)
%            .name     string containing the plot name
%
% Outputs:
%   results - Structure containing:
%             .detrended  - Detrended signal
%             .modelFit   - Fitted model signal
%             .linTrend    - Linear trend component
%             .slope      - Slope of the linear trend
%             .intersect  - Intercept of the linear trend
    [y_detrended, modelFit, params,linTrend,resnorm,residual] = detrend_damped_sinusoid(t,y,opts);

    % Plots
    if opts.showPlot
        figure('Name',opts.name,'Color','w');
        tiledlayout(2,2, 'Padding','compact','TileSpacing','compact');

        nexttile;
        plot(t,y,'g-','LineWidth',1); hold on;
        % linear fit line:
        plot(t,modelFit,'m.','LineWidth',1)
        plot(t, linTrend, 'k--', 'LineWidth', 1);
        yline(0);
        grid on
        xlabel('Time'); ylabel('Signal');
        title(sprintf('Original signal with linear trend fit'));
        legend('y','curve fit','linear fit','Location','best'); grid on;

        nexttile;
        plot(t, y_detrended, 'g-', 'LineWidth', 1); hold on
        yline(0, 'k-'); grid on;
        scatter(t, modelFit(:) - linTrend(:), 10, 'm', 'filled');
        xlabel('Time'); ylabel('Detrended signal');
        title('Detrended (should oscillate around 0)');
        legend('Detrended signal', 'Fit minus linear trend');

        % Display the linear function parameters
        linearFunction = sprintf('Linear function: y = %.2f + %.2f*t', params(5), params(6));
        annotation('textbox', [0.15, 0.1, 0.3, 0.1], 'String', linearFunction, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    end
    %--------- 7) Package results ----------
    results = struct();
    results.detrended        = y_detrended;
    results.modelFit         = modelFit;
    results.linTrend         = linTrend;
    results.frequency        = params(3);
    results.alpha            = params(2);
    results.slope            = params(6); % per unit time
    results.intersect        = params(5);
    results.resnorm          = resnorm;
    results.residual         = residual;
end




function [y_detrended, modelFit, params,linTrend,RESNORM,RESIDUAL] = detrend_damped_sinusoid(t, y,opts)
%    t = t(opts.skip:end);
%    y = y(opts.skip:end);
    
    % ----- Initial guesses -----
    % Frequency guess via FFT peak
    if num2str(opts.f0) == "Default"
        Y = fft(y);
        [~, idx] = max(abs(Y(2:floor(end/2)))); 
        f0 = (idx/length(t)) / (t(2)-t(1));  
    else
        f0 = opts.f0;
    end

    if num2str(opts.fL) == "Default"
        fL = 0.01;
    else
        fL = opts.fL;
    end

    if num2str(opts.fU) == "Default"
        fU = 0.1;
    else
        fU = opts.fU;
    end

    if num2str(opts.A0) == "Default"
        A0 = (max(y) - min(y)) / 2;   % amplitude guess
    else
        A0 = opts.A0;
    end

    if num2str(opts.alpha0) == "Default"
        alpha0 = -0.01;                % slow decay guess
    else 
        alpha0 = opts.alpha0;
    end

    if num2str(opts.alphaU) == "Default"
        alphaU = 0;                % upper bound for decay rate guess
    else
        alphaU = opts.alphaU;
    end

    if num2str(opts.alphaL) == "Default"
        alphaL = -1.0;                % lower bound for decay rate guess
    else
        alphaL = opts.alphaL;
    end

    if num2str(opts.phi0) == "Default"
        phi0 = 0;                     % phase guess
    else
        phi0 = opts.phi0;
    end

    if num2str(opts.c0) == "Default"
        c0 = mean(y);                 % intercept guess
    else
        c0 = opts.c0;
    end

    if num2str(opts.c1) == "Default"
        c1 = 0;                       % slope guess
    else
        c1 = opts.c1;
    end
    
    p0 = [A0, alpha0, f0, phi0, c0, c1];

    % ----- Model function -----
    modelFun = @(p,t) ...
        p(5) + p(6).*t + ... % trend
        p(1).*exp(p(2).*t) .* sin(2*pi*p(3).*t + p(4)); % damped or forced sine

    % ----- Nonlinear least squares fit -----
    opts = optimoptions('lsqcurvefit','Display','off');
    lb = [-inf,      alphaL,      fL,   -2*pi, -Inf, -1e4];
    ub = [Inf,  alphaU,  fU,    2*pi,  Inf,  1e4];

    [params,RESNORM,RESIDUAL] = lsqcurvefit(modelFun, p0, t(:), y(:), lb, ub, opts);

    % ----- Extract model fit -----
    modelFit = modelFun(params, t(:));
    linTrend = params(5)+params(6)*t;
    % ----- Detrended signal -----
    y_detrended = y(:) - linTrend(:);



end
