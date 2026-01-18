function data = preparePowerComparison(data,opts)
% This function aligns two power signals (xerom vs. mscnpp) by detecting
    % fix points and shifting destructively the most delayed signal.  
    %
    % Inputs:
    %   t_signal  - time vector for signal
    %   signal    - signal to be aligned
    %   t_ref     - time vector for reference
    %   reference - reference signal
    %
    % Outputs:
    %   t_shifted, signal_shifted       - shifted signal
    %   t_ref_shifted, reference_shifted- shifted reference

    %Extract data
    X = data.xerom.power;
    R = data.mscnpp.power_signal;
    tX = data.xerom.reduced_time_hours;
    tR = data.mscnpp.t_hours;

    % Output arrays
    data.modeplotting.xerom_modes  = zeros(size(X));
    data.modeplotting.mscnpp_modes = zeros(size(R));
    data.modeplotting.xerom_t      = zeros(size(X));
    data.modeplotting.mscnpp_t     = zeros(size(R));

    tx_skip = tX(opts.skip_xerom);
    tr_skip = tR(opts.skip_mscnpp);
    tx_seg = tX(opts.skip_xerom+1:end)- tx_skip;
    tr_seg = tR(opts.skip_mscnpp+1:end)-tr_skip;
    %x_seg = X(opts.skip_xerom+1:end);
    r_detrendResults = detrend(tr_seg,R(opts.skip_mscnpp+1:end),opts.powerfit_mscnpp);
    r_seg = r_detrendResults.detrended; % Store the detrended segment for further processing
    x_detrendResults = detrend(tx_seg,X(opts.skip_xerom+1:end),opts.powerfit_xerom);
    x_seg = x_detrendResults.detrended;
    
    [x_seg, tx_seg, r_seg, tr_seg] = alignModes(x_seg, tx_seg, r_seg, tr_seg, ...
                                                opts.powerplot.n_peak_xerom, opts.powerplot.n_peak_mscnpp);

    [x_seg, tx_seg] = cutAtZeroCrossing(x_seg, tx_seg, opts.powerplot.n_zero_xerom);
    [r_seg,tr_seg] = cutAtZeroCrossing(r_seg, tr_seg, opts.powerplot.n_zero_mscnpp);

        % Save peak lists for global comparison
    [xp,xl] = findpeaks(x_seg); 
    [rp,rl] = findpeaks(r_seg); 
    % Sort peaks
    [xp, sx] = sort(xp,'descend'); xl = xl(sx); txp = tx_seg(xl);
    [rp, sr] = sort(rp,'descend'); rl = rl(sr); trp = tr_seg(rl);

    L = min(length(x_seg),length(r_seg));
    data.powerplotting.mscnpp_power(1:L) = r_seg(1:L) / abs(rp(opts.powerplot.n_peak_mscnpp));
    data.powerplotting.xerom_power(1:L) = x_seg(1:L) / xp(opts.powerplot.n_peak_xerom);
    data.powerplotting.t_xerom(1:L) = tx_seg(1:L);
    data.powerplotting.t_mscnpp(1:L) = tr_seg(1:L);
    data.powerfit.xerom_frequency = x_detrendResults.frequency;
    data.powerfit.xerom_alpha     = x_detrendResults.alpha;
    data.powerfit.mscnpp_frequency = r_detrendResults.frequency;
    data.powerfit.mscnpp_alpha    = r_detrendResults.alpha;
end

