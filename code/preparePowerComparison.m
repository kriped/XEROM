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
    skip_xerom  = opts.skip_xerom;
    skip_mscnpp = opts.skip_mscnpp;
    skip_mscnpp = skip_mscnpp+30; % amount of points to be removed at the beginning of the signal.
    skip_end_mscnpp = 90; % Amount of point to be removed at the tail of the signal.
    %Extract data
    x_seg       = data.xerom.power(skip_xerom+1:end);
    r_seg       = data.mscnpp.power_signal(skip_mscnpp+1:end-skip_end_mscnpp);
    t_skip      = data.xerom.reduced_time_hours(skip_xerom+1);
    tX          = data.xerom.reduced_time_hours(skip_xerom+1:end)-t_skip;
    t_skip      = data.mscnpp.t_hours(skip_mscnpp+1);
    tR          = data.mscnpp.t_hours(skip_mscnpp+1:end-skip_end_mscnpp)-t_skip;

    %Normalise x_seg and r_seg
    % Normalize the segments
    %x_seg = x_seg / norm(x_seg);
    %r_seg = r_seg / norm(r_seg);

    x_detrendResults    = detrend(tX,x_seg,opts.powerfit_xerom);
    x_seg               = x_detrendResults.detrended;
    r_detrendResults    = detrend(tR,r_seg,opts.powerfit_mscnpp);
    r_seg               = r_detrendResults.detrended; % Store the detrended segment for further processing
  
    [x_seg, tx_seg]     = cutAtZeroCrossing(x_seg, tX, opts.powerplot.n_zero_xerom);
    [r_seg,tr_seg]      = cutAtZeroCrossing(r_seg, tR, opts.powerplot.n_zero_mscnpp);

    x_detrendResults    = detrend(tx_seg,x_seg,opts.powerfit_xerom);
    x_seg               = x_detrendResults.detrended;    
    r_detrendResults    = detrend(tr_seg,r_seg,opts.powerfit_mscnpp);
    r_seg               = r_detrendResults.detrended; % Store the detrended segment for further processing
    
    [x_seg, tx_seg, r_seg, tr_seg] = alignModes(x_seg, tx_seg, r_seg, tr_seg, ...
                                                opts.powerplot.n_peak_xerom, opts.powerplot.n_peak_mscnpp);

        % Save peak lists for global comparison
    [xp,~]     = findpeaks(x_seg); 
    [rp,~]     = findpeaks(r_seg); 
    % Sort peaks
    [xp, ~ ]    = sort(xp,'descend'); %xl = xl(sx); %txp = tx_seg(xl);
    [rp, ~]    = sort(rp,'descend'); %rl = rl(sr); %trp = tr_seg(rl);

    % L = min(length(x_seg),length(r_seg));
    data.powerplotting.xerom_power  = x_seg / xp(opts.powerplot.n_peak_xerom);
    data.powerplotting.mscnpp_power = r_seg / rp(opts.powerplot.n_peak_mscnpp);
    data.powerplotting.t_xerom      = tx_seg;
    data.powerplotting.t_mscnpp     = tr_seg;
    data.powerfit.xerom_frequency   = x_detrendResults.frequency;
    data.powerfit.xerom_alpha       = x_detrendResults.alpha;
    data.powerfit.xerom_amplitude   = x_detrendResults.amplitude;
    data.powerfit.xerom_phase       = x_detrendResults.phase;
    data.powerfit.mscnpp_frequency  = r_detrendResults.frequency;
    data.powerfit.mscnpp_alpha      = r_detrendResults.alpha;
    data.powerfit.mscnpp_amplitude  = r_detrendResults.amplitude;               
    data.powerfit.mscnpp_phase      = r_detrendResults.phase;
end

