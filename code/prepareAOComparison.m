function data = prepareAOComparison(data,opts)


    %Extract data
    X = data.xerom.AO;
    R = data.mscnpp.AO;
    tX = data.xerom.reduced_time_hours;
    tR = data.mscnpp.t_hours;

    % Output arrays
    %data.AOplotting.AO_xerom  = zeros(size(X));
    %data.AOplotting.AO_mscnpp = zeros(size(R));
    %data.AOplotting.xerom_t   = zeros(size(X));
    %data.AOplotting.mscnpp_t  = zeros(size(R));

    tx_skip = tX(opts.skip_xerom);
    tr_skip = tR(opts.skip_mscnpp);
    tx_seg = tX(opts.skip_xerom+1:end)- tx_skip;
    tr_seg = tR(opts.skip_mscnpp+1:end)-tr_skip;

    %x_seg = X(opts.skip_xerom+1:end);
    r_detrendResults = detrend(tr_seg,R(opts.skip_mscnpp+1:end),opts.AOfit_mscnpp);
    r_seg = r_detrendResults.detrended; % Store the detrended segment for further processing
    x_detrendResults = detrend(tx_seg,X(opts.skip_xerom+1:end),opts.AOfit_xerom);
    x_seg = x_detrendResults.detrended;

    [x_seg, tx_seg, r_seg, tr_seg] = alignModes(x_seg, tx_seg, r_seg, tr_seg, ...
                                                opts.AOplot.n_peak_xerom, opts.AOplot.n_peak_mscnpp);
    
    [x_seg, tx_seg] = cutAtZeroCrossing(x_seg, tx_seg, opts.AOplot.n_zero_xerom);
    [r_seg,tr_seg] = cutAtZeroCrossing(r_seg, tr_seg, opts.AOplot.n_zero_mscnpp);

        % Save peak lists for global comparison
    [xp,xl] = findpeaks(x_seg); 
    [rp,rl] = findpeaks(r_seg); 
    % Sort peaks
    [xp, sx] = sort(xp,'descend'); xl = xl(sx); txp = tx_seg(xl);
    [rp, sr] = sort(rp,'descend'); rl = rl(sr); trp = tr_seg(rl);

    L = min(length(x_seg),length(r_seg));
    % Trim the data to remove the trailing zeros
    x_seg = x_seg(x_seg ~= 0);
    r_seg = r_seg(r_seg ~= 0);
    tx_seg = tx_seg(1:length(x_seg));
    tr_seg = tr_seg(1:length(r_seg));
    data.AOplotting.AO_xerom = x_seg(1:L) / xp(opts.AOplot.n_peak_xerom);
    data.AOplotting.AO_mscnpp = r_seg(1:L) / rp(opts.AOplot.n_peak_mscnpp);
    data.AOplotting.xerom_t = tx_seg(1:L);
    data.AOplotting.mscnpp_t = tr_seg(1:L);
    data.AOfit.xerom_frequency = x_detrendResults.frequency;
    data.AOfit.xerom_alpha     = x_detrendResults.alpha;
    data.AOfit.mscnpp_frequency = r_detrendResults.frequency;
    data.AOfit.mscnpp_alpha    = r_detrendResults.alpha;
end

