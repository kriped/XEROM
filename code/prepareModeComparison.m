function data = prepareModeComparison(data,opts)
% --------------------------
% Extract input fields
% --------------------------
X = data.xerom.neutronModes;
R = data.mscnpp.a;
tX = data.xerom.reduced_time_hours;
tR = data.mscnpp.t_hours;

maxPeaks = opts.maxPeaks;
N_modes = size(X,2);

% Output arrays
data.modeplotting.xerom_modes  = zeros(size(X));
data.modeplotting.mscnpp_modes = zeros(size(R));
data.modeplotting.xerom_t      = zeros(size(X));
data.modeplotting.mscnpp_t     = zeros(size(R));
data.modeplotting.mscnpp_peaks = zeros(N_modes,maxPeaks);
data.modeplotting.xerom_peaks  = zeros(N_modes,maxPeaks);
data.modeplotting.xerom_peak_t = zeros(N_modes,maxPeaks);
data.modeplotting.mscnpp_peak_t= zeros(N_modes,maxPeaks);
% --------------------------
% Process each mode
% --------------------------
for i = 1:N_modes
    tx_skip = tX(opts.skip_xerom);
    tr_skip = tR(opts.skip_mscnpp);
    tx_seg = tX(opts.skip_xerom+1:end)- tx_skip;
    tr_seg = tR(opts.skip_mscnpp+1:end)-tr_skip;
    %x_seg = X(opts.skip_xerom+1:end,i);
    r_detrendResults = detrend(tr_seg,R(opts.skip_mscnpp+1:end,i),opts.modefit);
    r_seg = r_detrendResults.detrended; % Store the detrended segment for further processing
    x_detrendResults = detrend(tx_seg,X(opts.skip_mscnpp+1:end,i),opts.modefit);
    x_seg = x_detrendResults.detrended; % Store the detrended segment for further processing
    
    [x_seg, tx_seg] = cutAtZeroCrossing(x_seg, tx_seg, opts.modeplot.n_zero_xerom);
    [r_seg,tr_seg] = cutAtZeroCrossing(r_seg, tr_seg, opts.modeplot.n_zero_mscnpp);
    %[x_seg, tx_seg, r_seg, tr_seg] = alignModes(X(:,i), tX, R(:,i), tR, ...
    %                                            opts.n_peak_xerom, opts.n_peak_mscnpp);
    [x_seg, tx_seg, r_seg, tr_seg] = alignModes(x_seg, tx_seg, r_seg, tr_seg, ...
                                                opts.modeplot.n_peak_xerom, opts.modeplot.n_peak_mscnpp);
    % Save peak lists for global comparison
    [xp,xl] = findpeaks(x_seg); 
    [rp,rl] = findpeaks(r_seg); 
    % Sort peaks
    [xp, sx] = sort(xp,'descend'); xl = xl(sx); txp = tx_seg(xl);
    [rp, sr] = sort(rp,'descend'); rl = rl(sr); trp = tr_seg(rl);

    % Normalize and store
    L = min(length(x_seg),length(r_seg));
    data.modeplotting.xerom_modes (1:L,i) = x_seg(1:L) / xp(opts.modeplot.n_peak_xerom);
    data.modeplotting.mscnpp_modes(1:L,i) = r_seg(1:L) / abs(rp(opts.modeplot.n_peak_mscnpp));
    data.modeplotting.xerom_t(1:L,i)      = tx_seg(1:L);
    data.modeplotting.mscnpp_t(1:L,i)     = tr_seg(1:L);
    data.modeplotting.xerom_peaks(1:length(xp),i)   = xp;
    data.modeplotting.mscnpp_peaks(1:length(rp),i)  = rp;
    data.modeplotting.xerom_peak_t(1:length(txp),i) = txp;
    data.modeplotting.mscnpp_peak_t(1:length(trp),i)= trp;
end



end






