function data = prepareIodineXenoncomparison(data,opts)
%% define variables
XXe(:) = data.xerom.mean_xenon_concentration;
XI(:)  = data.xerom.mean_iodine_concentration;
tX  = data.xerom.reduced_time_hours;
RXe(:) = data.mscnpp.mean_xenon;
RI(:)  = data.mscnpp.mean_iodine;
tR  = data.mscnpp.t_hours;

tx_skip = tX(opts.skip_xerom);
tr_skip = tR(opts.skip_mscnpp);
tx_segX = tX(opts.skip_xerom+1:end)- tx_skip;
tr_segX = tR(opts.skip_mscnpp+1:end)-tr_skip;
tx_segI = tx_segX;
tr_segI = tr_segX;

XXe_detrendResults  = detrend(tx_segX,XXe(opts.skip_xerom+1:end)  , opts.xenonfit_xerom);
XI_detrendResults   = detrend(tx_segI,XI(opts.skip_xerom+1:end)   , opts.iodinefit_xerom);
XXe_seg             = XXe_detrendResults.detrended;
XI_seg              = XI_detrendResults.detrended;
RXe_detrendResults  = detrend(tr_segX,RXe(opts.skip_mscnpp+1:end) , opts.xenonfit_mscnpp);
RI_detrendResults   = detrend(tr_segI,RI(opts.skip_mscnpp+1:end)  , opts.iodinefit_mscnpp);
RXe_seg             = RXe_detrendResults.detrended;
RI_seg              = RI_detrendResults.detrended;

%start the signal a zero crossing
[XXe_seg, tx_segX]  = cutAtZeroCrossing(XXe_seg  , tx_segX, opts.xenonplot.n_zero_xerom);
[XI_seg,tx_segI]    = cutAtZeroCrossing(XI_seg   , tx_segI, opts.iodineplot.n_zero_xerom);
[RXe_seg, tr_segX]  = cutAtZeroCrossing(RXe_seg  , tr_segX, opts.xenonplot.n_zero_mscnpp);
[RI_seg,tr_segI]    = cutAtZeroCrossing(RI_seg   , tr_segI, opts.iodineplot.n_zero_mscnpp);
% Align Xenon peaks
[XXe_seg, tx_segX, RXe_seg, tr_segX]= alignModes(XXe_seg, tx_segX, RXe_seg, tr_segX, ...
    opts.xenonplot.n_peak_xerom, opts.xenonplot.n_peak_mscnpp);
% Align Iodine peaks
[XI_seg, tx_segI, RI_seg, tr_segI]  = alignModes(XI_seg, tx_segI, RI_seg, tr_segI, ...
    opts.iodineplot.n_peak_xerom, opts.iodineplot.n_peak_mscnpp);


%Find peaks
    [XXep,~] = findpeaks(XXe_seg); 
    [RXep,~] = findpeaks(RXe_seg); 
    [XIp,~]  = findpeaks(XI_seg); 
    [RIp,~]  = findpeaks(RI_seg); 
    % % Sort peaks
    % [XXep, XXes] = sort(XXep,'descend'); XXel = XXel(XXes); tXXep = tx_segX(XXel);
    % [RXep, RXes] = sort(RXep,'descend'); RXel = RXel(RXes); tRXep = tr_segX(RXel);
    % [XIp, XIs] = sort(XIp,'descend'); XIl = XIl(XIs); tXIp = tx_segI(XIl);
    % [RIp, RIs] = sort(RIp,'descend'); RIl = RIl(RIs); tRIp = tr_segI(RIl);
    % Save Xenon data
    data.xenonplotting.mscnpp_xenon = RXe_seg / RXep(opts.xenonplot.n_peak_mscnpp);
    data.xenonplotting.xerom_xenon  = XXe_seg / XXep(opts.xenonplot.n_peak_xerom);
    data.xenoplotting.t_xerom       = tx_segX;
    data.xenoplotting.t_mscnpp      = tr_segX;
    data.xenonfit.xerom_frequency   = XXe_detrendResults.frequency;
    data.xenonfit.xerom_alpha       = XXe_detrendResults.alpha;
    data.xenonfit.mscnpp_frequency  = RXe_detrendResults.frequency;
    data.xenonfit.mscnpp_alpha      = RXe_detrendResults.alpha;
    
    % Save Iodine data
    data.iodineplotting.mscnpp_xenon= RI_seg / RIp(opts.iodineplot.n_peak_mscnpp);
    data.iodineplotting.xerom_xenon = XI_seg / XIp(opts.iodineplot.n_peak_xerom);
    data.iodineplotting.t_xerom     = tx_segI;
    data.iodineplotting.t_mscnpp    = tr_segI;
    data.iodinefit.xerom_frequency  = XI_detrendResults.frequency;
    data.iodinefit.xerom_alpha      = XI_detrendResults.alpha;
    data.iodinefit.mscnpp_frequency = RI_detrendResults.frequency;
    data.iodinefit.mscnpp_alpha     = RI_detrendResults.alpha;
end