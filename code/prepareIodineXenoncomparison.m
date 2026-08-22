function data = prepareIodineXenoncomparison(data,opts)

extra_skip_mscnpp = 65;
%% define variables
XXe_seg(:) = data.xerom.mean_xenon_concentration(opts.skip_xerom+1:end);
XI_seg(:)  = data.xerom.mean_iodine_concentration(opts.skip_xerom+1:end);
tX  = data.xerom.reduced_time_hours;
RXe_seg(:) = data.mscnpp.mean_xenon(opts.skip_mscnpp+1+extra_skip_mscnpp:end);
RI_seg(:)  = data.mscnpp.mean_iodine(opts.skip_mscnpp+1+extra_skip_mscnpp:end);
tR  = data.mscnpp.t_hours;

tx_skip = tX(opts.skip_xerom+1);
tr_skip = tR(opts.skip_mscnpp+1+extra_skip_mscnpp);
tx_segX = tX(opts.skip_xerom+1:end)- tx_skip;
tr_segX = tR(opts.skip_mscnpp+1+extra_skip_mscnpp:end)-tr_skip;
tx_segI = tx_segX;
tr_segI = tr_segX;

XXe_detrendResults  = detrend(tx_segX,XXe_seg  , opts.xenonfit_xerom);
XI_detrendResults   = detrend(tx_segI,XI_seg   , opts.iodinefit_xerom);
XXe_seg             = XXe_detrendResults.detrended;
XI_seg              = XI_detrendResults.detrended;

RXe_detrendResults  = detrend(tr_segX,RXe_seg , opts.xenonfit_mscnpp);
RI_detrendResults   = detrend(tr_segI,RI_seg  , opts.iodinefit_mscnpp);
RXe_seg             = RXe_detrendResults.detrended;
RI_seg              = RI_detrendResults.detrended;

%start the signal a zero crossing
[XXe_seg, tx_segX]  = cutAtZeroCrossing(XXe_seg  , tx_segX, opts.xenonplot.n_zero_xerom);
[XI_seg,tx_segI]    = cutAtZeroCrossing(XI_seg   , tx_segI, opts.iodineplot.n_zero_xerom);
[RXe_seg, tr_segX]  = cutAtZeroCrossing(RXe_seg  , tr_segX, opts.xenonplot.n_zero_mscnpp);
[RI_seg,tr_segI]    = cutAtZeroCrossing(RI_seg   , tr_segI, opts.iodineplot.n_zero_mscnpp);

% normalise XXe_seg, XI_seg, RXe_seg, RI_seg by their maximum value
% Normalize the segments by their maximum values
%XXe_seg = XXe_seg / max(XXe_seg);
%XI_seg  = XI_seg  / max(XI_seg);
%RXe_seg  = RXe_seg / max(RXe_seg);
%RI_seg   = RI_seg  / max(RI_seg);

XXe_detrendResults  = detrend(tx_segX,  XXe_seg, opts.xenonfit_xerom);
XI_detrendResults   = detrend(tx_segI,  XI_seg, opts.iodinefit_xerom);
XXe_seg             = XXe_detrendResults.detrended;
XI_seg              = XI_detrendResults.detrended;
RXe_detrendResults  = detrend(tr_segX,  RXe_seg, opts.xenonfit_mscnpp);
RI_detrendResults   = detrend(tr_segI,  RI_seg, opts.iodinefit_mscnpp);
RXe_seg             = RXe_detrendResults.detrended;
RI_seg              = RI_detrendResults.detrended;
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
    data.xenonplotting.mscnpp = RXe_seg / RXep(opts.xenonplot.n_peak_mscnpp);
    data.xenonplotting.xerom  = XXe_seg / XXep(opts.xenonplot.n_peak_xerom);
    data.xenonplotting.t_xerom       = tx_segX;
    data.xenonplotting.t_mscnpp      = tr_segX;
    data.xenonfit.xerom_frequency   = XXe_detrendResults.frequency;
    data.xenonfit.xerom_alpha       = XXe_detrendResults.alpha;
    data.xenonfit.xerom_phase       = XXe_detrendResults.phase;
    data.xenonfit.xerom_amplitude   = XXe_detrendResults.amplitude;
    data.xenonfit.mscnpp_frequency  = RXe_detrendResults.frequency;
    data.xenonfit.mscnpp_alpha      = RXe_detrendResults.alpha;
    data.xenonfit.mscnpp_phase      = RXe_detrendResults.phase;
    data.xenonfit.mscnpp_amplitude  = RXe_detrendResults.amplitude;
    data.xenonfit.xerom_relL2percent     = sqrt(XXe_detrendResults.resnorm)/norm(XXe_seg)*100;
    data.xenonfit.mscnpp_relL2percent    = sqrt(RXe_detrendResults.resnorm)/norm(RXe_seg)*100;
    
    % Save Iodine data
    data.iodineplotting.mscnpp= RI_seg / RIp(opts.iodineplot.n_peak_mscnpp);
    data.iodineplotting.xerom = XI_seg / XIp(opts.iodineplot.n_peak_xerom);
    data.iodineplotting.t_xerom     = tx_segI;
    data.iodineplotting.t_mscnpp    = tr_segI;
    data.iodinefit.xerom_frequency  = XI_detrendResults.frequency;
    data.iodinefit.xerom_alpha      = XI_detrendResults.alpha;
    data.iodinefit.mscnpp_frequency = RI_detrendResults.frequency;
    data.iodinefit.mscnpp_alpha     = RI_detrendResults.alpha;
    data.iodinefit.mscnpp_amplitude = RI_detrendResults.amplitude;
    data.iodinefit.xerom_phase      = XI_detrendResults.phase;
    data.iodinefit.xerom_amplitude  = XI_detrendResults.amplitude;
    data.iodinefit.mscnpp_phase     = RI_detrendResults.phase;
    data.iodinefit.xerom_resnorm    = XI_detrendResults.resnorm;
    data.iodinefit.mscnpp_resnorm   = RI_detrendResults.resnorm;
    data.iodinefit.xerom_relL2percent     = sqrt(XI_detrendResults.resnorm)/norm(XI_seg)*100;
    data.iodinefit.mscnpp_relL2percent    = sqrt(RI_detrendResults.resnorm)/norm(RI_seg)*100;
end