function data = computePowerScaling(data,opts)

% computePowerScaling - Computes the scaling factor between two power signals.
%
% Syntax:
%   data = computePowerScaling(data, opts)
%
% Inputs:
%   data - A struct containing the following fields:
%       .xerom.Dpower - Power data from the XEROM signal.
%       .mscnpp.Dpower - Power data from the MScNPP signal.
%       .xerom.reduced_time_hours - Time vector for the XEROM signal in hours.
%       .mscnpp.t_hours - Time vector for the MScNPP signal in hours.
%   opts - A struct containing options for processing, including:
%       .skip_xerom - Number of initial data points to skip for XEROM.
%       .skip_mscnpp - Number of initial data points to skip for MScNPP.
%       .powerfit - Parameters for detrending the signals.
%       .powerplot.n_zero_xerom - Number of zero crossings to consider for XEROM.
%       .powerplot.n_zero_mscnpp - Number of zero crossings to consider for MScNPP.
%       .powerplot.n_peak_xerom - Number of peaks to consider for XEROM.
%       .powerplot.n_peak_mscnpp - Number of peaks to consider for MScNPP.
%
% Outputs:
%   data - The input struct with an additional field:
%       .scaling.xerom_to_mscnpp - The computed scaling factor from XEROM to MScNPP.
%
% Description:
%   This function processes the power data from two signals (XEROM and MScNPP),
%   detrends the data, identifies peaks, and computes a scaling factor based on
%   the specified peaks. The scaling factor is then stored in the output data struct.

xerom_Dpower = data.xerom.Dpower;
mscnpp_Dpower = data.mscnpp.Dpower;
tX = data.xerom.reduced_time_hours;
tR = data.mscnpp.t_hours;
X = xerom_Dpower;
R = mscnpp_Dpower;
tx_seg = tX(opts.skip_xerom+1:end)-tX(opts.skip_xerom);
tr_seg = tR(opts.skip_mscnpp+1:end)-tR(opts.skip_mscnpp);
%x_seg = X(opts.skip_xerom+1:end);
r_detrendResults = detrend(tr_seg,R(opts.skip_mscnpp+1:end),opts.powerfit);
r_seg = r_detrendResults.detrended; % Store the detrended segment for further processing
x_detrendResults = detrend(tx_seg,X(opts.skip_xerom+1:end),opts.powerfit);
x_seg = x_detrendResults.detrended;    
[x_seg, tx_seg] = cutAtZeroCrossing(x_seg, tx_seg, opts.powerplot.n_zero_xerom);
[r_seg,tr_seg] = cutAtZeroCrossing(r_seg, tr_seg, opts.powerplot.n_zero_mscnpp);
[xerom_Dpower, ~, mscnpp_Dpower, ~] = alignModes(x_seg, tx_seg, r_seg, tr_seg, ...
    opts.powerplot.n_peak_xerom, opts.powerplot.n_peak_mscnpp);
% --- 1. Find peaks in both signals ---
[xerom_peaks, ~] = findpeaks(xerom_Dpower);
[mscnpp_peaks, ~] = findpeaks(mscnpp_Dpower);

if isempty(xerom_peaks) || isempty(mscnpp_peaks)
    error('No peaks found in one or both signals.');
end

% --- 2. Sort peaks ---
[xerom_peaks_sorted, ~] = sort(xerom_peaks, 'descend');
[mscnpp_peaks_sorted, ~] = sort(mscnpp_peaks, 'descend');

% --- 3. Pick fix points ---
if numel(mscnpp_peaks_sorted) < 2
    warning('Reference signal has fewer than 2 peaks. Using highest peak as fix point.');
    mscnpp_fix_point_mag = mscnpp_peaks_sorted(1);
else
    mscnpp_fix_point_mag = mscnpp_peaks_sorted(opts.powerplot.n_peak_mscnpp);
end
xerom_fix_point_mag = xerom_peaks_sorted(opts.powerplot.n_peak_xerom);

% --- 4. Compute scaling factor ---
scale_factor = mscnpp_fix_point_mag / xerom_fix_point_mag;

% --- 5. Store in data struct ---
data.scaling.xerom_to_mscnpp = scale_factor;
fprintf('Computed scaling factor (MScNPP/XEROM) = %.2f\n', scale_factor);
end
