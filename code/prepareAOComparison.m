function data = prepareAOComparison(data,opts)
% This function aligns two Axial offset signals (xerom vs. mscnpp) by detecting
    % fix points and shifting destructively the most delayed signal.  
    %
    % Inputs:
    %   t_xerom  - time vector for signal
    %   xerom_AO    - signal to be aligned
    %   t_mscnpp    - time vector for reference
    %   mscnpp_AO - reference signal
    %
    % Outputs:
    %   t_shifted, xerom_AO_shifted       - shifted signal
    %   t_ref_shifted, mscnpp_AO_shifted- shifted reference
    
    % p = inputParser;
    % opts.addParameter("skip",1) % Determines the number of index points to skip in the output of the signal. 1 means everything is kept.
    % opts.addParameter("mscnppFixpoint",2) % Determines which peak to use as a fix point for the syncronisation (ordered by size)
    % opts.addParameter("xeromFixpoint",1) % Determines which peak to use as a fix point for the syncronisation (ordered by size)
    % opts.parse(varargin{:})
    % skip = opts.Results.skip; mscnppFixpoint= opts.Results.mscnppFixpoint; xeromFixpoint = opts.Results.xeromFixpoint;

    %Extract data
    xerom_AO = data.xerom.AO;
    mscnpp_AO = data.mscnpp.AO;
    t_xerom = data.xerom.reduced_time_hours;
    t_mscnpp = data.mscnpp.t_hours;

    shift_mscnpp = mean(mscnpp_AO);
    mscnpp_AO = mscnpp_AO - shift_mscnpp;
    shift_xerom = mean(xerom_AO);
    xerom_AO = xerom_AO-shift_xerom;

    
    % --- 1. Find peaks in both signals ---
    [xerom_peaks, xerom_loc] = findpeaks(xerom_AO);
    [mscnpp_peaks, mscnpp_loc] = findpeaks(mscnpp_AO);

     % --- 2. Handle missing peaks ---
    if isempty(xerom_peaks) || isempty(mscnpp_peaks)
        error('No peaks found in one or both signals.');
    end
        % --- 3. Sort peaks ---
    [xerom_peaks_sorted, xerom_sorted_idx] = sort(xerom_peaks, 'descend');
    xerom_loc_sorted = xerom_loc(xerom_sorted_idx);
    [mscnpp_peaks_sorted, mscnpp_sorted_idx] = sort(mscnpp_peaks, 'descend');
    mscnpp_loc_sorted = mscnpp_loc(mscnpp_sorted_idx);

     % --- 4. Pick fix point (second highest peak) ---
    if numel(mscnpp_peaks_sorted) < mscnppFixpoint
        warning('Reference signal has fewer than 2 peaks. Using highest peak as fix point.');
        mscnpp_fix_point_mag = mscnpp_peaks_sorted(1);
        mscnpp_fix_point_loc = mscnpp_loc_sorted(1);
    else
        mscnpp_fix_point_mag = mscnpp_peaks_sorted(mscnppFixpoint);
        mscnpp_fix_point_loc = mscnpp_loc_sorted(mscnppFixpoint);
    end
    xerom_fix_point_mag = xerom_peaks_sorted(xeromFixpoint);
    xerom_fix_point_loc = xerom_loc_sorted(xeromFixpoint);
    % --- 5.Normalise AO signals relative to fix points ---
    xerom_AO_normalised = xerom_AO/xerom_fix_point_mag;
    mscnpp_AO_normalised = mscnpp_AO/mscnpp_fix_point_mag;
    
    % --- 7. Get times of fix points ---
    t_peak_xerom = t_xerom(xerom_fix_point_loc);
    t_peak_mscnpp    = t_mscnpp(mscnpp_fix_point_loc);
    time_difference = t_peak_xerom - t_peak_mscnpp;

    % --- 8. Compute index shift ---
    if time_difference < 0
        % MScNPP is ahead → trim reference by interpolating
        idx_mscnpp_fix   = mscnpp_fix_point_loc;
        idx_xerom_fix   = interp1(t_mscnpp, 1:length(t_mscnpp), t_peak_xerom, 'nearest', 'extrap');
        index_shift   = abs(idx_mscnpp_fix - idx_xerom_fix);

        mscnpp_AO_trimmed = mscnpp_AO_normalised(index_shift+1:end);
        t_mscnpp_trimmed     = t_mscnpp(1:end-index_shift);
        xerom_AO_trimmed  = xerom_AO_normalised;
        t_xerom_trimmed      = t_xerom;
    elseif time_difference > 0
        % XEROM is ahead → trim reference by interpolating
        idx_mscnpp_fix   = interp1(t_xerom, 1:length(t_xerom), t_peak_mscnpp, 'nearest', 'extrap');
        idx_xerom_fix   = xerom_fix_point_loc;
        index_shift   = abs(idx_mscnpp_fix - idx_xerom_fix);

        mscnpp_AO_trimmed = mscnpp_AO_normalised;
        t_mscnpp_trimmed     = t_mscnpp;
        xerom_AO_trimmed  = xerom_AO_normalised(index_shift+1:end);
        t_xerom_trimmed      = t_xerom(1:end-index_shift);
    else
        % There is no time difference
        mscnpp_AO_trimmed = mscnpp_AO_normalised;
        t_mscnpp_trimmed     = t_mscnpp;
        xerom_AO_trimmed  = xerom_AO_normalised;
        t_xerom_trimmed      = t_xerom;
    end
     % Trim and store results
    mscnpp_AO_shifted = mscnpp_AO_trimmed(skip:end);
    xerom_AO_shifted  = xerom_AO_trimmed(skip:end);
    t_xerom_shifted   = t_xerom_trimmed(skip:end);
    t_mscnpp_shifted = t_mscnpp_trimmed(skip:end);

    data.AOplotting.mscnpp_AO = mscnpp_AO_shifted;
    data.AOplotting.xerom_AO = xerom_AO_shifted;
    data.AOplotting.t_xerom = t_xerom_shifted;
    data.AOplotting.t_mscnpp = t_mscnpp_shifted;
end