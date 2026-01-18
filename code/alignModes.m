
function [x, tx, r, tr] = alignModes(x, tx, r, tr, kx, kr)
    [xp, xl] = findpeaks(x);
    [rp, rl] = findpeaks(r);

    if isempty(xp) || isempty(rp)
        error("No peaks found for alignment.")
    end

    % Sort peaks
    [xp, sx] = sort(xp,'descend'); xl = xl(sx);
    [rp, sr] = sort(rp,'descend'); rl = rl(sr);

    tX = tx(xl(kx));
    tR = tr(rl(kr));
    dt = tX - tR;

    if dt < 0
        idx = interp1(tr,1:length(r),tX,'nearest','extrap');
        shift = rl(kr) - idx;
        r  = r (shift+1:end);
        tr = tr(1:end-shift);
    elseif dt > 0
        idx = interp1(tx,1:length(x),tR,'nearest','extrap');
        shift = xl(kx) - idx;
        x  = x (shift+1:end);
        tx = tx(1:end-shift);
    end

    L = min(length(x),length(r));
    x = x(1:L); tx = tx(1:L);
    r = r(1:L); tr = tr(1:L);
end