function [seg, tseg] = cutAtZeroCrossing(vec, tvec, n_zero)
    z = find(vec(1:end-1).*vec(2:end) < 0);
    if isempty(z)
        start = 1;
    else
        start = z(min(n_zero,length(z))) + 1;
    end
    seg = vec(start+1:end);
    tseg = tvec(start+1:end) - tvec(start);
end
