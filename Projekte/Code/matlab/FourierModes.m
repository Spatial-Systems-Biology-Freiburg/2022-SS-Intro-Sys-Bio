function dpq = FourierModes(xmax, ymax)
dpq = zeros(1,xmax*ymax);
idx = 0;
% loop over all pairs of p < xmax,q < ymax to determine the possible modes:
for p = 1:xmax
    for q = 1:ymax
        idx = idx +1;
        dpq(idx) = (sin(pi*p/xmax)^2 + sin(pi*q/ymax)^2);
    end
end
dpq = unique(dpq);
end