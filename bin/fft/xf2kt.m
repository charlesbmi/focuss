function kt = xf2kt(xf)
% Transforms zero-centered x-f space data into zero-centered k-t space data
% Charles Guan (charles.guan@stanford.edu)
if ndims(xf) == 2
    kt = fftc(ifftc(xf,2),1);
elseif ndims(xf) == 3
    kt = fftc(fftc(ifftc(xf,3),2),1);
end
