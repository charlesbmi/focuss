function xf = kt2xf(kt)
% Transforms zero-centered k-t space data to zero-centered x-f space data
% Charles Guan (charles.guan@stanford.edu)
if ndims(kt) == 2
    xf = fftc(ifftc(kt,1),2);
elseif ndims(kt) == 3
    xf = fftc(ifftc(ifftc(kt,1),2),3);
end
