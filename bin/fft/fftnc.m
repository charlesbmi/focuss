function res = fftnc(x)
% res = fftnc(x)
% orthonormal forward N-D FFT
% Charles Guan (charles.guan@stanford.edu)
res = 1/sqrt(length(x(:)))*fftshift(fftn(ifftshift(x)));
