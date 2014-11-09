function res = ifftnc(x)
% res = ifftnc(x)
% orthonormal reverse N-D FFT
% Charles Guan (charles.guan@stanford.edu)
res = sqrt(length(x(:)))*ifftshift(ifftn(fftshift(x)));
