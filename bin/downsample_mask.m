function mask = downsample_mask(nx, ny, nt, R, num_low_freq, pat)
% Creates downsampling mask k-t data in one dimension given an acceleration factor
% Inputs
%       nx, ny, nt - dimensions of mask
%       R - acceleration (undersampling) factor
%       num_low_freq - number of lines above and below the center to fully sample
%       pat - undersampling pattern. options: 
%           0 - uniform random in ky lines
%           1 - gaussian random in ky lines
%           2 - uniform random across kx, ky
% Output
%       du - undersampled data
%       mask - sampling locations
%
% Charles Guan (charles.guan@stanford.edu)

assert(mod(num_low_freq,2)==0, 'num_low_freq must be even')

rand_vals = rand(nx,1,nt);
if pat == 0 % uniform random in ky lines
    mask = rand_vals < 1/R;
    mask = repmat(mask,[1,ny,1]);
elseif pat == 1
    mask = Random_DownsamplingMASK(nx,ny,nt,R,0); % can fftshift this to center if necessary
else % uniform random across kx, ky
    mask = rand(nx,ny,nt) < 1/R;
end
% low-frequency full sampling
mask(1:num_low_freq/2,:,:) = 1;
mask(end-num_low_freq/2+1:end,:,:)=1;
