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

% number of lines other than the central lines to sample
num_high_freq = round(ny / R) - num_low_freq;
nlf = num_low_freq;
nhf = num_high_freq;

assert (nhf >= 0, 'too many low frequency lines sampled');


mask = zeros(nx,ny,nt);

rand_vals = rand(ny,1,nt);
if pat == 0 % uniform random in ky lines
    for t = 1:nt
        rand_inds = randsample(ny-nlf, nhf)+nlf/2;
        mask(rand_inds,:,t) = 1;
    end
    mask(1:nlf/2,:,:) = 1;
    mask(end-nlf/2+1:end,:,:)=1;
elseif pat == 1
    mask = Random_DownsamplingMASK(nx,ny,nt,R,nlf); % can fftshift this to center if necessary
else % uniform random across kx, ky
    mask = rand(ny,nx,nt) < 1/R;
end
