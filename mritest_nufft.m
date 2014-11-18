% NUFFT test on MRI data
close all;
clear all;
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
load('2D_data.mat');
load('k_radial.mat');
disp('Loaded');

nl = 10;
kxy = k(:,1:nl);
wxy = wi(:,1:nl);
xc = real(kxy(:));
yc = imag(kxy(:));

img = data(:,:,1);

% create NUFFT structure
%N = [1 1]*2^8;
N = size(img);
J = [5 5];	% interpolation neighborhood
K = N*2;	% two-times oversampling
om = [xc yc];	% 'frequencies' are locations here!

% the following line probably needs to be changed
% to get the proper scaling/dimensions in the pattern
% but for now i just make it fill up [-pi/20,pi/20]
% in hopes of getting a nice 'picture' of pattern
%om = (pi/20) * om / max(om(:));
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');

%weights = ones(size(xc)); % equal weights on each element; could change
weights = nufft(img, st);

% call the *adjoint* NUFFT; this is what does "the FT of unequal data"
pattern = nufft_adj(weights, st);

figure(1)
imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('pattern')

pattern = nufft_adj(weights.*wxy(:), st);
figure(2)
imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('weights')
