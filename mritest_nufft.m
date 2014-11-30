% NUFFT test on MRI data
close all;
clear all;
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
load('2D_data.mat');
load('k_radial.mat');
disp('Loaded');

img = data(:,:,1:20);
[nx ny nt] = size(img);
nl = 8; % radial lines per time point
for t = 1:nt
    kxy = k(:,(t*nl-nl+1):t*nl);
    xc = real(kxy(:));
    yc = imag(kxy(:));

    % create NUFFT structure
    N = [nx ny];
    J = [5 5];	% interpolation neighborhood
    K = N*2;	% two-times oversampling
    om = [xc yc];	% 'frequencies' are locations here!
    
    st(t) = nufft_init(om, N, J, K, N/2, 'minmax:kb'); % use 'table' to save memory

end

wxy = wi(:,1:nt*nl);
wxy = reshape(wxy, [size(wxy,1)*nl size(wxy,2)/nl]); % reshape into time lines?

A = @(img) nufft3(img, st);
AT = @(ks) nufft3_adj(wxy.*ks,st);

% call the *adjoint* NUFFT
weights = A(img);
patterns = AT(weights);
pattern = patterns(:,:,5);

figure(1)
imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('pattern')
