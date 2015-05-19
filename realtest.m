% NUFFT test on MRI data
close all;
%clear all;
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
load('meas_vismotor_z24.mat');
load('sens_z24.mat');
load('traj_vismotor.mat');
disp('Loaded');
fftw('planner','patient');
rng(0);

total_lines = 6000;
nc = size(z24,1);
z = squeeze(z24(1,:,1:total_lines));
[np nl] = size(z);
num_bin = 16;
nt = nl / num_bin;
nl = num_bin;

% convert sensitivties standard to recon size
% image size
sens = permute(sens_z24, [2 3 1]);
[nx,ny,nc] = size(sens);

% convert sensitivies to 


%k = repmat(k,[1 ceil(nl/8)]);
%wi = repmat(k,[1 ceil(nl/8)]);

tic
for t = nt:-1:1 % reverse index to pre-allocate st
    kx = k(:,(t*nl-nl+1):t*nl,1);
    ky = k(:,(t*nl-nl+1):t*nl,2);

    % create NUFFT structure
    N = [nx ny];
    J = [6 6];	% interpolation neighborhood
    K = N*2;	% two-times oversampling
    om = [kx(:) ky(:)];	% 'frequencies' are locations here!
    
    st(t) = nufft_init(om, N, J, K, N/2, 'minmax:kb'); % use 'table' to save memory

end
toc
disp('Pre-allocated')

wxy = wi(:,1:nt*nl);
% first should be points per time point. 
% second should be time points
wxy = reshape(wxy, [np nl nt]); % reshape into time lines?

% A = @(img,~) nufft3(ifft(img,[],3), st); % z is kt, A(xf) = kt
% AT = @(ks,~) fft(nufft3_adj(wxy.*ks,st),[],3); % A(kt) = xf

A = @(img,~) nufft3(img, st, nl); % z is kt, A(xf) = kt
AT = @(ks,~) nufft3_adj(wxy.*ks,st); % A(kt) = xf

% call the *adjoint* NUFFT
zr = reshape(z, [np nl nt]);
patterns = AT(zr,0);
pattern = patterns(:,:,5);

figure(1)
imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('pattern')

% Reconstruct all coils together
A = @(img,~) mcnufft3(img, st, nl,sens); % z is kt, A(xf) = kt
AT = @(ks,~) mcnufft3_adj(wxy.*ks,st,sens); % A(kt) = xf
%Y = reshape(z24(,:,1:total_lines), [np nl nt]);
recon = KTFOCUSS(A,AT,Y,Y,ones(size(Y)),factor,lambda,Minner,Mouter);

% Recontstruct each coil separately before combining
% for coil = nc:-1:1
% 	tic
% 	Y = reshape(z24(coil,:,1:total_lines), [np nl nt]);
% 	Low_Y = Y;
% 	mask = ones(size(Y));
% 	Mouter = 2;
% 	Minner = 40;
% 	factor = 0.5;
% 	lambda = 0.1;
% 	X_FOCUSS(:,:,:,coil) = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda, Minner, Mouter);
% 	toc % time focuss
% end

% figure(2);
% imagesc(abs(X_FOCUSS(:,:,1,1))); axis off; axis equal; colormap gray; colorbar;
% title('reference reconstruction');

% recon = multicoil_recon(X_FOCUSS,sens);

figure(3);
imagesc(abs(recon(:,:,1))); axis off; axis equal; colormap gray; colorbar;
title('all-coil reconstruction');
