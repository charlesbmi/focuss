 % k-t FOCUSS
% Based off Jung et al, 2007

%% Add path  recursively
%clear all;
close all;
addpath(genpath('bin'));
addpath(genpath('data'));
im = @(x) imshow(mat2gray(abs(x(:,:,1))));

%% Load full measurement 
rng(0);
filename = ['2D_data.mat']; % load full x-y-t data, and coils
% filename = ['fmrib_blur.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]); % (20,47,:) for F, (25,45,:) for M, (31,44,:) for R, (37,44,:) for I, (39,44,:) for B
load(filename);
full_sample_img = func_data; % func_data = data + mask
% full_sample_img = mask;
orig_img = data;
fmrib_img = mask; % fmrib image
disp('Loaded');

%% Downsampling rate 
ds_rate = 6;
ds_pat = 1; % downsampling pattern, 0 for uniform random in ky, 1 for gaussian in ky, 2 for random in kx and ky % option to change

switch ds_pat
case 0
    ds_pat_str = 'cart';
case 1
    ds_pat_str = 'gauss';
otherwise
    ds_pat_str = 'uniform random in kx and ky';
end
disp(['Downsample rate: ',num2str(ds_rate)]);
disp(['Downsample pattern: ',ds_pat_str, ' random in ky']);

% low frequency full sampling number (1~num_phase/2), (end-num_phase/2~end) -> full sampling
num_low_freq = 2; % option to change
disp(['Low frequency full sampling number: ', num2str(num_low_freq)]);

% Optimize FFT
fftw('planner','patient');

[nx ny nt] = size(full_sample_img);
kt_data = fft(fft(full_sample_img,[],1),[],2);

mask = downsample_mask(nx,ny,nt,ds_rate,num_low_freq,1);
kt_data_ds = kt_data.*mask;

ref = ifft(ifft(kt_data_ds,[],1),[],2);
%ref = ifft(ifft(kt_data_ds(:,:,1),[],1),[],2);
%figure(5);
%imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
%title('reference undersampled');

% % function setting
F = @(x,mask) fft(x,[],1).*mask;
FT = @(x,mask) ifft(x.*mask,[],1);
tic
X_FOCUSS = ifft(kt_focuss(F,FT,kt_data_ds,mask,num_low_freq),[],2);
toc % time focuss

xt_data_ds = ifft(ifft(kt_data_ds,[],1),[],2);

% KLT with FFT first
xe = X_FOCUSS;
kt_data_ds = fft(fft(X_FOCUSS,[],1),[],2);

% Perhaps we should take svd of kx-y-t space. 
% KLT PCA FOCUSS
% This is correct. the svd. xyklt 
[~,~,V] = svd(reshape(xe,size(xe,1)*size(xe,2),size(xe,3)));
% Original
%F = @(xyklt,mask) ifft(fft(fft(iklt3(xyklt,V),[],1),[],2),[],3).*mask; % v = F (rho), rho is xyklt, v is xkyt
%FT = @(kt,mask) klt3(ifft(ifft(fft(kt.*mask,[],3),[],1),[],2),V);  % rho = FT (v)
% flipped
%A = @(xyklt,mask) ifft(fft(fft(iklt3(xyklt,V),[],1),[],2),[],3).*mask; % v = F (rho), rho is xyklt, v is xkyt
%AT = @(kt,mask) klt3(ifft(ifft(fft(kt.*mask,[],3),[],1),[],2),V);  % rho = FT (v)

A = @(xyklt,mask) fft(fft(iklt3(ifft(xyklt,[],3),V),[],1),[],2); % v = F (rho), rho is xyklt, v is xkyt
AT = @(kt,mask) fft(klt3(ifft(ifft(kt,[],1),[],2),V),[],3);  % rho = FT (v)

tic
%X_FOCUSS = kt_focuss(A,AT,kt_data_ds,ones(size(mask)),64); % will give back ifft(rho)? or A(. No, returns FT (kt) == rho, without the fft
toc


% KT-FOCUSS using rho as KLT support basis. but rho should stay in kx mode. It's not needed to do it the other way.
%v = ifft(kt_data_ds,[],1); % basis is now in x-ky-t
%[~,~,V] = svd(reshape(v,size(v,1)*size(v,2),size(v,3)));
%A = @(kxyklt,mask) ifft(fft(iklt3(kxyklt,V),[],1),[],3); % x ky klt
%AT = @(kt,mask) klt3(ifft(fft(kt,[],3),[],1),V);

%X_FOCUSS = ifft(kt_focuss(A,AT,kt_data_ds,mask,num_low_freq),[],2);

%imagesc(abs(X_FOCUSS(:,:,2))); axis off; axis equal; colormap gray; colorbar;
%title('reference reconstruction');
%figure
err = norm(full_sample_img(:) - X_FOCUSS(:))
em = err_map(X_FOCUSS, full_sample_img);
%imagesc(em);
%axis off; axis equal; colormap gray; colorbar; title('error map')

norm_av = err_map(full_sample_img,0); % todo instead of dividing by norm map, just divide by total norm / 
%figure
%imagesc(em./norm_av);
%axis off; axis equal; colormap gray; colorbar; title('error map normalized by norm of each voxel'); caxis([0 0.05]); % normalize caxis to 5%

%figure
norm_all = norm(full_sample_img(:) / prod(size(em))); % or perhaps mean of mean of norm. to get norm per pixel
%imagesc(em/(norm_all));
%axis off; axis equal; colormap gray; colorbar; title('error map normalized by average norm of entire'); caxis([0 0.05]); % normalize caxis to 5%

%ets = err_plot(X_FOCUSS, full_sample_img);
%figure;
%plot(ets);

figure;
im([xe, full_sample_img]);
title('left: FFT recon; right: fully sampled original')

figure;
hold on;
plot(abs(squeeze(X_FOCUSS(20,47,:))));
plot((squeeze(full_sample_img(20,47,:))),'r');
legend('recon','mask or func_data');
hold off;

%xx = [];
%figure;
%for t = 1:512
    %xx(t) = norm(X_FOCUSS(:,:,t),'fro');
%end
%plot(xx);

%save results.mat ets em err ds_rate ds_pat num_low_freq ref
%save (sprintf('recon_results/%s_%dx_%dlowfreq_klt.mat', ds_pat_str, ds_rate, num_low_freq), 'recon','mask');

% the F in fMRIB is at (20,47,:)
