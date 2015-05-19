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
%filename = ['fmrib_blur.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]); % (20,47,:) for F, (25,45,:) for M, (31,44,:) for R, (37,44,:) for I, (39,44,:) for B
load(filename);
full_sample_img = func_data; % func_data = data + mask
%full_sample_img = mask;
orig_img = data;
fmrib_img = mask; % fmrib image
disp('Loaded');

% Optimize FFT
fftw('planner','patient');

%% Downsampling rate 
ds_rate = 2;
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

[nx ny nt] = size(full_sample_img);
kt_data = fft(fft(full_sample_img,[],1),[],2);

mask = downsample_mask(nx,ny,nt,ds_rate,num_low_freq,1);
kt_data_ds = ifft(kt_data.*mask,[],2);

ref = ifft(kt_data_ds,[],1);
%ref = ifft(ifft(kt_data_ds(:,:,1),[],1),[],2);
%figure(5);
%imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
%title('reference undersampled');

% % function setting
F = @(x) ifft(fft(x,[],1),[],3);
FT = @(x) ifft(fft(x,[],3),[],1);
tic
focuss_ft_recon = ifft(kt_focuss_old(F,FT,kt_data_ds,mask,num_low_freq),[],3);
toc % time focuss

% KT-FOCUSS using rho as KLT support basis.
v = ifft(focuss_ft_recon,[],3);
% principal components were calculated from low-resolution x–f support obtained using only
% low-frequency k-space samples.
[~,~,V] = svd(reshape(v,size(v,1)*size(v,2),size(v,3)));
A = @(xklt,mask) iklt3(fft(xklt,[],1),V); % x ky klt
AT = @(kt,mask) ifft(klt3(kt,V),[],1);

% A = @(xklt,mask) fft(fft(iklt3(xklt,V),[],1),[],2); % x ky klt
% AT = @(kt,mask) klt3(ifft(ifft(kt,[],1),[],2),V);

tic
%focuss_klt_recon = iklt3(kt_focuss_old(A,AT,fft(focuss_ft_recon,[],1),mask,num_low_freq),V);
toc

%imagesc(abs(focuss_ft_recon(:,:,2))); axis off; axis equal; colormap gray; colorbar;
%title('reference reconstruction');
%figure
err = norm(full_sample_img(:) - focuss_ft_recon(:))
em = err_map(focuss_ft_recon, full_sample_img);
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

%ets = err_plot(focuss_ft_recon, full_sample_img);
%figure;
%plot(ets);

figure;
im([focuss_ft_recon, full_sample_img]);
title('left: FFT recon; right: fully sampled original')

x = 47; y = 20;
figure;
hold on;
plot(abs(squeeze(focuss_ft_recon(y,x,:))));
plot((squeeze(full_sample_img(y,x,:))),'r');
legend('recon','mask or func_data');
title('signal ')
hold off;