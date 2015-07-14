 % k-t FOCUSS
% Based off Jung et al, 2007

%% Add path  recursively
%clear all;
close all;
addpath(genpath('bin'));
addpath(genpath('data'));
im = @(x) imshow(mat2gray(abs(x(:,:,1))));

%% Load full measurement 
filename = ['2D_data.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]); % (20,47,:) for F, (25,45,:) for M, (31,44,:) for R, (37,44,:) for I, (39,44,:) for B
load(filename);

maskname = ['sampling_masks.mat'];
disp(['Loading masks from: ', maskname]);
load(maskname);

full_sample_img = func_data; % func_data = data + mask
orig_img = data;
fmrib_img = mask; % fmrib image
disp('Loaded');

% Optimize FFT
fftw('planner','patient');

kt_data = fft(fft(full_sample_img,[],1),[],2);

masks = {cart_sampling_mask_2x_4low_freq, cart_sampling_mask_4x_4low_freq, cart_sampling_mask_8x_4low_freq};
filenames = {'recon_results/cart_2x_4lowfreq_focuss.mat', 'recon_results/cart_4x_4lowfreq_focuss.mat', 'recon_results/cart_8x_4lowfreq_focuss.mat'};

[nx ny nt] = size(full_sample_img);
kt_data = fft(fft(full_sample_img,[],1),[],2);

idx = 3;
mask = masks{idx};
kt_data_ds = ifft(kt_data.*mask,[],2);

% % function setting
F = @(x) ifft(fft(x,[],1),[],3).*mask;
FT = @(x) ifft(fft(x.*mask,[],3),[],1);
tic
focuss_ft_recon = ifft(kt_focuss_old(F,FT,kt_data_ds,mask),[],3);
toc % time focuss

% KT-FOCUSS using rho as KLT support basis.
xtd = ifft(kt_data_ds,[],1);
% principal components were calculated from low-resolution x–f support obtained using only
% low-frequency k-space samples.
[~,~,V] = svd(reshape(v,size(v,1)*size(v,2),size(v,3)));
A = @(xklt,mask) iklt3(fft(xklt,[],1),V); % x ky klt
AT = @(kt,mask) ifft(klt3(kt,V),[],1);

a = @(xklt,mask) fft(fft(iklt3(xklt,v),[],1),[],2); % x ky klt
at = @(kt,mask) klt3(ifft(ifft(kt,[],1),[],2),v);

tic
%focuss_klt_recon = iklt3(kt_focuss_old(A,AT,fft(fft(focuss_ft_recon,[],1),[],2),mask,
toc

figure;
im([focuss_ft_recon, full_sample_img]);
title('left: FFT recon; right: fully sampled original')

x = 47; y = 20;
figure;
hold on;
plot(abs(squeeze(focuss_ft_recon(y,x,:))));
plot((squeeze(full_sample_img(y,x,:))),'r');
legend('recon','mask or func_data');
hold off;
