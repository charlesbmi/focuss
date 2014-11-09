% k-t FOCUSS
% Based off Jung et al, 2007

%% Add path  recursively
clear all;
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
filename = ['2D_data.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]);
load(filename);
full_sample_img = func_data; % func_data = data + mask
orig_img = data;
fmrib_img = mask; % fmrib image
disp('Loaded');

%% Downsampling rate 
ds_rate = 4; % option to change
ds_pat = 1; % downsampling pattern, 0 for uniform random in ky, 1 for gaussian in ky, 2 for random in kx and ky % option to change
switch ds_pat
case 0
    ds_pat_str = 'uniform random in ky';
case 1
    ds_pat_str = 'gaussian random in ky';
otherwise
    ds_pat_str = 'uniform random in kx and ky';
end
disp(['Downsample rate: ',num2str(ds_rate)]);
disp(['Downsample pattern: ',ds_pat_str]);

% low frequency full sampling number (1~num_phase/2), (end-num_phase/2~end) -> full sampling
num_low_freq = 4; % option to change
disp(['Low frequency full sampling number: ', num2str(num_low_freq)]);

% Optimize FFT
fftw('planner','patient');

[nx ny nt] = size(full_sample_img);
kt_data = fft(fft(full_sample_img,[],1),[],2);

mask = downsample_mask(nx,ny,nt,ds_rate,num_low_freq,1);
kt_data_ds = kt_data.*mask;

ref = ifft(ifft(kt_data_ds(:,:,1),[],1),[],2);
figure(5);
imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
title('reference undersampled');

% % function setting
A = @(x,mask)  fft(fft(x,[],1),[],2).*mask;
AT = @(x,mask) ifft(ifft(x.*mask,[],1),[],2);

tic
X_FOCUSS = kt_focuss(A,AT,kt_data_ds,mask,num_low_freq);
toc % time focuss

err = norm(full_sample_img(:) - X_FOCUSS(:))
em = err_map(X_FOCUSS, full_sample_img);
imshow(mat2gray(em));
title('error map')
ets = err_plot(X_FOCUSS, full_sample_img);
figure;
plot(ets);

% save for further analysis
save results.mat ets em err ds_rate ds_pat num_low_freq ref
