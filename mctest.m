% k-t FOCUSS
% Based off Jung et al, 2007

%% Add path  recursively
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
%filename = '2D_data.mat'; % load full x-y-t data, and coils
filename = 'old_test_data.mat'; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]);
if exist('data_2D','var') == 0
    load(filename);
end
data = data_2D;
filename = '8coil_map.mat'; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]);
if exist('ecmap_3D','var') == 0
    load(filename);
end
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

[nx ny nt] = size(data);
nc = size(ecmap_3D,4);
coil_data = repmat(data, [1 1 1 nc]).*ecmap_3D; % coil-weighted observations
kt_data = fft(fft(coil_data,[],1),[],2);

mask = downsample_mask(nx,ny,nt,ds_rate,num_low_freq,1);
kt_data_ds = kt_data.*repmat(mask, [1 1 1 nc]); % replicate mask across coils

ref = ifft(ifft(kt_data_ds(:,:,1,1),[],1),[],2);
figure(5);
imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
title('reference');

% % function setting
A = @(x,mask)  fft(fft(x,[],1),[],2).*mask;
AT = @(x,mask) ifft(ifft(x.*mask,[],1),[],2);

tic
X_FOCUSS = mc_kt_focuss(A,AT,kt_data_ds,mask,num_low_freq,ecmap_3D);
toc % time focuss

err = norm(data(:) - X_FOCUSS(:))
