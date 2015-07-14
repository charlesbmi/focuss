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

%% Downsampling rate 
% Optimize FFT
fftw('planner','patient');

kt_data = fft(fft(full_sample_img,[],1),[],2);

masks = {cart_sampling_mask_2x_4low_freq, cart_sampling_mask_4x_4low_freq, cart_sampling_mask_8x_4low_freq, cart_sampling_mask_12x_4low_freq, cart_sampling_mask_16x_4low_freq};
filenames = {'recon_results/cart_2x_4lowfreq_kaist.mat', 'recon_results/cart_4x_4lowfreq_kaist.mat', 'recon_results/cart_8x_4lowfreq_kaist.mat', 'recon_results/cart_12x_4lowfreq_kaist.mat', 'recon_results/cart_16x_4lowfreq_kaist.mat'};
for idx = 1:length(masks)

    mask = masks{idx};

    kt_data_ds = kt_data.*mask;

    % % function setting
    F = @(x,mask) fft(x,[],1).*mask;
    FT = @(x,mask) ifft(x.*mask,[],1);
    tic
    X_FOCUSS = ifft(kt_focuss(F,FT,kt_data_ds,mask,4),[],2);
    toc % time focuss

    err = norm(full_sample_img(:) - X_FOCUSS(:))
    em = err_map(X_FOCUSS, full_sample_img);
    norm_av = err_map(full_sample_img,0); % todo instead of dividing by norm map, just divide by total norm / 
    norm_all = norm(full_sample_img(:) / prod(size(em))); % or perhaps mean of mean of norm. to get norm per pixel

    figure;
    im([X_FOCUSS, full_sample_img]);
    title('left: FFT recon; right: fully sampled original')

    fmrib_recon = X_FOCUSS;
    figure;
    hold on;
    plot(abs(squeeze(fmrib_recon(20,47,:))));
    plot((squeeze(full_sample_img(20,47,:))),'r');
    legend('recon','mask or func_data');
    hold off;

    save(filenames{idx}, 'fmrib_recon');
end
