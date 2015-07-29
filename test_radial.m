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
load(filename, 'func_data');

maskname = ['radial_sampling_masks.mat'];
disp(['Loading masks from: ', maskname]);
load(maskname);

full_sample_img = func_data; % func_data = data + mask
orig_img = data;
fmrib_img = mask; % fmrib image
disp('Loaded');

filenames = {'recon_results/radial_2lines_focuss.mat', 'recon_results/radial_4lines_focuss.mat', 'recon_results/radial_8lines_focuss.mat', 'recon_results/radial_16lines_focuss.mat', 'recon_results/radial_32lines_focuss.mat'};
ks = {radial_sampling_mask_2readout_lines_k, radial_sampling_mask_4readout_lines_k, radial_sampling_mask_8readout_lines_k, radial_sampling_mask_16readout_lines_k, radial_sampling_mask_32readout_lines_k};
wis = {radial_sampling_mask_2readout_lines_wi, radial_sampling_mask_4readout_lines_wi, radial_sampling_mask_8readout_lines_wi, radial_sampling_mask_16readout_lines_wi, radial_sampling_mask_32readout_lines_wi};
nls = [2, 4, 8, 16, 32];

% Optimize FFT
fftw('planner','patient');

[nx ny nt] = size(full_sample_img);
np = size(ks{1},1);

for idx = 3:length(filenames)
    filename = filenames{idx};
    k = ks{idx};
    wi = wis{idx};
    nl = nls(idx);
    display(filenames{idx});
    mask = masks{idx};
    kt_data_ds = ifft(kt_data.*mask,[],2);

tic
nx = dims(1); ny = dims(2); nt = dims(4);
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

wxy = wi(:,1:nt*nl);
wxy = reshape(wxy, [np nl nt]);

F = @(img,~) fft(nufft3(img, st, nl),[],3); % z is kt, A(xf) = kt
FT = @(ks,~) nufft3_adj(wxy.*ifft(ks,[],3),st); % A(kt) = xf
radial_data = nufft3(func_data, st, nl);

    % % function setting
    tic
    focuss_ft_recon = ifft(kt_focuss_old(F,FT,radial_data,ones(size(radial_data))),[],3);
    toc % time focuss
    recon = focuss_ft_recon;
    save(filenames{idx}, 'recon');

    figure;
    im([focuss_ft_recon, full_sample_img]);
    title('left: FFT recon; right: fully sampled original')

    x = 47; y = 20;
    figure;
    hold on;
    plot(abs(squeeze(focuss_ft_recon(y,x,:))));
    plot((squeeze(full_sample_img(y,x,:))),'r');
    legend('FT recon','KLT recon', 'mask or func_data');
    hold off;
    drawnow;
end
