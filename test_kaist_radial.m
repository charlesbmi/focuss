%   Example 2D multi-coil radial undersampled k-t FOCUSS reconstruction
%
%   Charles Guanw
%

%   ============================================================================
%   Load data and set paths
%   ============================================================================
addpath(genpath('bin'));
addpath(genpath('data'));
im = @(x) imshow(mat2gray(abs(x(:,:,1))));

load('sens_data');
load('2D_data.mat','func_data');
func_data_512 = func_data;
load('radial_sampling_masks.mat');

filenames = {'recon_results/radial_2lines_kaist.mat', 'recon_results/radial_4lines_kaist.mat', 'recon_results/radial_8lines_kaist.mat', 'recon_results/radial_16lines_kaist.mat', 'recon_results/radial_32lines_kaist.mat'};
ks = {radial_sampling_mask_2readout_lines_k, radial_sampling_mask_4readout_lines_k, radial_sampling_mask_8readout_lines_k, radial_sampling_mask_16readout_lines_k, radial_sampling_mask_32readout_lines_k};
wis = {radial_sampling_mask_2readout_lines_wi, radial_sampling_mask_4readout_lines_wi, radial_sampling_mask_8readout_lines_wi, radial_sampling_mask_16readout_lines_wi, radial_sampling_mask_32readout_lines_wi};
nls = [2, 4, 8, 16, 32];

%   ============================================================================
%   Transform 32 coils into 4 compressed coils (using pre-computed transform,
%   see Buehrer et al., MRM 2007)
%   ============================================================================

%for idx = 3:length(filenames)
for idx = [3,4]
  filename = filenames{idx};
  k = ks{idx};
  wi = wis{idx};
  nl = nls(idx);

% Generate data.
dims    =   [64, 64, 1, 512];
func_data_512 = func_data_512(1:dims(1),1:dims(2),1:dims(4));
ncoils  =   1;
psens   =   coil_compression_xfm(1:ncoils, :)*reshape(sens2D,32,[]);
psens   =   reshape(psens, [ncoils, dims(1:2)]);

%   ============================================================================
%   Generate multi-coil nufft transform (and sampling, contained in the
%   k_radial_16.mat file) operator
%   ============================================================================

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
toc
disp('Pre-allocated')

np = size(k,1);
wxy = wi(:,1:nt*nl);
% first should be points per time point. 
% second should be time points
wxy = reshape(wxy, [np nl nt]); % reshape into time lines?

A = @(img,~) nufft3(img, st, nl); % z is kt, A(xf) = kt
AT = @(ks,~) nufft3_adj(wxy.*ks,st); % A(kt) = xf

% call the *adjoint* NUFFT
radial_data = A(func_data,0);
patterns = AT(radial_data,0);
pattern = patterns(:,:,5);

% figure(1)
% imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('pattern')

% FOCUSS parameters
Mouter = 3;
Minner = 40;
factor = 0.5;
lambda = 0.1;

    tic
    Y = reshape(radial_data, [np nl nt]);
    LOW_Y = Y;
    mask = ones(size(Y));
    X_FOCUSS = KTFOCUSS(A,AT,Y,LOW_Y,mask,factor,lambda, Minner, Mouter);
    toc

    err = norm(func_data(:) - X_FOCUSS(:))

    %figure;
    %im([X_FOCUSS, func_data]);
    %title('left: FFT recon; right: fully sampled original')

    recon = X_FOCUSS;
    %figure;
    %hold on;
    %plot(abs(squeeze(fmrib_recon(20,47,:))));
    %plot((squeeze(func_data(20,47,:))),'r');
    %legend('recon','mask or func_data');
    %hold off;
    %drawnow;

    save(filenames{idx}, 'recon');
end
