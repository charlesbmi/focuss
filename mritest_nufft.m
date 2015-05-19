% NUFFT test on MRI data
close all;
%clear all;
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
load('2D_data.mat');
load('k_radial.mat');
disp('Loaded');
rng(0);

%img = func_data(:,:,1:100);
img = func_data;
[nx ny nt] = size(img);
nl = 2; % radial lines per time point
%for nl = [4,32]
ppl = size(k,1); % k-space sampling points per line
mask = zeros(size(img)); % for pseudo-radial, utilize a mask
k = repmat(k,[1 ceil(nl/8)]);
wi = repmat(k,[1 ceil(nl/8)]);

xmax = max(real(k(:)));
xmin = min(real(k(:)));
xr = xmax-xmin;
ymax = max(imag(k(:)));
ymin = min(imag(k(:)));
yr = ymax-ymin;

for t = nt:-1:1 % reverse index to pre-allocate st
    kxy = k(:,(t*nl-nl+1):t*nl);
    xc = real(kxy(:));
    yc = imag(kxy(:));

    % create NUFFT structure
    N = [nx ny];
    J = [5 5];	% interpolation neighborhood
    K = N*2;	% two-times oversampling
    om = [xc yc];	% 'frequencies' are locations here!
    
    %st(t) = nufft_init(om, N, J, K, N/2, 'minmax:kb'); % use 'table' to save memory

    xc = round((xc-xmin)/xr*(nx-1)+1); % expand to fit 1-to-64 indexed k-space
    yc = round((yc-ymin)/yr*(ny-1)+1);

    tc = t*ones(size(xc));
    sampl_idx = sub2ind(size(img),xc,yc,tc);
    mask(sampl_idx) = 1;

end

wxy = wi(:,1:nt*nl);
wxy = reshape(wxy, [ppl*nl nt]); % reshape into time lines?

A = @(img) nufft3(img, st);
AT = @(ks) nufft3_adj(wxy.*ks,st);

% call the *adjoint* NUFFT
%weights = A(img);
%patterns = AT(weights);
%pattern = patterns(:,:,5);

%figure(1)
%imagesc(abs(pattern)); axis off; axis equal; colormap gray; colorbar; title('pattern')

kt_data = fftc(fftc(img,1),2);
kt_data_ds = kt_data.*mask;

ref = ifftc(ifftc(kt_data_ds(:,:,1),1),2);
%figure(5);
%imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
title('reference undersampled');

A = @(x,mask)  fftc(fftc(x,1),2).*mask;
AT = @(x,mask) ifftc(ifftc(x.*mask,1),2);
%A = @(x,mask)  ifft(fftc(fftc(x,1),2).*mask,[],3);
%AT = @(x,mask) ifftc(ifftc(fft(x,[],3).*mask,1),2);

tic
Y = kt_data_ds;
Low_Y = Y;
Low_Y([1:32-2,33+2:end],:,:)=0;
Low_Y(:,[1:32-2,33+2:end],:)=0;
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda = 1;
X_FOCUSS = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda, Minner, Mouter);
%X_FOCUSS = kt_focuss(A,AT,kt_data_ds,mask,10);
toc % time focuss

%imagesc(abs(X_FOCUSS(:,:,1))); axis off; axis equal; colormap gray; colorbar;
title('reference reconstruction');

%figure(3)
plot(abs(squeeze(X_FOCUSS(20,47,:))))
hold on
plot(abs(squeeze(func_data(20,47,:))),'r')

err = norm(func_data(:) - X_FOCUSS(:))
em = err_map(X_FOCUSS, func_data);
%figure
%imagesc(em);
axis off; axis equal; colormap gray; colorbar; title('error map')

norm_av = err_map(func_data,0); % todo instead of dividing by norm map, just divide by total norm / 
%figure
%imagesc(em./norm_av);
axis off; axis equal; colormap gray; colorbar; title('error map normalized by norm of each voxel'); caxis([0 0.05]); % normalize caxis to 5%

%figure
norm_all = norm(func_data(:) / prod(size(em))); % or perhaps mean of mean of norm. to get norm per pixel
%imagesc(em/(norm_all));
axis off; axis equal; colormap gray; colorbar; title('error map normalized by average norm of entire'); caxis([0 0.05]); % normalize caxis to 5%

ets = err_plot(X_FOCUSS, func_data);
%figure;
plot(ets);

recon = X_FOCUSS;
save (sprintf('recon_results/radial_%dlines_lam1.mat', nl), 'recon','mask','k','wi');
%end
