% Non-unform k-t FOCUSS test
% Based off Jung et al, 2007

%% Add path  recursively
clear all;
addpath(genpath('bin'));
addpath(genpath('data'));
% take undersampling, by eg taking every other line and undersampling by just using less lines so you have less lines per bin but still 512 lines

%% Load full measurement 
load('2D_data.mat'); % load full x-y-t data, and coils
load('k_radial.mat'); % load k-space sampling locations
disp('Loaded');

nl = 32; % number of sampling lines. or you can say R = 8; compare nl = 64 / 8 and Cartesian
% each sampling line contains 
% todo for more sampling, we can extend k
ppl = size(k,1); % k-space sampling points per line

img = mask(:,:,1:512); % todo starting small
[nx ny nt] = size(img);

k = repmat(k,[1 ceil(nl/8)]);
wi = repmat(k,[1 ceil(nl/8)]);

tic
%parfor t = 1:nt
for t = min(20,nt):-1:1 % go backwards to pre-allocate st
    kxy = k(:,(t*nl-nl+1):t*nl);
    kx = real(kxy(:)); % kx frequencies of each sampling point
    ky = imag(kxy(:)); % ky frequencies of each sampling point

    % create NUFFT structure
    N = [nx ny];
    J = [6 6];	% interpolation neighborhood
    K = N*2;	% two-times oversampling
    om = [kx ky];	% frequencies are locations here!
    
    st(t) = nufft_init(om, N, J, K, N/2, 'minmax:kb'); % use 'table' to save memory
end
disp('structures initialized')
toc

wxy = wi(:,1:nt*nl);
wxy = reshape(wxy, [ppl*nl nt]); % reshape into time lines?

A = @(img,~) nufft3_20(img, st);
AT = @(ks,~) nufft3_20_adj(wxy.*ks,st);
ATA = @(img,~) AT(A(img));

% TODO why does adj size go to 64 64 20 and not 1024 20?

% undersample radially by NUFFT to undersampled radial lines and NUFFT_adj back for each line
patterns = ATA(img);
imshow(mat2gray(abs(patterns(:,:,4))))

tic
%X_FOCUSS = kt_focuss(A,AT,A(img),ones(size(A(img))),4);
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda = 0.1;
%X_FOCUSS = KTFOCUSS(A,AT,patterns,patterns,ones(size(patterns)),factor,lambda, Minner, Mouter);
% todo maybe it is just better to take that data and then apply Cartesian... I don't think NUFFT is the sparsifying transform we want is it? Does it sparsify? Maybe in time... in k-space it takes the near-zero entries.
toc % time focuss

%imagesc(X_FOCUSS(:,:,5))
%axis off; axis equal; colormap gray; colorbar; title('recon')

% calculate mask if you do it the cartesian way. But then you lose to inaccuracy of assigning it a discrete frequency. 
% oversample?
mask = zeros(size(img));
for t = nt:-1:1
    linearInd = sub2ind(size(img),round(st(t).om*31/pi)+32);
    mask(linearInd,t) = 1; % todo about this this is 3d, so may want to build this up or get all indexing
end
