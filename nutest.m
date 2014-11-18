% Non-unform k-t FOCUSS test
% Based off Jung et al, 2007

%% Add path  recursively
clear all;
addpath(genpath('bin'));
addpath(genpath('data'));
% take undersampling, by eg taking every other line and undersampling by just using less lines so you have less lines per bin but still 512 lines

%% Load full measurement 
filename = ['2D_data.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]); % (20,47,:) for F, (25,45,:) for M, (31,44,:) for R, (37,44,:) for I, (39,44,:) for B
load(filename);
filename = ['k_radial.mat']; % load k-space sampling locations
disp(['Loading data from: ',filename]); % (20,47,:) for F, (25,45,:) for M, (31,44,:) for R, (37,44,:) for I, (39,44,:) for B
load(filename);

[nx,ny,nt] = size(data);

R = 1; % either 2^0, 2^1, 2^2, or 2^3.
nl = 8 / R; % number of consecutive lines to bin in 1 time. even additions
[pts_per_line, num_total_lines] = size(k); % set of sequential radial spokes, ordered by an angle increment of 111.25º (pi/golden ratio).
pts_per_time = pts_per_line*nl;
%k = reshape(k, [pts_per_line*nl,num_total_lines/(nl*R)]); % each column is a radial trajectory and each row is a time sequence? not sure if that's right
kx = real(k);
ky = imag(k);
t = ceil([1:R:num_total_lines*pts_per_line]/pts_per_time).'; % bin (dividing and taking ceiling) consecutive lines into the same time

% create NUFFT structure. Start off 2D
N = size(data);
J = [6 6 6];	% interpolation neighborhood
K = N*2;	% two-times oversampling
om = [kx(:) ky(:) t];	% discrete k-space frequencies to sample on
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');

weights = nufft(img, st);
