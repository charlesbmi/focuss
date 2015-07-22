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

%   ============================================================================
%   Transform 32 coils into 4 compressed coils (using pre-computed transform,
%   see Buehrer et al., MRM 2007)
%   ============================================================================

dims    =   [64, 64, 1, 8];
func_data_512 = func_data_512(1:dims(1),1:dims(2),1:dims(4));
ncoils  =   1;
psens   =   coil_compression_xfm(1:ncoils, :)*reshape(sens2D,32,[]);
psens   =   reshape(psens, [ncoils, dims(1:2)]);

%   ============================================================================
%   Generate multi-coil nufft transform (and sampling, contained in the
%   k_radial_16.mat file) operator
%   ============================================================================
opts    =   {%'coils',   psens,...
             'k',       reshape([real(k(:)) imag(k(:))], [], dims(4),2),...
             'wi',      reshape(abs(k), [], dims(4))};

xfm     =   transform(dims, 'NUFFT', opts{:});

%   ============================================================================
%   Generate under-sampled multi-coil radial data. This is a bit of a short-cut,
%   since in reality we would have 32 coil measurements and would need to
%   transform the undersampled k-space data with the same coil compression xfm.
%   Here, we just shortcut this by using the compressed sensitivities directly
%   to generate the data.
%   ============================================================================
data    =   xfm*func_data_512;
