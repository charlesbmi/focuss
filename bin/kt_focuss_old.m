function [ rho ] = kt_focuss(F, FT, vu, mask, nlf)
% Attempts to solve vu = F*rho for rho, with a soft sparsity constraint on rho,
% using the k-t FOCUSS image reconstruction algorithm for kx-t and kx-ky-t data
% k-t FOCUSS method taken from Jung et al. (2009) and Zong (2014):
% (0) Initialize parameters and image estimate
% (1) Initialize rho to zero-filled inverse of v for initial estimate
% (2) Repeat until converged
% (2a) Compute weighting matrix W
% (2b) Compute the nth FOCUSS estimate by solving l2 optimization
% Inputs
%       F - Sparsifying operator, such as 3D Fourier transform
%       FT - Transpose operator of F
%       vu - raw undersampled x-ky-t space image, to be pruned to v
%       mask - sampling locations of vu, 1 where sampled, 0 where
%           missing data
% Outputs
%       v - v(k,t), the converged k-t image estimate
%
% Charles Guan (charles.guan@stanford.edu)

% todo: make these flexible parameters
disp('k-t FOCUSS start:')
L = 1e-2;
p = 0.5;
updates = 2;

% Optimize FFT
fftw('planner','patient');

% Main k-t FOCUSS loop
rho_initial = FT(vu);
rho = rho_initial;
for iter = 1:updates
    disp(sprintf('kt-FOCUSS iteration #: %d',iter))
    W = abs(rho).^p;
    % Choose between using Nielsen conj_grad and MATLAB pcg
    rho = focuss_est(vu,mask,W,L,F,FT);
end
baseline_flux_norm = norm(rho(:)) % TODO norm should be frobenius norm, not norm(x,2)
diff_eh = abs(rho) - abs(rho_initial);
diff_test = norm(diff_eh(:))
