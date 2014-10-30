%function [ v ] = kt_focuss(vu, mask)
% In progress: Implements the k-t FOCUSS image reconstruction algorithm 
% k-t FOCUSS method taken from Jung et al. (2009) and Zong (2014):
% (0) Initialize parameters and image estimate
% (1) Initialize rho to zero-filled inverse of v for initial estimate
% (2) Repeat until converged
% (2a) Compute weighting matrix W
% (2b) Compute the nth FOCUSS estimate by solving l2 optimization
% Inputs
%       vu - raw undersampled (kx,ky,t) space image
%       mask - sampling locations of vu, 1 where sampled, 0 where
%           missing data
%       ve - initial estimate of k-t space image
%       opt - sparsifying transform option, 'ft' or 'klt'
% Outputs
%       v - v(k,t), the converged k-t image estimate

% Idea: we wish to minimize
% MF multibojective function min ||q||_2 (frobenius norm) + lambda ||v-FWq||
%   i.e. solve th
% vs rho = rho(x,y,f)
% Assumes p - 0.5, such that rho = Wq = q.*q
% uses all 

clear;
disp('Initializing data');
load('test_data/2D_data.mat');
%load('test_data/mask.mat');
addpath('fft'); % Lustig FFT library
addpath('optimization_tools'); % Stanford Systems Optimization Lab least squares

%v = fftc(fftc(data,2),1);
v = fft2c(data);
mask = repmat(rand([size(v,1),1,size(v,3)]) < 0.5,[1,size(v,2),1]);
vu = v.*mask;

disp('Initializing k-t FOCUSS')
L = 1e-3;   % lambda weighting multiobjective
updates = 1;
p = 0.5;    % from paper
% Optimize FFT
fftw('planner','patient');

% Resolve function handles
disp(sprintf('sparsifying transform selected: %s','Fourier'))
%ft = @(xf) ifft(fft(xf,[],1),[],2); % uncentered
%ift = @(kt) fft(ifft(xf,[],1),[],2); % centered

%q = kt2xf(vu); % initial estimate
q = fft2c(vu);
q = q./sqrt(abs(q));
for iter = 1:updates
    prev_q = q;
    disp(sprintf('kt-FOCUSS iteration #: %d',iter))
    tic
    %AFUN = @(x,topt) afun(x,prev_q,size(vu),mask,L,topt);
    % AFUN(x,'notransp') should give A*x. AFUN(x,'transp') should give A'*x)
    %q = lsqr(AFUN, [vu(logical(mask)); zeros(numel(q),1)],[],5,[],[],prev_q(:));
    AFUN = @(x,topt) lfun(x,prev_q,size(vu),mask,topt);
    % AFUN(x,1) should give A*x. AFUN(x,2) should give A'*x)
%
%function [x, istop, itn, normr, normAr, normA, condA, normx]...
   %= lsmr(A, b, lambda, atol, btol, conlim, itnlim, localSize, show)

    b = vu(logical(mask));
    atol = 1e-6;
    btol = 1e-6;
    conlim = []; 
    itnlim = 20;
    localSize = 0;
    show = 1;
    [q,ISTOP,ITN,NORMR,NORMAR,NORMA,CONDA,NORMX] = lsmr(AFUN,b,atol,btol,conlim,itnlim,localSize,show);
    %[ q,  istop, itn ,r1norm,r2norm] =lsqrSOL(sum(sum(sum(mask))), prod(size(vu)), AFUN, vu(logical(mask)),L,1e-6,1e-6,0,20,1);
    toc
end
q = reshape(q,size(vu));
prev_q = reshape(prev_q,size(vu));
%v_est = xf2kt(abs(prev_q).*q); % v_est = Aq, no mask
v_est = fftc(abs(prev_q).*q,1);
err = L*norm(q(:)) + norm(reshape(vu-v_est.*mask,[],1),'fro')
disp('kt-FOCUSS completed');
