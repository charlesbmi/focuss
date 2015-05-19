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
%       mask - logical sampling locations of vu, 1 where sampled, 0 where
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

% TODO maybe key is to just do FFT in undersampling domain just like the original implementation does
clear;
disp('Initializing data');
addpath(genpath('bin'));
addpath(genpath('data'));

%% Load full measurement 
filename = ['old_test_data.mat']; % load full x-y-t data, and coils
disp(['Loading data from: ',filename]);
load(filename);
disp('Loaded');

v = fftc(fftc(data_2D,2),1);
mask = repmat(rand([size(v,1),1,size(v,3)]) < 0.5,[1,size(v,2),1]);
vu = v.*mask;
kt_avg = mean(vu,3); % DC component of X-F support
kt_avg = repmat(kt_avg, [1,1,size(v,3)]);
%vd = vu-kt_avg; % subtract DC component for X-F

disp('Initializing k-t FOCUSS')
L = 1e-1;   % lambda weighting multiobjective
updates = 1;
p = 0.5;    % from paper
% Optimize FFT
fftw('planner','patient');

% Resolve function handles
disp(sprintf('sparsifying transform selected: %s','Fourier'))
%ft = @(xf) ifft(fft(xf,[],1),[],2); % uncentered
%ift = @(kt) fft(ifft(xf,[],1),[],2); % centered

q = kt2xf(vu); % initial estimate
q = q./sqrt(abs(q));
for iter = 1:updates
    prev_q = reshape(q,size(vu)); % reshape if not iteration 1
    disp(sprintf('kt-FOCUSS iteration #: %d',iter))
    tic
    % lsqr
    AFUN = @(x,topt) afun(x,prev_q,size(vu),mask,L,topt);
    b = [vu(mask); zeros(numel(q),1)];
    tol = [];
    maxit = 20;
    M1 = [];
    M2 = [];
    x0 = prev_q(:);
    %q = lsqr(AFUN,b,tol,maxit,M1,M2,x0);

    % lsmr
    m = sum(sum(sum(mask)));
    n = numel(vu);
    AFUN = @(x,topt) lfun(x,prev_q,size(vu),mask,topt);
    b = vu(mask);
    atol = 1e-6;
    btol = 1e-6;
    conlim = []; % default
    itnlim = 100;
    localSize = []; % default
    show = true;
    x0 = prev_q(:);
    %[q,ISTOP,ITN,NORMR,NORMAR,NORMA,CONDA,NORMX] = lsmr(AFUN,b,L,atol,btol,conlim,itnlim,localSize,show,prev_q(:));
    [q,ISTOP,ITN] = lsmr(AFUN,b,L,atol,btol, conlim,itnlim,localSize,show,x0);
    %q = lsqrSOL(m, n, AFUN, b, L, atol, btol, conlim, itnlim, show);
    toc

v_est = lfun(q,prev_q,size(vu),ones(size(mask)),1); % v_est = Aq, no mask
%v_est = afun(q,prev_q,size(vu),ones(size(mask)),L,'notransp'); % v_est = Aq, no mask
% residual: AFUN(q)-vu
%res = L*norm(q(:)) + norm(vu(:)-v_est(1:numel(v)).*mask(:),'fro')
res = L*norm(q(:))^2 + norm(b-AFUN(q,1))^2;
res = norm(b-AFUN(q,1))
% error: AFUN(q)-v
err = norm(v(:)-v_est(:),'fro')
end
% add back in kt_avg dc component
disp('kt-FOCUSS completed');
