function [ rho_n, costs, steps, step_norms] = focuss_est_new(kt, mask, W, L, F, FT)
% Computes the nth FOCUSS estimate for k-t FOCUSS algorithm by minimizing
% Cost(q) = ||kt - F*W*q||^2_2 + L*||q||^2_2, as
% described in "Improved k-t BLAST and k-t SENSE using FOCUSS (Jung et al., 2007)
% Inputs
%       kt - raw undersampled k-t space data of object
%       mask - sampling locations of v, 1 where sampled, 0 where
%            missing data
%       W - updated weighting vector
%       L - Lagrangian parameter, which weights cost function
%       F - function handle for sparsifying transform
%       FT - function handle for inverse of sparsifying transform
% Outputs
%       rho_n - nth iterative estimate of pruned sparse result
%
% Charles Guan (charles.guan@stanford.edu)

% descent parameters
ALPHA = 0.1; %0.05;
BETA = 0.5; %0.1;
MAXITERS = 50;
NTTOL = 1e-8;
GRADTOL = 1e-4;
%q = W; % last estimate chosen as initial guess. For some reason this leads to stagnant, slow-slope convergence.
q = zeros(size(W));
vals = []; steps = []; step_norms = [];t=1;
dq = 0; prev_grad_norm = Inf;
rho = FT(kt);
errnorm = norm(rho(:));

for iter = 1:MAXITERS
    [cost, grad] = focuss_cost(kt, mask, W, L, q, F, FT);
    costs = [costs, cost];
    dq = -grad + dq*(norm(grad(:))/prev_grad_norm)^2;
    fprime = grad(:)'*dq(:);

    if norm(grad(:)) < GRADTOL * errnorm, break; end;
    %t = 1; % speed up by "remembering" previous t. heuristic
    t = min(1, t/BETA^2);
    % todo is focuss cost and alpha t and fprime all real?
    while ( focuss_cost(kt, mask, W, L, q + t*dq, F, FT) > ...
            val + ALPHA*t*fprime )
    % linesearch
            t = BETA*t;
    end;
    q = q+t*dq;
    step_size = norm(t*dq(:));
    steps = [steps,t];
    disp(sprintf('Iter: %03i, Step: %f, Grad: %f, Cost: %f',iter,step_size,norm(grad(:))/errnorm,sqrt(focuss_cost(kt, mask, W, L, q, F, FT))/errnorm));
end;

rho_n = q.*W;
recon = ifft(rho_n,[],3);
%recon = iklt3(rho_n, V);
