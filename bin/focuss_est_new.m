function [ rho_n, vals, steps, step_norms] = focuss_est_new(kt, mask, W, L, F, FT)
% Computes the nth FOCUSS estimate for k-t FOCUSS algorithm by minimizing
% Cost(q) = ||kt − F*W*q||^2_2 + L*||q||^2_2, as
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

tic

% descent parameters
ALPHA = 0.01;
BETA = 0.5;
MAXITERS = 50;
NTTOL = 1e-8;
GRADTOL = 1e-3;
%q = W; % last estimate chosen as initial guess. For some reason this leads to stagnant, slow-slope convergence.
q = zeros(size(W));
vals = []; steps = []; step_norms = [];t=1;

% conjugate gradient method
% Hessian is F*W + L*I
for iter = 1:MAXITERS
	iter
	[val, grad] = focuss_cost(kt, mask, W, L, q, F, FT);
	vals = [vals, val];
	v = -grad; % Using conjugate gradient
	%v = -grad./hessian;% Using Newton's method
	fprime = grad(:)'*v(:);
	val
	if norm(grad(:)) < GRADTOL, break; end;
	t = min(1,t/BETA^2); % Start with lower t to help faster convergence
	% t = 1;
	while ( focuss_cost(kt, mask, W, L, q + t*v, F, FT) > ...
		val + ALPHA*t*fprime )
		t = BETA*t;
	end;
	q = q+t*v;
	steps = [steps,t];
	step_norms = [step_norms, norm(t*v(:))];
end;
toc

rho_n = W.*q;