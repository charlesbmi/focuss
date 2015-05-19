function [ rho_n, vals, steps, step_norms] = focuss_est_new(kt, mask, Wn, L, F, FT)
% Computes the nth FOCUSS estimate for k-t FOCUSS algorithm by minimizing
% the cost function using Nielsen's conj_grad or MATLAB function lsqnonlin, as
% described in "Improved k-t BLAST and k-t SENSE using FOCUSS (Jung et al., 2007)
% Inputs
%       kt - raw undersampled k-t space data of object
%       mask - sampling locations of v, 1 where sampled, 0 where
%            missing data
%       Wn - updated weighting vector
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
MAXITERS = 100;
NTTOL = 1e-8;
GRADTOL = 1e-3;

% gradient method
vals = []; steps = [];
for iter = 1:MAXITERS
	iter
	[val, grad] = focuss_cost(kt, mask, W, L, q, F, FT);
	vals = [vals, val];
	v = -grad;
	fprime = grad(:)'*v(:);
	val
	gradient_norm = norm(grad(:))
	if norm(grad(:)) < GRADTOL, break; end;
	t = 1;
	while ( focuss_cost(kt, mask, W, L, q + t*v, F, FT) > ...
		val + ALPHA*t*fprime )
		t = BETA*t;
	end;
	q = q+t*v;
	steps = [steps,t];
	step_norms = [step_norms, norm(t*v(:))];
end;
toc