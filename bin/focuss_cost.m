function [ cost, grad ] = focuss_cost(v, mask, W, L, q, A, AT)
% Computes the cost of a solution q for rho = W*q, as described in 
% "Improved k-t BLAST and k-t SENSE using FOCUSS" by Jung, 2007. Page 3206
% Parameters and variables have been named as in the paper
% Note: W is a real scaling vector, instead of a diagonal matrix
% mask is required in norm calculation because cost function is
% only computed for acquired data points
% Inputs
%       v - acquired, undersampled k-t space data
%       mask - sampling locations of v, 1 where sampled, 0 where
%           missing data
%       vsize - original size of v, mask, W, and q,
%           needed because CG algorithms frequently only take vectors as I/O
%       W - weighting vector
%       L - Lagrangian parameter, which weights cost function
%       q - iterative estimate. Cost is a function of q
%       A - function handle for sparsifying transform
%       AT - function handle for inverse of sparsifying transform
% Outputs
%       cost - cost of q for optimization problem (Equation 12, Jung)
%       grad - dC/dq, the gradient (partial derivative) of C with respect to q (Equation 38, Jung)
%
% Charles Guan (charles.guan@stanford.edu)

% Reshaping necessary because conj_grad takes 1-D vectors as input and output
%q = reshape(q,vsize);
err_term = (v-A(W.*q)).*mask;
cost = norm(err_term(:),2)^2 + L*norm(q(:),2)^2;
if nargout > 1
	grad = -W.*AT(err_term) + L*q;
end
%grad = grad(:);
