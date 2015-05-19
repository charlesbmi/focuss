function [ FUN, J ] = focuss_lsq_cost(v, sample_loc, W, L, q, A, AT)
% Computes the function FUN(q) whose sum of squares is minimized as
% cost of a solution q for rho = W*q, as described in 
% "Improved k-t BlAST and k-t SENSE using FOCUSS" by Jung, 2007. Page 3206
% Parameters and variables have been named as in the paper
% Note: W is a real scaling vector, instead of a diagonal matrix
% sample_loc is required in norm calculation because cost function is
% only computed for acquired data points
% Inputs
%       v - acquired, undersampled k-t space data
%       sample_loc - sampling locations of v, 1 where sampled, 0 where
%           missing data
%       W - weighting vector
%       L - Lagrangian parameter, which weights cost function
%       q - iterative estimate. Cost is a function of q
%       A - function handle for sparsifying transform
%       AT - function handle for transpose of A
% Outputs
%       FUN - cost of q for optimization problem (Equation 12, Jung)
%       J - jacobian of FUN
%
% Charles Guan (charles.guan@stanford.edu)

err_term = (v-A(W.*q)).*sample_loc;
FUN = [norm(err_term(:)); norm(sqrt(L)*q(:))];

if nargout > 1 % Two output arguments
	J1 = -W.*AT(err_term);
    J(1,:) = J1(:);
    J(2,:) = L*q(:);
end
