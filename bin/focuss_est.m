function rho_n = focuss_est(v, sample_loc, Wn, L, F, FT)
% Computes the nth FOCUSS estimate for k-t FOCUSS algorithm by minimizing
% the cost function using Nielsen's conj_grad or MATLAB function lsqnonlin, as
% described in "Improved k-t BLAST and k-t SENSE using FOCUSS (Jung et al., 2007)
% Inputs
%       v - raw undersampled k-t space data of object
%       sample_loc - sampling locations of v, 1 where sampled, 0 where
%           missing data
%       Wn - updated weighting vector
%       L - Lagrangian parameter, which weights cost function
%       F - function handle for sparsifying transform
%       FT - function handle for inverse of sparsifying transform
% Outputs
%       rho_n - nth iterative estimate of pruned sparse result
%
% Charles Guan (charles.guan@stanford.edu)

% Using Nielsen's conj_grad
% cost_fun = @(q,~)focuss_cost(v,size(v),sample_loc,Wn,L,q,F,FT);
% [qn, infor] = conj_grad(cost_fun,[],Wn(:),[2,2,1,1e-10,1e-10,1000,1e-10,1e-10,100]);
% switch infor(6)
%     case 1
%         disp(sprintf('conj_grad converged at evaluation %d by small gradient ||g||_inf = %d',infor(5),infor(2)));
%     case 2
%         disp(sprintf('conj_grad converged at evaluation %d by small x-step ||dx||_2 = %d',infor(5),infor(3)));
%     case 3
%         disp(sprintf('conj_grad failed to converge: maximum evaluation of %d reached. gradient residual: %d',infor(5),infor(2))); 
% end
% disp(sprintf('Final function evaluation: %g',infor(1)));
% disp(sprintf('Final gradient evaluation: %g',infor(2)));
% disp(sprintf('Final x-step: %g',infor(3)));
% qn = reshape(qn, size(v));

% Using MATLAB's optimization toolbox lsqnonlin
options = optimset('lsqnonlin');
options = optimset(options,'Jacobian','on');
options = optimset(options,'DerivativeCheck','on');
option = optimset(options,'FunValCheck','on');
options = optimset(options,'TolFun',1e-3,'TolX',1e-3);
lsq_cost_fun = @(q)focuss_lsq_cost(v,sample_loc,Wn,L,q,F,FT);
qn = lsqnonlin(lsq_cost_fun,Wn,[],[],options);

rho_n = Wn.*qn;
