function rho_n = focuss_est(kt, mask, Wn, L, F, FT)
% Computes the nth FOCUSS estimate for k-t FOCUSS algorithm by minimizing
% the cost function using Nielsen's conj_grad or MATLAB function lsqnonlin, as
% described in "Improved k-t BLAST and k-t SENSE using FOCUSS (Jung et al., 2007)
% Inputs
%       kt - raw undersampled k-t space data of object
%       mask - sampling locations of kt, 1 where sampled, 0 where
%           missing data
%       Wn - updated weighting vector
%       L - Lagrangian parameter, which weights cost function
%       F - function handle for sparsifying transform
%       FT - function handle for inverse of sparsifying transform
% Outputs
%       rho_n - nth iterative estimate of pruned sparse result
%
% Charles Guan (charles.guan@stanford.edu)

% Initialize variables
rho = FT(kt);
W = sqrt(abs(rho));
x = zeros(size(W));

kt = mask.*fft(fft(func_data,[],1),[],2);
F = @(rho) ifft(fft(fft(rho,[],1),[],2),[],3);
FT = @(kt) ifft(ifft(fft(kt,[],3),[],1),[],2);

%[U,S,V] = svd(reshape(recon, [64*64 512]));
%F = @(rho) iklt3(fft(fft(rho,[],1),[],2),V);
%FT = @(kt) ifft(ifft(klt3(kt,V),[],1),[],2);

rho = FT(kt);
W = sqrt(abs(rho));
x = zeros(size(W));

% descent parameters
ALPHA = 0.05;
BETA = 0.1;
MAXITERS = 100;
NTTOL = 1e-8;
GRADTOL = 1e-4;

% gradient method
vals = []; steps = [];
errnorm = norm(rho(:))
t = 1;
dx = 0;
prev_grad_norm = Inf;
for iter = 1:MAXITERS
	[val, grad] = focuss_cost(kt, mask, W, L, x, F, FT);
	vals = [vals, val];
	dx = -grad + dx*(norm(grad(:))/prev_grad_norm)^2;
	fprime = grad(:)'*dx(:);

        disp(sprintf('Iter: %03i, Grad: %f, Cost Funct: %f',iter,norm(g(:))/errnorm,value/errnorm));

	%t = 1; % speed up by "remembering" previous t. heuristic
        t = min(1, t/BETA^2);
        % todo is focuss cost and alpha t and fprime all real?
	while ( focuss_cost(kt, mask, W, L, x + t*dx, F, FT) > ...
		val + ALPHA*t*fprime )
        % linesearch
		t = BETA*t;
	end;
	x = x+t*dx;
	step_size = norm(t*dx(:))
	steps = [steps,t];

	if norm(t*dx(:)) < GRADTOL * norm(x(:)), break; end;
end;

rho_n = x.*W;
