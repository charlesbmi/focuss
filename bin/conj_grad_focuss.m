%function rho_n = focuss_est(v, sample_loc, Wn, L, F, FT)
% Focuss parameters
load('2D_data.mat')
tic
L=0.1;
%q = zeros(size(kt));
F = @(x) ifft(fft(x,[],1),[],3);
FT = @(x) ifft(fft(x,[],3),[],1);
rho = fft(func_data,[],3);
kt = fft(fft(func_data,[],1),[],2);
W = sqrt(abs(rho));
im = @(x) imshow(mat2gray(abs(x(:,:,1))));
q = zeros(size(kt));
q = W;

% descent parameters
ALPHA = 0.05;
BETA = 0.1;
MAXITERS = 30;
NTTOL = 1e-8;
GRADTOL = 1e-3;

% gradient method
vals = []; steps = [];
for iter = 1:MAXITERS
	iter
    % todo check 
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
	step_size = norm(t*v(:))
	steps = [steps,t];
end;
toc

rho_n = q.*W;
