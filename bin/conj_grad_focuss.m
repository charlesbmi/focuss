%function rho_n = focuss_est(v, sample_loc, Wn, L, F, FT)
% Focuss parameters
tic
L=0.1;
%q = zeros(size(kt));
F = @(x) ifft(fft(x,[],1),[],3);
FT = @(x) ifft(fft(x,[],3),[],1);
F = @(x) ifft(fft(fft(x,[],1),[],2),[],3);
FT = @(x) ifft(ifft(fft(x,[],3),[],1),[],2);
nlf = 6;
rho = fft(ifft(ifft(kt,[],1),[],2),[],3);
W = sqrt(abs(rho));
im = @(x) imshow(mat2gray(abs(x(:,:,1))));
q = W;

% descent parameters
ALPHA = 1e-4;
BETA = 0.1;
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
	step_size = norm(t*v(:))
	steps = [steps,t];
end;
toc