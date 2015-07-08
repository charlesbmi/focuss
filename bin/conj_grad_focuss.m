%function rho_n = focuss_est(dx, sample_loc, Wn, L, F, FT)
% Focuss parameters
clear all
load('2D_data.mat')
load('sampling_masks.mat')
tic

%x = zeros(size(kt));
mask = cart_sampling_mask_8x_4low_freq;
% todo maybe just fft in 1-D is fine.
kt = mask.*fft(fft(func_data,[],1),[],2);
F = @(x) ifft(fft(fft(x,[],1),[],2),[],3);
FT = @(x) ifft(ifft(fft(x,[],3),[],1),[],2);
rho = FT(kt);
%W = sqrt(abs(rho));
W = sqrt(abs(rho));
x = W;
x = zeros(size(W));

im = @(x) imshow(mat2gray(abs(x(:,:,1))));
L=0.1;
% descent parameters
ALPHA = 0.05;
BETA = 0.1;
MAXITERS = 50;
NTTOL = 1e-8;
GRADTOL = 1e-4;

% gradient method
vals = []; steps = [];
errnorm = norm(rho(:))
for iter = 1:MAXITERS
    % todo check that focuss cost is working
	[val, grad] = focuss_cost(kt, mask, W, L, x, F, FT);
	vals = [vals, val];
	dx = -grad;
	fprime = grad(:)'*dx(:);

    %%% Value printing
	iter
	val
	gradient_norm = norm(grad(:))
    %%% End Value printing

	if norm(grad(:)) < GRADTOL * errnorm, break; end;
	t = 1;
    % todo is focuss cost and alpha t and fprime all real?
	while ( focuss_cost(kt, mask, W, L, x + t*dx, F, FT) > ...
		val + ALPHA*t*fprime )
        % linesearch
		t = BETA*t;
	end;
	x = x+t*dx;
	step_size = norm(t*dx(:))
	steps = [steps,t];
end;
toc

rho_n = x.*W;
recon = ifft(rho_n,[],3);

x = 47; y = 20;
figure;
hold on;
plot(abs(squeeze(recon(y,x,:))));
plot((squeeze(func_data(y,x,:))),'r');
legend('recon','mask or func_data');
hold off;

