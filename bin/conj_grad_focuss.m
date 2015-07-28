%function rho_n = focuss_est(dx, sample_loc, Wn, L, F, FT)
% Focuss parameters
%clear all
load('2D_data.mat')
load('cart_sampling_masks.mat')
tic

mask = cart_sampling_mask_8x_4low_freq;
kt = mask.*fft(fft(func_data,[],1),[],2);
F = @(rho) ifft(fft(fft(rho,[],1),[],2),[],3);
FT = @(kt) ifft(ifft(fft(kt,[],3),[],1),[],2);

[U,S,V] = svd(reshape(recon, [64*64 512]));
F = @(rho) iklt3(fft(fft(rho,[],1),[],2),V);
FT = @(kt) ifft(ifft(klt3(kt,V),[],1),[],2);

rho = FT(kt);
W = sqrt(abs(rho));
x = zeros(size(W));

im = @(x) imshow(mat2gray(abs(x(:,:,1))));
L=0.1;
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

    %%% Value printing and use
	iter
	val
	prev_grad_norm = norm(grad(:))
    %%% End Value printing

	if norm(grad(:)) < GRADTOL * errnorm, break; end;
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
end;
toc

rho_n = x.*W;
recon = ifft(rho_n,[],3);
recon = iklt3(rho_n, V);

x = 47; y = 20;
figure;
hold on;
plot(abs(squeeze(recon(y,x,:))));
plot((squeeze(func_data(y,x,:))),'r');
legend('recon','mask or func_data');
hold off;

