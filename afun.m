function est = afun(x, prev_q, img_size, mask, L, topt)
% Function handle with extracting 
% x - n x 1 flattened x-f image rho = Wq if topt == 'notransp'
% if topt (transpose option) == 'transp, x = flattened image v
% prev_q - previous estimate of q
% img_size - 3-D size of x
% trans - Sparsifying transform, ie A is the function trans along with reshaping
% mask - 1s and 0s sampling matrix
% L - square of lambda weighting of objective functions (ie sparsity penalty)
% tranpose option

% A = | Sampling * Fourier * W |
%     | lambda                 |
% todo comment this
A = @(q) mask.*xf2kt(abs(prev_q).*q);
AT = @(v) kt2xf(v.*mask).*abs(prev_q);
ATA = @(q) AT(A(q));

if strcmp(topt,'transp') % compute A'*x
    nsamples = sum(sum(sum(mask)));
    v(mask) = x(1:nsamples);
    num_trailing_zeros = prod(img_size)-numel(v);
    trailing_zeros = zeros([1,num_trailing_zeros]);
    v = reshape([v,trailing_zeros], img_size);
    q_est = AT(v);
    est = q_est(:) + L*x(nsamples+1:end);
else % compute A*x
    q = reshape(x,img_size);
    v_est = A(q);
    v_est = v_est(logical(mask));
    est = [v_est; L*x];
end
