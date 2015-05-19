function est = lfun(x, prev_q, img_size, mask, topt)
% Function handle for LSMR 
% x - n x 1 flattened x-f image rho = Wq if topt == 'notransp'
% if topt (transpose option) == 'transp, x = flattened image v
% prev_q - previous estimate of q
% img_size - 3-D size of x
% trans - Sparsifying transform, ie A is the function trans along with reshaping
% mask - logical sampling matrix
% L - square of lambda weighting of objective functions (ie sparsity penalty)
% tranpose option

A = @(q) mask.*xf2kt(abs(prev_q).*q); % Sampling * Fourier * W
AT = @(v) kt2xf(v.*mask).*abs(prev_q);

if topt == 2
    v(mask) = x;
    num_trailing_zeros = prod(img_size)-numel(v);
    trailing_zeros = zeros([1,num_trailing_zeros]);
    v = reshape([v,trailing_zeros], img_size);
    q_est = AT(v);
    est = q_est(:);
else
    q = reshape(x,img_size);
    v_est = A(q);
    est = v_est(logical(mask));
end
