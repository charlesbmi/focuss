  function img = nufft3_adj(adj, st)
% function img = nufft3_adj(adj, st)
% transforms non-uniform k-space data indexed at st values to cartesian image

% todo pre-allocate
nt = numel(st);
for t = 1:nt
    img(:,:,t) = nufft_adj(adj(:,t),st(t));
end
