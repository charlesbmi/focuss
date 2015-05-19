  function img = nufft3_adj(adj, st)
% function img = nufft3_adj(adj, st)
% transforms non-uniform k-space data indexed at st values to cartesian image

nt = numel(st);
for t = nt:-1:1
	adj_frame = adj(:,:,t);
    img(:,:,t) = nufft_adj(adj_frame(:),st(t));
end
