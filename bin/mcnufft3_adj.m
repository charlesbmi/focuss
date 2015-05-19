  function img = mcnufft3_adj(adj, st, coils)
% function img = mcnufft3_adj(adj, st, coils)

%	Modified by M. Chiew
%	Based on nufft3_adj.m 

% transforms non-uniform multi-coil k-space data indexed at st values to cartesian image

nt = numel(st);
for c = 1:size(coils,3)
	for t = nt:-1:1
		adj_frame = adj(c,:,:,t);
	    img(:,:,t,c) = nufft_adj(adj_frame(:),st(t));
	end
end
img = multicoil_recon(img, coils);