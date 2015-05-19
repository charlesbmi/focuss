  function adj = mcnufft3(img, st, nl, coils)
% function adj = mcnufft3(img, st,nl, coils)

%	Modified by M. Chiew
%	Based on nufft3.m 
%
% Transforms img from x-y-t space to coil-nufft2-t space
% Returns in coils, num points, num lines, num time space

[nx,ny,nt] = size(img);
M = st(1).M;
for c = 1:size(coils,3)
	for t = nt:-1:1
		adj_frame = nufft(coils(:,:,c).*img(:,:,t),st(t));
	    adj(c,:,:,t) = reshape(adj_frame,[M/nl nl]);
	end
end
