  function adj = nufft3(img, st, nl)
% function adj = nufft3(img, st)
% Transforms img from x-y-t space to nufft2-t space
% Returns in num points, num lines, num time space

[nx,ny,nt] = size(img);
M = st(1).M;
for t = nt:-1:1
	adj_frame = nufft(img(:,:,t),st(t));
    adj(:,:,t) = reshape(adj_frame,[M/nl nl]);
end
