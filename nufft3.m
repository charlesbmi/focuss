  function adj = nufft3(img, st)
% function adj = nufft3(img, st)
% Transforms img from x-y-t space to nufft2-t space

% todo pre-allocate and temporal transform
[nx,ny,nt] = size(img);
for t = 1:nt
    adj(:,t) = nufft(img(:,:,t),st(t));
end
