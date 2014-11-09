function im_mr = multicoil_recon(im, map)
% Performs a multicoil reconstruction using coil sensitivity data
% This assumes that the noise correlation matrix is identity I
% Inputs 
%   im - multicoil image data in (x,y,t/f,coil)-domain
%   map - complex sensitivities for each channel,
%         map is the same size as im if time-varying (x,y,t,coil)
%         map is 1 dim less if not (x,y,coil)
% Outputs
%   im_mr - full multicoil reconstruction of im, 1 dim less than im
% Note: also works for single-frame images and maps in (x,y,coil) domain
%
% Charles Guan (charles.guan@stanford.edu)

if ndims(im) == 4 && ndims(map) == 3
    % Replicate time-constant sensitivity map across time frames
    [nx,ny,nt,nc] = size(im);
    map = reshape(map,[nx ny 1 nc]);
    map = repmat(map,[1 1 nt 1]);
end
im_mr = sum(im.*conj(map),ndims(im))./sum(map.*conj(map),ndims(im));
