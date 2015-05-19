function [ xt ] = iklt3(xyklt, V)
% Transforms 3D x-y-klt data (via kx-ky-klt) into kx-ky-t data
% Inputs
%		xyklt - (x,y,klt) data to transform
%		V - right-singular vectors of 2-D reshaped kt estimate, as given by [~,~,V] = svd(kt_estimate_reshaped)
% Output
%		t - data in (x,y,t) basis, with KLT basis from input
%
% Charles Guan (charles.guan@stanford.edu)

xkltr = reshape(xyklt,[size(xyklt,1)*size(xyklt,2),size(xyklt,3)]);
xtr = xkltr*V';
xt = reshape(xtr,size(xyklt));
