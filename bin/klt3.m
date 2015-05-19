function [ xyklt ] = klt3(xyt, V)
% Transforms 3D x-y-t data (via x-y-klt) to  x-y-klt data
% Inputs
%		t - (x,y,t) data to transform
%		V - right-singular vectors of 2-D reshaped xyt estimate, as given by [~,~,V] = svd(kt_estimate_reshaped)
% Output
%		xyklt - data in (x,y,klt) basis, with KLT basis from input
%
% Charles Guan (charles.guan@stanford.edu)

xytr = reshape(xyt,[size(xyt,1)*size(xyt,2),size(xyt,3)]);
xykltr = xytr*V;
xyklt = reshape(xykltr, size(xyt));
