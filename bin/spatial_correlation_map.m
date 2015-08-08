function [corr_map] = spatial_correlation_map(ts1, ts2)
%function [corr_map] = spatial_correlation_map(ts1, ts2)
%
% 
%   Charles Guan
%   August 2015
%
%   Inputs:     ts1           size: [nx, ny, nt]
%               ts2           size: [nx, ny, nt]
%   Outputs:    corr_map      size: [nx, ny]

assert(all(size(ts1) == size(ts2)));
[nx, ny, nt] = size(ts1);
