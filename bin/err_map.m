function map = err_map(ts_est, ts)
% Creates an image of the norm error
% Inputs
%   ts_est - estimated time-series
%   ts - original time-series
% Returns
%   map - image of 2-norm error, same size as one frame of ts

err = ts-ts_est;
map = sqrt(sum(abs(err.^2),3));
