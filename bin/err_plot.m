function ets = err_plot(ts_est, ts)
% Creates a time series for the frobenius error norm
% Inputs
%   ts_est - estimated time-series
%   ts - original time-series
% Returns
%   ets - time-series of frobenius norm error at each time slice

[nx,ny,nt] = size(ts);
err = reshape(ts_est-ts,nx*ny,nt);
ets = sqrt(sum(abs(err.^2),1));
