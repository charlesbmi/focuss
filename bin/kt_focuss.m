function X_FOCUSS = kt_focuss(A,AT,Y,mask,num_low_freq);
% Wrapper for more KTFOCUSS defining default options
% Solves A*X = Y for X
% num_low_freq - number of low frequencies lines fully sampled, to use for low-freq estimate
% % parameter setting
% Y is old DownSino = Mask.*kt_data
Low_Y = Y;
Low_Y(num_low_freq/2+1:end-num_low_freq/2,:,:) = 0; % Low_Y definition only works for cartesian case
focuss_iters = 6; % paper states that more than 5 iterations provides minimal improvements
inner_loop_iters = 100;
factor = 0.5;
lambda = 0.1;
X_FOCUSS = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda, inner_loop_iters, focuss_iters);
