clear all;
load('data/2D_data.mat')
addpath('bin')
%% Cartesian downsampling rates
rng(0)
[nx ny nt] = size(func_data);
for ds_rate = [2, 4, 8, 12, 16]
    for ds_pat = [0, 1]; % downsampling pattern, 0 for uniform random in ky, 1 for gaussian in ky, 2 for random in kx and ky % option to change
        for num_low_freq = [2, 4]; % option to change

            if (ds_pat == 1) & (num_low_freq == 4) & (ds_rate == 16)
              continue
            end

            switch ds_pat
            case 0
                ds_pat_str = 'cart';
            case 1
                ds_pat_str = 'gauss';
            otherwise
                ds_pat_str = 'uniform random in kx and ky';
            end

            disp(['Downsample rate: ', num2str(ds_rate)]);
            disp(['Downsample pattern: ', ds_pat_str, ' random in ky']);

            % low frequency full sampling number (1~num_phase/2), (end-num_phase/2~end) -> full sampling
            disp(['Low frequency full sampling number: ', num2str(num_low_freq)]);

            mask = downsample_mask(nx,ny,nt,ds_rate,num_low_freq,ds_pat);

            var_name = sprintf('%s_sampling_mask_%dx_%dlow_freq', ds_pat_str, ds_rate, num_low_freq)
            S.(var_name) = mask;
        end
    end
end

filename = 'data/cart_sampling_masks.mat';
display(['masks to file: ', filename])
save(filename, '-struct', 'S');

clearvars S
% radial sampling patterns
N_readout_points = 128;
for nl = [2,4,8,16,32,64]
    disp(['Radial sampling with lines / time frame: ', num2str(nl)]);
    [k, wi] = gen_radial(0, N_readout_points, nl*nt, 1);
    var_name = sprintf('radial_sampling_mask_%dreadout_lines', nl)
    S.(sprintf('%s_k', var_name)) = k;
    S.(sprintf('%s_wi', var_name)) = wi;
end

filename = 'data/radial_sampling_masks.mat';
display(['masks to file: ', filename])
save(filename, '-struct', 'S');
