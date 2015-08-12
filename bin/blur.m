image1, image2, roi

filename = '';
%load(filename);

% assert
% nx,ny,nt = size(
recon_auto_corr_map = {};
true_vs_recon_corr_map = {};
nr = nl = length(roi);

for idx = 1:length(roi)
    letter = roi{idx};
    recon_auto_corr_map{idx} = zeros(nx, ny);
    true_vs_recon_corr_map{idx} = zeros(nx, ny);

    for x = 1:nx
        for y = 1:ny
            letter_x = mod(roi{4}, nx); % todo may need to switch x and y. also todo not sure if all recon letter are correct given the 0 - I. Check manually.
            letter_y = floor(roi{4} / ny);
            recon_auto_corr_map{idx} = corr(func_data(letter_x,letter_y,:), recon(x,y,:));
            true_vs_recon_corr_map{idx} = corr(recon(letter_x,letter_y,:), recon(x,y,:));
        end
    end

end

% save(filename
% save back to same location with new variables

% Tweak this and send it over to Mark with some results tonight.
