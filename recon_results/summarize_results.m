format compact;

result_files = dir('~/focuss/recon_results/');
filenames = { result_files.name };

addpath(genpath('..'));
load('2D_data.mat')
TRUTH = func_data; % func_data = data + mask
TNORM = norm(TRUTH(:));

for idx = 1:length(filenames)
  clearvars -except TRUTH TNORM idx filenames

  filename = filenames{idx};
  [path, name, ext] = fileparts(filename);
  if (strcmp(ext, '.mat'))
    vars = whos('-file', filename);
    if ~ismember('recon', {vars.name}) % todo, and not recon_z24.mat
        continue
    end
    display(name)
    load(filename,'recon');

    display(filename);
    err = norm(recon(:)-TRUTH(:))/TNORM

    x = 47; y = 20;
    figure(idx);
    hold on;
    plot(abs(squeeze(recon(y,x,:))));
    plot((squeeze(TRUTH(y,x,:))),'r');
    legend('recon','truth');
    title(strcat('flashing fmrib 2D_data recon: ', filename), 'interpreter', 'none');
    hold off;
    drawnow;
    plot_filename = strcat(fullfile(path, 'result_plots', name),'.pdf');
    saveas(idx, plot_filename, 'pdf');
  end
end

save('results_summary.mat', 'TRUTH');
