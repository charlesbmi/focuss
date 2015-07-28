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
  disp(path);
  if (strcmp(ext, '.mat'))
    display(name)
    fileobj = matfile(filename);
    filedata = whos(fileobj);
    filedata.name;
    load(filename,'recon');

    display(filename);
    err = norm(recon(:)-TRUTH(:))/TNORM

    x = 47; y = 20;
    figure(idx);
    hold on;
    plot(abs(squeeze(recon(y,x,:))));
    plot((squeeze(full_sample_img(y,x,:))),'r');
    legend('recon','truth');
    title(strcat('flashing fmrib 2D_data recon: ', filename));
    hold off;
    plot_filename = strcat(fullfile(path, name),'.pdf');
    saveas(idx, plot_filename, 'pdf');
  end
end

save('results_summary.mat', 'TRUTH');
