result_files = dir('~/focuss/recon_results/');
filenames = { result_files.name};

load('2D_data.mat')
TRUTH_DATA = func_data; % func_data = data + mask
TD_NORM = norm(TRUTH_DATA(:));

for idx = 1:length(filenames)
    filename = filenames{idx};
    [path, name, ext] = fileparts(filename);
    if (strcmp(ext, '.mat'))
      display(name)
      fileobj = matfile(filename);
      filedata = whos(fileobj);
      filedata.name;
      load(filename,'recon');
    end
end
