function X_FOCUSS = mc_kt_focuss(A,AT,Y,mask,num_low_freq,map);
% multicoil k-t FOCUSS
% map - coil sensitivity map
% sparse solution for A*X_FOCUSS = Y

% parallel map function instead of 

nc = size(map,4);
recon = zeros(size(map)); % pre-allocate

disp('Parallel multi-coil k-t FOCUSS loop');
parfor coil = 1:nc
    disp(sprintf('Applying k-t FOCUSS to coil #%d',coil));
    recon(:,:,:,coil) = kt_focuss(A,AT,Y(:,:,:,coil),mask,num_low_freq);
    disp(sprintf('Finished k-t FOCUSS on coil #%d',coil));
end
X_FOCUSS = multicoil_recon(recon,map);
