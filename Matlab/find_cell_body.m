pixelSize=[0.4 0.4 2];
PSF = [.4 .4 2.5];
A=eye(4)*0.5;
A(4,4)=1;
R=makeresampler('cubic','bound');
tform=maketform('affine',A);
%%
structure = distStrel3D(5, pixelSize);
%%

cell_vol=readtiff('V:\users\Aaron\Database\BMWR70\FOV2',[],'Stack_32bit');
%%
cell_vol2 = cell_vol(301:800, 301:800, 51:100);
%%
cell_vol2(isnan(cell_vol2))=0;
si=size(cell_vol2);
disp('upsample')
tic
cell_vol_up2 = tformarray(cell_vol2, tform, R, [1 2 3], [1 2 3], ...
                         [si(1)*0.5 si(2)*0.5, si(3)*0.5],[],0);
toc
%%
writetiff(cell_vol_up2,'cell_vol_down')
%%
[center_img,sphere_img,cent,radi]=SphericalHough(cell_vol_up2,[4 15], 0.001, 4, 0.1, 0);

%%
for ind=1:length(radi)
fprintf('Obj %d - ',ind)
fprintf('cen: %d,%d,%d; ',round(cent(ind,:)))
fprintf('rad: %1g pixels\n',radi(ind))
end