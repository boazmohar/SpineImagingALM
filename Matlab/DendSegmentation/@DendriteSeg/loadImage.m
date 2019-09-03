function loadImage(obj)
% loading an image and preprossesing it
set(obj.handles.ustatH,'String','Loading image')
[Name,Path] = uigetfile('*.tif','Select an image');
image=readtiff(Path,[],Name);
%% getting user input for pixels size and radiuos
prompt = {'X size:','Y size:','Z size:', 'Radius'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0.4','0.4','2','0.5'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
pixelSize = [0,0,0];
pixelSize(1) = str2double(answer{1});
pixelSize(2) = str2double(answer{2});
pixelSize(3) = str2double(answer{3});
radius = str2double(answer{4});
%%
set(obj.handles.ustatH,'String','Upsampleing image')
image(isnan(image))=0;
si=size(image);
A=eye(4)*2;
A(4,4)=1;
R=makeresampler('cubic','bound');
tform=maketform('affine',A);
imageUp = tformarray(image, tform, R, [1 2 3], [1 2 3], [si(1)*2 si(2)*2, si(3)*2],[],0);
%% top hat transform
set(obj.handles.ustatH,'String','Enhancing image')
pixelSize=pixelSize./2;
imageUp=imtophat(imageUp,distStrel3D(radius, pixelSize));
obj.display.volume = imageUp;
uisave('imageUp')
set(obj.handles.ustatH,'String','Ready')
end

