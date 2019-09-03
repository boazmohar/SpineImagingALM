expdir='I:\BMWR35\FOV2';
stackName='Stack2Viewtarget2';
stack=readtiff(expdir,[],[stackName,'.tif']);
fullRangeScale=1.035;

siOrig=size(stack);
rescaledStack=stack;
scaleStep=(fullRangeScale-1)./siOrig(3);
for i=1:siOrig(3)
    resized=imresize(stack(:,:,i),(1+scaleStep*i));
    si=size(resized);
    rescaledStack(:,:,i)=imcrop(resized,[round((si(1)-siOrig(1))./2),round((si(2)-siOrig(2))./2),siOrig(1)-1,siOrig(2)-1]);
end

writetiff(rescaledStack,[expdir,'\',stackName,'Rescaled.tif']);
