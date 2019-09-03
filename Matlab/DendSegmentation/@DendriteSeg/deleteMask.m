function deleteMask(obj,X,Y,Z)
obj.saving.changed = 1;
ind = obj.cells.labelSpine(Y,X,Z);
obj.cells.labelSpine(obj.cells.labelSpine == ind) = 0;
obj.cells.maskSpine = obj.cells.labelSpine>0;
if obj.cells.labelBorder(Y,X,Z) > 0
    obj.cells.labelBorder(obj.cells.labelBorder == ind) = 0;
    obj.cells.maskBorder = obj.cells.labelBorder>0;
    fprintf(1, 'Deleted border %d: %d pix\n', ind, obj.cells.all{ind}.pxNum);
else
    fprintf(1, 'Deleted obj %d: %d pix\n', ind, obj.cells.all{ind}.pxNum);
end
obj.updateImage();