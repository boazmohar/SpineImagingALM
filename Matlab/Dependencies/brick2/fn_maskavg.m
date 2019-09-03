function a = fn_maskavg(a,mask)
% function a = fn_maskavg(a,mask)
%---
% region-wise averaging of a according to a 2-dimensional mask

% Thomas Deneux
% Copyright 2007-2012

labels = unique(mask(:));
nreg = length(labels);

s = size(a);
[nx ny nt] = size(a);

a = reshape(a,[nx*ny nt]);

for ireg=1:nreg
    ind = find(mask==labels(ireg));
    % average in region
    ak = mean(a(ind,:),1);
    a(ind,:) = repmat(ak,length(ind),1);
end
a = reshape(a,s);
