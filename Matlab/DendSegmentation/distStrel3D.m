function SE=distStrel3D(radius, pixelSize)
w=ceil(radius./pixelSize.*2);
[X,Y,Z]=meshgrid(-w(1):w(1),-w(2):w(2),-w(3):w(3));
d=(((X.*pixelSize(1)).^2+(Y.*pixelSize(2)).^2+(Z.*pixelSize(3)).^2).^(1/2));
C=exp(-(d.^2/(2.*radius)));
C=C-min(C(:));
C=C./max(C(:));
SE=strel('arbitrary',ones(size(C)),C);