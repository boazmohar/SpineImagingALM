function nhood=distNhood3D(radius, pixelSize)
w=round((radius./pixelSize));
[X,Y,Z]=meshgrid(-w(1):w(1),-w(2):w(2),-w(3):w(3));
d=(((X.*pixelSize(1)).^2+(Y.*pixelSize(2)).^2+(Z.*pixelSize(3)).^2).^(1/2));
nhood=d<radius;