function [X, Y, Z] = c_PadZDense(obj,X,Y,Z)
if length(X) > 3
    dx1 = X(2) - X(1);
    dx2 = X(end) - X(end-1);
    dy1 = Y(2) - Y(1);
    dy2 = Y(end) - Y(end-1);
    dz1 = Z(2) - Z(1);
    dz2 = Z(end) - Z(end-1);
    X = [X(1)-dx1*2;    X(1)-dx1;  X; X(end)+dx2;   X(end)+dx2*2;];
    Y = [Y(1)-dy1*2;    Y(1)-dy1;  Y; Y(end)+dy2;   Y(end)+dy2*2;];
    Z = [Z(1)-dz1*2;    Z(1)-dz1;  Z; Z(end)+dz2;   Z(end)+dz2*2;];
else
    X = [X(1); X(1); X; X(end); X(end)];
    Y = [Y(1); Y(1); Y; Y(end); Y(end)];
    Z = [Z(1)-2; Z(1)-2; Z; Z(end)+2; Z(end)+4];
end


