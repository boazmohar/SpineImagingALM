function fx = fn_fit(x,y,m,startpoint)
% function fx = fn_fit(x,y,model,startpoint)

% Thomas Deneux
% Copyright 2008-2012

if nargin==0, help fn_fit, return, end

if ~any(size(x)==1), error('x must be a vector'), end
x =  x(:);
if size(y,1)==1, y=y'; end

if ischar(m)
    m = fittype(m);
end

opt = fitoptions('method','NonlinearLeastSquares','startpoint',startpoint);

nk = size(y,2);
fx = cell(1,nk);
for k=1:nk
    yk = y(:,k);
    fx{k} = fit(x,yk,m,opt);
end

if nk==1, fx = fx{1}; end