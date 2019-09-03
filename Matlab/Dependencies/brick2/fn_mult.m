function y=fn_mult(u,v)
% function y=fn_mult(u,v)
%----
% tool to multiply a matrix row- or column-wise
% ex: y = fn_mult(rand(3,4),(1:3)')
%     y = fn_mult(rand(5,2,5),shiftdim(ones(5,1),-2))
%
% See also fn_add

% Thomas Deneux
% Copyright 2002-2012

s1 = size(u); s2 = size(v);
if length(s2)>length(s1), y = fn_mult(v,u); return, end
if length(s1)>length(s2), for i=length(s2)+1:length(s1), s2(i)=1; end, end
if ~all(s1==1 | s2==1 | s1==s2), error('dimension mismatch'), end
s = max(s1,s2);

u = repmat(u,s-s1+1);
v = repmat(v,s-s2+1);
y = u.*v;
        
            
            
    
    