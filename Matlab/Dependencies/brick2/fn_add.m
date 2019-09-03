function y = fn_add(u,v)
% function y = fn_add(u,v)
%---
% tool to add matrices and vectors
% ex: y = fn_add(rand(4,5),(1:4)')
%     y = fn_add(1:5,(1:4)')
%
% See also fn_mult

% Thomas Deneux
% Copyright 2002-2012


s1 = size(u); 
s2 = size(v);
n = max(length(s1),length(s2));
s1(end+1:n) = 1;
s2(end+1:n) = 1;
if ~all(s1==1 | s2==1 | s1==s2), error('dimension mismatch'), end
s = s1; s(s1==1)=s2(s1==1);

u = repmat(u,s-s1+1);
v = repmat(v,s-s2+1);
y = u+v;
    
% A FRUITLESS ATTEMPT TO OPTIMIZE THE CODE IN THE CASE OF LARGE ARRAYS:
% % sizes
% s1 = size(u); 
% s2 = size(v);
% n = max(length(s1),length(s2));
% s1(end+1:n) = 1;
% s2(end+1:n) = 1;
% if ~all(s1==1 | s2==1 | s1==s2), error('dimension mismatch'), end
% s = s1; s(s1==1)=s2(s1==1);
% 
% % check memory "consumption"
% type = class(u([])+v([]));
% test = zeros(1,type);
% test = whos('test');
% nbytes = test.bytes * prod(s);
% 
% % use a for loop only if data is big
% if nbytes < 5e8
%     u = repmat(u,s-s1+1);
%     v = repmat(v,s-s2+1);
%     y = u+v;
% else
%     % we want u to be the "big guy"
%     if prod(s1)<prod(s2), y = fn_add(v,u); return, end 
%     y = repmat(u,s-s1+1);
%     % now perform the loop it the way that minimizes the number of
%     % iterations!
%     if prod(s2)<(prod(s)/prod(s2))
%         subs = struct('type','()','subs',repmat({':'},1,n));
%         dims = find(s2>1);
%         for k=1:prod(s2)
%             ind = fn_indices(s2,k);
%             for i=dims, subs.subs{i}=ind(i); end
%             y = su
%         end
%     else
%     end
% end
    
    