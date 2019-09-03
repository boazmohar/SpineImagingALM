function B = fn_map(fun,A,mode)
% function B = fn_map(fun,A,mode)
% 
% map function 'fun' 
% to elements(mode=0) / columns(mode=1)[default] / rows(mode=2) of A

% Thomas Deneux
% Copyright 2006-2012

if nargin<2, mode=1; end

if mode==0
    B=A;
    for i=1:prod(size(A)), B(i)=feval(fun,A(i)); end
elseif mode==1
    B=zeros(prod(size(feval(fun,A(:,1)))),size(A,2));
    for i=1:size(A,2), B(:,i)=feval(fun,A(:,i)); end
elseif mode==2
    B=zeros(size(A,1),prod(size(feval(fun,A(1,:)))));
    for i=1:size(A,1), B(i,:)=feval(fun,A(i,:)); end
end
    