function fn_plotscale(tlabel,ylabel)
% function fn_plotscale(tlabel,ylabel)
%---
% Creates a 2 directional scale bar for graph display
%
% Input:
% - tlabel      string with %f%s format, for example '1s'
% - ylabel      string with %f%s format, for example '10mV'
%
% See also fn_scale

% Thomas Deneux
% Copyright 2007-2012

if nargin==0, help fn_plotscale, return, end
    
set(gca,'visible','off')
delete(findobj(gca,'tag','fn_plotscale'))

tsize = sscanf(tlabel,'%f');
ysize = sscanf(ylabel,'%f');

% axis tight
ax = axis; 
posa = get(gca,'position'); posf = get(gcf,'position');

% scale in the bottom left
axis([ax(1) ax(2) ax(3)-2*ysize ax(4)])
ax = axis; fact = [ax(2)/(posa(3)*posf(3)) ax(4)/(posa(4)*posf(4))];
orig = ax([1 3])+fact;

% % scale in the left and vertically centered
% axis([ax(1)-2*tsize ax(2) ax(3) ax(4)])
% ax = axis; fact = [ax(2)/(posa(3)*posf(3)) ax(4)/(posa(4)*posf(4))];
% orig = [ax(1) (ax(3)+ax(4)-ysize)/2]+fact;


line(orig(1)+[0 tsize],orig(2)+[0 0],'color','black','linewidth',3, ...
    'tag','fn_plotscale')
text(orig(1)+tsize/2,orig(2)-2*fact(2),tlabel, ...
    'horizontalalignment','center','verticalalignment','top', ...
    'tag','fn_plotscale')
line(orig(1)+[0 0],orig(2)+[0 ysize],'color','black','linewidth',3, ...
    'tag','fn_plotscale')
text(orig(1)-2*fact(1),orig(2)+ysize/2,ylabel, ...
    'horizontalalignment','right','verticalalignment','middle', ...
    'tag','fn_plotscale')

