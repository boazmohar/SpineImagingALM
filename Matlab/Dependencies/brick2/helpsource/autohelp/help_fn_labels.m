%% fn_labels

%% Syntax
%  fn_labels(x_label,y_label[,legend_[,title[,figure name]])

%% Description
%  Shortcut for:
%  xlabel(x_label)
%  ylabel(y_label)
%  fn_nicegraph
%  if nargin>2, legend(legend_{:},'location','northwest'), end
%  if nargin>3, title(title_), end
%  if nargin>4, set(gcf,'name',name_), end
% 
%  See also fn_nicegraph, fn_scale, fn_plotscale

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
