function boxplotx(data,varargin)
% boxplotx  Plot boxplots with whiskers
%  2016-03-20  Matlab  Copyright (c) 2016, W J Whiten  BSD License
%
% boxplot2(points,xpos,width,optn)
%  data     data for boxplot matrix of columns or cell array of vectors
%  xpos     Position vector for boxplots  (optional)
%  width    Width vector for boxplots (optional after xpos)
%  optn     Optional struct or name value pairs, see optndfts
%            .xpos   Position vector for boxplots (default 1:n)
%            .width  Width vector for boxplots (default 0.5)
%            .tightx Set x limits to tight fit (default true)
%            .lines  Code or cells for boxplot lines (default 'b-')
%            .points Code or cells for boxplot points (deflaut 'k+')
%            .labels Cell vector of labels (default none)
%            .rotation Rotation for X labels (default 60)
%            .qpct   Box upper and lower percentile (default 25)
%            .bpct   Bar upper & lower percentile based on normal 
%                     distribution and qfrct (default 1)
%
% Simple useage:
%  boxplot(data)  % boxplots from columns of data
%
%  Output is to plot:
%
% Top bar        Value of largest point less than estimated 
%                 100-bpct percentile value
% Top of box     100-qpct percentile (default 3rd quartile)
% Middle of box  Median value
% Bottom of box  qpct percentile (default 1st quartile)
% Bottom bar     Value of smallest point less than estimated 
%                 bpct percentile value
%
% Median, 100-qpct percentile (default 3rd quartile) and 
%  qpct percentile (default 1st quartile) are interpolated from data.
% 100-bpct percentile is estimated from median and 100-qpct percentile
%  using the Gausian (normal) distribution
% bpct percentile is estimated from median and qpct percentile
%  using the Gausian (normal) distribution

if(size(data,1)==1)
    data=data(:);
end
if(isreal(data))
    n=size(data,2);
else
    n=length(data);
end

[optn,optn2]=optndfts(varargin,{'xpos','width'},'xpos',1:n,  ...
    'width',NaN,'lines','b-','points','k+','tightx',true,  ...
    'labels',{},'rotation',60);

optn.xpos=optn.xpos(:)';

if(isnan(optn.width))
    if(n==1)
        optn.width=0.55;
    else
        optn.width=repmat((max(optn.xpos)-min(optn.xpos))/(n-1)*0.55,1,n);
    end
end
if(length(optn.width)==1)
    optn.width=repmat(optn.width,1,n);
end
if(ischar(optn.lines))
    optn.lines=repmat({optn.lines},1,n);
end
if(ischar(optn.points))
    optn.points=repmat({optn.points},1,n);
end

% check if adding to existing plot
if(ishold())
    userdata=get(gca,'UserData');
    if(isfield(userdata,'xpos'))
        userdata=userdata.xpos;
        t1=[userdata,optn.xpos];
        t2=unique(t1);
        if(length(t1)>length(t2))
            optn.xpos=1+max(userdata)-min(optn.xpos)+optn.xpos;
        end
        userdata=[userdata,optn.xpos];
    else
        userdata=optn.xpos;
    end
    xtick=get(gca,'xtick');
    xticklabel=get(gca,'xticklabel');
else
    userdata=optn.xpos;
    xtick=[];
    xticklabel={};
end

% repeat for each data column
for i=1:n
    if(iscell(data))
        datai=data{i};
    else
        datai=data(:,i);
    end
    width=optn.width(i);
    xpos=optn.xpos(i);

    [bars,outliers,pcts]=boxplotinfo(datai,optn2);
    if(isempty(bars))
        if(n==1)
            plot(NaN,NaN)
            warning('boxplotx:  No data')
        end
        continue
    end

    w2=width/2;
    w4=width/4;

    x1=xpos-w4;
    x2=xpos+w4;
    x3=xpos-w2;
    x4=xpos+w2;

    % horizontal bars
    xx=[x1,x2;x3,x4;x3,x4;x3,x4;x1,x2]';
    yy=[bars;bars];

    plot(xx,yy,optn.lines{i})
    hold on

    % vertical lines
    xx=[xpos,xpos;x3,x3;x4,x4;xpos,xpos]';
    yy=[bars(1),bars(2);bars(2),bars(4);bars(2),bars(4);bars(4),bars(5)]';

    plot(xx,yy,optn.lines{i})

    % additional points, spread duplicate points
    if(~isempty(outliers))
        xs=xspread(outliers,xpos,width/2,bars(5)-bars(1));
        plot(xs,outliers,optn.points{i})
    end
end

% set x axis range to close to range used
if(optn.tightx && n>1)
    x1=min(optn.xpos);
    x2=max(optn.xpos);
    xd=(x2-x1)/(n-1);
    xlim([x1-xd,x2+xd]);
end

% add labels if given
if(~isempty(optn.labels))
    [xtick,indx]=sort([xtick,optn.xpos(:)']);
    xticklabel=[xticklabel;optn.labels(:)];
    xticklabel=xticklabel(indx);
    set(gca,'xtick',xtick,'xticklabel',xticklabel,  ...
        'XTickLabelRotation',optn.rotation);    
end
    
% set x axis range to close to range used
if(optn.tightx && length(userdata)>1)
    x1=min(userdata);
    x2=max(userdata);
    xd=(x2-x1)/(length(userdata)-1);
    xlim([x1-xd,x2+xd]);
end

% save xpos for updating graph
set(gca,'UserData',struct('xpos',userdata));

ylabel(['Bars at  min>',num2str(pcts(1)),',  ',num2str(pcts(2)),  ...
    ',  ',num2str(pcts(3)),',  ',num2str(pcts(4)),',  max<',  ...
    num2str(pcts(5)),'  percentiles'])
hold off

return
end



function [bars,outliers,pcts]=boxplotinfo(data,varargin)
% boxplotinfo  Bar positions and extreme points for boxplot1
%  2016-03-20  Matlab  Copyright (c) 2016, W J Whiten  BSD License
%
% [bars,points]=boxplotinfo(data,optn)
%  data     Vector of data values for boxplot
%  optn     Optional struct or name value pairs, see optndfts
%            .qpct  Box upper and lower percentile (default 25)
%            .bpct   Bar upper & lower percentile based on normal 
%                     distribution and qfrct (default 1)
%
%  bars     Bar positions [low,25,50,75,high] or []
%  outliers Points to be plotted, those outside [low,high]
%  pcts     Percent values for horizontal bars

optn=optndfts(varargin,'qpct',25,'bpct',1);
bar=erfinv(optn.bpct/50-1)/erfinv(optn.qpct/50-1);

data=sort(data);
data=data(isfinite(data));

n=length(data);
n1=n-1;

if(n<=1)
    bars=repmat(data,1,5);
    outliers=data;
    pcts=[optn.bpct,optn.qpct,50,100-optn.qpct,100-optn.bpct];

    return
end

% median
n50a=floor(n1/2)+1;
n50b=ceil(n1/2)+1;
bars3=(data(n50a)+data(n50b))/2;

if(n==2)
    bars=repmat(bars3,1,5);
    outliers=data;
    pcts=[optn.bpct,optn.qpct,50,100-optn.qpct,100-optn.bpct];


    return
end

% lower quartile
n25=n1*optn.qpct/100+1;
n25a=floor(n25);
n25b=ceil(n25);
if(n25a==n25b)
    bars2=data(n25);
else
    bars2=data(n25a)*(n25b-n25)+data(n25b)*(n25-n25a);
end

% upper quartile
n75=n1*(1-optn.qpct/100)+1;
n75a=floor(n75);
n75b=ceil(n75);
if(n75a==n75b)
    bars4=data(n75);
else
    bars4=data(n75a)*(n75b-n75)+data(n75b)*(n75-n75a);
end

% lower bar and upper bar
bars1=data(sum(data<bars3-bar*(bars3-bars2))+1);
bars5=data(sum(data<bars3+bar*(bars4-bars3)));

% points outside bars
outliers=data(data<bars1 | data>bars5);

bars=[bars1,bars2,bars3,bars4,bars5];
pcts=[optn.bpct,optn.qpct,50,100-optn.qpct,100-optn.bpct];


return
end



function x=xspread(y,xpos,width,height)
% xspread  Spread x values for plots of nearly equal y values
%  20116-03-20 Matlab Copyright (c) 2016, W J Whiten  BSD License
%
% x=xspread(y,xpos,width)
%  y      Y values to be plotted in sorted order
%  xpos   X postion to plot
%  width  Width to spread over
%  height Hight of box plot
%
%  x     X values to be used in plotting y

% faction of range for equal values
ystep=0.01;

n=length(y);
x=repmat(xpos,size(y));

% step size for near equal values
dy=ystep*(max(y)-min(y)+height);

% find repeated values and adjust x values to separate
i=1;
i1=1;
while(i<=n)
    while(i<=n && abs(y(i)-y(i1))<=dy)
        i=i+1;
    end
    if(i-1>i1)
        id=i-i1;
        w1=width*(1-1/id);
        x(i1:i-1)=(0:id-1)/(id-1)*w1+xpos-w1/2;
    end
    i1=i;
end

return
end
