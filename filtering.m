function [output] = filtering(x,y)
%FILTERING Summary of this function goes here
%   Detailed explanation goes here
xmin = min(x);
xmax = max(x);
n = 4*length(x);
xorder = linspace(xmin,xmax,n);
dx = (xmax-xmin)/(n-1);
yorder = interp1(x,y,xorder);
xorder = [xorder+xmin-xmax-dx,xorder,xorder-xmin+xmax+dx];
yorder = [yorder,yorder,yorder];
[B,A] = butter(4,0.2,'low');
yorder = filter(B,A,yorder);
output = interp1(xorder,yorder,x);
end

