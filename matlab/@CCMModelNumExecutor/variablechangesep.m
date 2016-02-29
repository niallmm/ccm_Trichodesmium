function [r,h,c] = variablechangesep(obj,u)
x = obj.x;
% now we need to convert u into two vectors h and c
[times,sizeu] = size(u);
h = zeros(times, sizeu/2);
c = h;
h=u(:,1:2:sizeu);
c=u(:,2:2:sizeu);

r = (3*x).^(1/3);