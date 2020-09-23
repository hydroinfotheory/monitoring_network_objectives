function draw_arc(S,oncircle,bycircle,LR,long)
% draw_arc(S,oncircle,bycircle,LR,long)
% Draws circle arcs between intersections of two circles 
% Inputs:
% S= struct output from function venn that draws venn diagrams
% oncircle = nr of circle along which to draw
% bycircle = nr of circle whose intersections to use
% LR = 0 or 1, flip
% long 0 or 1, draw long rout between intersections?


hold on;
x1=S.Position(oncircle,1);
x2=S.Position(bycircle,1);
y1=S.Position(oncircle,2);
y2=S.Position(bycircle,2);
r1=S.Radius(oncircle);
r2=S.Radius(bycircle);
r=r1;
% get intersection of 2 circles
[xout,yout] = circcirc(x1,y1,r1,x2,y2,r2)
x1=xout(1);
x2=xout(2);
y1=yout(1);
y2=yout(2);
%calc stuff for arcs
d = sqrt((x2-x1)^2+(y2-y1)^2); % Distance between points
if LR==0
    a = atan2(-(x2-x1),y2-y1); % Perpendicular bisector angle
else
    a = atan2(x2-x1,-(y2-y1));
end
b = asin(d/2/r); % Half arc angle
if long
    b=pi-b;
end
c = linspace(a-b,a+b,1000); % Arc angle range
e = sqrt(r^2-d^2/4); % Distance, center to midpoint
if long
    e=-e;
end
x = (x1+x2)/2-e*cos(a)+r*cos(c); % Cartesian coords. of arc
y = (y1+y2)/2-e*sin(a)+r*sin(c);
plot(x,y,'r-','linewidth',3)