function hCirc = drawCirclesvenn(xc, yc, r, c, fa, tag)
% Draw circles for Venn Diagram
% given centre coordinates xc and yc, radius r, color c, 
% transparency alpha fa, and label tag
% takes vectors to draw multiple circles
% adapted from: https://github.com/brian-lau/matutils/blob/master/%2Bfig/venn.m
% see copyright notice.

    
    %P and V are cell arrays of patch parameter/values
    xc = xc(:); yc = yc(:);     %Circle centers
    r = r(:);                   %Radii
    n = length(r);              
    
    %Independent parameter
    dt = 0.05;
    t = 0:dt:2*pi;

    %Origin centered circle coordinates
    X = r*cos(t);
    Y = r*sin(t);
    
    hCirc = zeros(1,n);

    for i = 1:n
        xx = X(i,:)+xc(i);  
        yy = Y(i,:)+yc(i);
        hCirc(i) = patch (xx, yy, c{i}, 'FaceAlpha', fa{i});
        text(xc(i)+1.5*r(i),yc(i),tag{i});

    end

    
end %plotCircles