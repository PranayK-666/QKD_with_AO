function [cx, cy] = circle_coord(x,y,r)

    th = 0:pi/50:2*pi;
    cx = r*cos(th) + x;
    cy = r*sin(th) + y;

end