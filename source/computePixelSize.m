function [xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(FOV, N, offset, implementation)
%COMPUTEPIXELSIZE Computes the pixel size and distance from origin
%   Utility function for OMEGA

% Pixel boundaries
etaisyys = -(FOV) / 2 + offset;
if implementation == 2 || implementation == 3 || implementation == 5
    xx = single(linspace(etaisyys(1) + offset(1), -etaisyys(1) + offset(1), N(1) + 1));
    yy = single(linspace(etaisyys(2) + offset(2), -etaisyys(2) + offset(2), N(2) + 1));
    zz = single(linspace(etaisyys(3) + offset(3), -etaisyys(3) + offset(3), N(3) + 1));
else
    xx = double(linspace(etaisyys(1) + offset(1), -etaisyys(1) + offset(1), N(2) + 1));
    yy = double(linspace(etaisyys(2) + offset(2), -etaisyys(2) + offset(2), N(1) + 1));
    zz = double(linspace(etaisyys(3) + offset(3), -etaisyys(3) + offset(3), N(3) + 1));
end

% Distance of adjacent pixels
dx = diff(xx(1:2));
dy = diff(yy(1:2));
dz = diff(zz(1:2));

% Distance of image from the origin
bx = xx(1);
by = yy(1);
bz = zz(1);
end

