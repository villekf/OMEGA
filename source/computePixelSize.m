function [xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(R, FOVax, FOVay, Z, axial_fov, Nx, Ny, Nz, implementation)
%COMPUTEPIXELSIZE Computes the pixel size and distance from origin
%   Utility function for OMEGA

% Pixel boundaries
etaisyys_x = (R - FOVax) / 2;
etaisyys_y = (R - FOVay) / 2;
etaisyys_z = (Z - axial_fov) / 2;
if implementation == 2 || implementation == 3 || implementation == 5
    zz = single(linspace(etaisyys_z, Z - etaisyys_z, Nz + 1));
    xx = single(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = single(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
else
    zz = double(linspace(etaisyys_z, Z - etaisyys_z, Nz + 1));
    xx = double(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = double(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
end
%     zz = zz(2*block1 + 1 : 2*blocks + 2);

% Distance of adjacent pixels
dx = diff(xx(1:2));
dy = diff(yy(1:2));
dz = diff(zz(1:2));

% Distance of image from the origin
bx = xx(1);
by = yy(1);
bz = zz(1);
end

