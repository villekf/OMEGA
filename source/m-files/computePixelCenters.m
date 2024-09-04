function [x_center,y_center,z_center] = computePixelCenters(xx,yy,zz,dx,dy,dz,TOF,options)
%COMPUTEPIXELCENTERS Computes the coordinates for the voxel centers
%   Utility function
% Coordinates of the centers of the voxels, mainly for projector types 2
% and 3
if options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33 ...
        || options.projector_type == 12 || options.projector_type == 13 || options.projector_type == 21 || options.projector_type == 31 ...
        || options.projector_type == 42 || options.projector_type == 43 || options.projector_type == 24 || options.projector_type == 34
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
    z_center = zz(1 : end - 1)' + dz/2;
elseif (options.projector_type == 1 && TOF)
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
else
    x_center = xx(1);
    y_center = yy(1);
    if isempty(zz)
        z_center = 0;
    else
        z_center = zz(1);
    end
end
end

