function [x_center,y_center,z_center] = computePixelCenters(xx,yy,zz,dx,dy,dz,TOF,options)
%COMPUTEPIXELCENTERS Computes the coordinates for the voxel centers
%   Utility function
% Coordinates of the centers of the voxels
if options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33 ...
        || options.projector_type == 12 || options.projector_type == 13 || options.projector_type == 21 || options.projector_type == 31 ...
        || options.projector_type == 42 || options.projector_type == 43 || options.projector_type == 24 || options.projector_type == 34
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
    z_center = zz(1 : end - 1)' + dz/2;
%     temppi = min([options.FOVa_x / options.Nx, options.axial_fov / options.Nz]);
%     if options.tube_width_z > 0
%         temppi = max([1,round(options.tube_width_z / temppi)]);
%     else
%         temppi = max([1,round(options.tube_width_xy / temppi)]);
%     end
%     temppi = temppi * temppi * 4;
%     if options.apply_acceleration || options.implementation == 4
%         if options.tube_width_z == 0
%             dec = uint32(sqrt(options.Nx^2 + options.Ny^2) * temppi);
%         else
%             dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * temppi);
%         end
%     else
%         dec = uint32(0);
%     end
elseif (options.projector_type == 1 && TOF)
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
%     if (options.apply_acceleration || options.implementation == 4) && options.n_rays_transaxial * options.n_rays_axial == 1
%         dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * 2);
%     else
%         dec = uint32(0);
%     end
else
    x_center = xx(1);
    y_center = yy(1);
    if isempty(zz)
        z_center = 0;
    else
        z_center = zz(1);
    end
%     if (options.apply_acceleration) && options.n_rays_transaxial * options.n_rays_axial == 1
%         dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * 2);
%     else
%         dec = uint32(0);
%     end
end
end

