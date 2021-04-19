function [V,Vmax,bmin,bmax] = computeVoxelVolumes(dx,dy,dz,options)
%COMPUTEVOXELVOLUMES Computes the voxel volumes for the spherical voxels in
%THOR algorithm
%   Utility function
% Compute the various volumes available for the spherical voxels in
% volume-based ray tracer. Basically establishes a look-up-table.
if options.projector_type == 3
    dp = max([dx,dy,dz]);
    options.voxel_radius = sqrt(2) * options.voxel_radius * (dp / 2);
    bmax = options.tube_radius + options.voxel_radius;
    b = linspace(0, bmax, 10000)';
    b(options.tube_radius > (b + options.voxel_radius)) = [];
    b = unique(round(b*10^3)/10^3);
    V = volumeIntersection(options.tube_radius, options.voxel_radius, b);
    diffis = [diff(V);0];
    b = b(diffis <= 0);
    V = abs(V(diffis <= 0));
    Vmax = (4*pi)/3*options.voxel_radius^3;
    bmin = min(b);
else
    V = 0;
    Vmax = 0;
    bmin = 0;
    bmax = 0;
end
if options.implementation == 2 || options.implementation == 3
    V = single(V);
    Vmax = single(Vmax);
    bmin = single(bmin);
    bmax = single(bmax);
else
    V = double(V);
    Vmax = double(Vmax);
    bmin = double(bmin);
    bmax = double(bmax);
end
end

