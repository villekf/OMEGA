function [V,Vmax,bmin,bmax] = computeVoxelVolumes(dx,dy,dz,options)
%COMPUTEVOXELVOLUMES Computes the voxel volumes for the spherical voxels in
%THOR algorithm
%   Utility function
% Compute the various volumes available for the spherical voxels in
% volume-based ray tracer. Basically establishes a look-up-table.
if options.projector_type == 3 || options.projector_type == 33 || options.projector_type == 13 || options.projector_type == 31 ...
        || options.projector_type == 34 || options.projector_type == 43
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

V = cast(V, options.cType);
end