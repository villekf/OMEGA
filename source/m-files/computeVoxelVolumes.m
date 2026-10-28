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
    % Aligned look-up table: a regular grid with the spacing 1/CC = 1e-3 starting at bmin (the kernel reads the index
    % round((b - bmin) * CC)), long enough for the index of bmax. No entries are dropped, so that the table index always
    % corresponds to the same distance as in the kernel
    bmin = max(0, floor((options.tube_radius - options.voxel_radius) * 1e3) / 1e3);
    nV = ceil((bmax - bmin) * 1e3) + 3;
    b = bmin + (0:nV-1)' * 1e-3;
    V = zeros(nV, 1);
    inside = b < bmax;
    V(inside) = volumeIntersection(options.tube_radius, options.voxel_radius, b(inside));
    % b == tube radius is a removable singularity of the closed form (NaN), interpolate over it
    bad = inside & ~isfinite(V);
    if any(bad)
        V(bad) = interp1(b(inside & ~bad), V(inside & ~bad), b(bad), 'linear');
    end
    % The volume cannot increase with the distance
    V = cummin(abs(V));
    Vmax = (4*pi)/3*options.voxel_radius^3;
else
    V = 0;
    Vmax = 0;
    bmin = 0;
    bmax = 0;
end

V = cast(V, options.cType);
end