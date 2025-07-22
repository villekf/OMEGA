function [options] = computeProjectorScalingValues(options)
%COMPUTEPROJECTORSCALINGVALUES Computes scaling values required by
%projector types 4 and 5
%   Utility function, used by projector types 4 and 5 only

options.dScaleX4 = 1 ./ (options.dx .* cast(options.Nx, options.cType));
options.dScaleY4 = 1 ./ (options.dy .* cast(options.Ny, options.cType));
options.dScaleZ4 = 1 ./ (options.dz .* cast(options.Nz, options.cType));
if options.projector_type == 5 || options.projector_type == 15 || options.projector_type == 45 || options.projector_type == 54 || options.projector_type == 51
    options.dSizeY = 1 ./ (options.dy .* cast(options.Ny, options.cType));
    options.dSizeX = 1 ./ (options.dx .* cast(options.Nx, options.cType));
    options.dScaleX = 1 ./ (options.dx .* cast(options.Nx + 1, options.cType));
    options.dScaleY = 1 ./ (options.dy .* cast(options.Ny + 1, options.cType));
    options.dScaleZ = 1 ./ (options.dz .* cast(options.Nz + 1, options.cType));
    options.dSizeZBP = (options.nColsD + 1) * options.dPitchX;
    options.dSizeXBP = (options.nRowsD + 1) * options.dPitchY;
end
if options.projector_type == 4 || options.projector_type == 14 || options.projector_type == 54 || options.projector_type == 45 ...
         || options.projector_type == 42 || options.projector_type == 43 || options.projector_type == 24 || options.projector_type == 34
    options.kerroin = cast((options.dx .* options.dy .* options.dz) ./ (options.dPitchX * options.dPitchY * options.sourceToDetector), options.cType);
    if options.dL == 0
        options.dL = options.FOVa_x(1) / double(options.Nx(1));
    elseif options.dL ~= options.FOVa_x(1) / double(options.Nx(1))
        options.dL = options.dL * (options.FOVa_x(1) / double(options.Nx(1)));
    end
else
    options.kerroin = 0;
end
if (options.projector_type == 4 || options.projector_type == 5 || options.projector_type == 14 || options.projector_type == 15 || options.projector_type == 45 ...
        || options.projector_type == 54) && options.CT
    options.use_64bit_atomics = false;
    options.use_32bit_atomics = false;
end
end