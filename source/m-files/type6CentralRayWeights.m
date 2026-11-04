function weights = type6CentralRayWeights(options, views)
%TYPE6CENTRALRAYWEIGHTS Voxel-volume / clipped central-ray length per detector bin.
% Uses the same panel geometry and elliptical-cylinder or box support as
% type6CentralRayPathWeights in the native reconstruction projector.
rows = double(options.nRowsD);
cols = double(options.nColsD);
rowOffset = ((0:rows-1)' - (rows-1)/2) * double(options.dPitchX);
axialOffset = ((0:cols-1) - (cols-1)/2) * double(options.dPitchY);
center = double([options.ellipseCenterX, options.ellipseCenterY, options.ellipseCenterZ]);
radius = double([options.ellipseRadiusX, options.ellipseRadiusY, options.ellipseRadiusZ]);
if any(~isfinite(radius)) || any(radius <= 0)
    error('type6CentralRayWeights:InvalidSupport', 'Type 6 support radii must be finite and positive.');
end
power = double(options.ellipsePower);
if ~(power == 2 || (isinf(power) && power > 0))
    error('type6CentralRayWeights:InvalidSupport', 'Type 6 support must be an ellipse or box.');
end
voxelVolume = double(options.dx(1)) * double(options.dy(1)) * double(options.dz(1));
weights = zeros(rows, cols, numel(views), options.cType);
for kk = 1:numel(views)
    view = views(kk);
    angle = double(options.swivelAngles(view));
    nx = cosd(angle); ny = sind(angle);
    shift = double(options.blurPlanes2F(view)) * double(options.dx(1));
    qx = -ny * (rowOffset + shift);
    qy = nx * (rowOffset + shift);
    valid = repmat(abs(axialOffset - center(3)) <= radius(3), rows, 1);
    sourceN = double(options.radiusPerProj(view)) * cosd(double(options.angles(view)) - angle) + double(options.CORtoDetectorSurface);
    if isinf(power)
        lo = -inf(rows, 1); hi = inf(rows, 1);
        for axis = 1:2
            if axis == 1, direction = nx; q = qx; else, direction = ny; q = qy; end
            if abs(direction) < 1e-8
                valid = valid & (q >= center(axis)-radius(axis)) & (q <= center(axis)+radius(axis));
            else
                a = (center(axis)-radius(axis)-q) / direction;
                b = (center(axis)+radius(axis)-q) / direction;
                lo = max(lo, min(a,b)); hi = min(hi, max(a,b));
            end
        end
        chord = max(hi - max(lo, sourceN), 0);
    else
        px = qx-center(1); py = qy-center(2);
        invx = 1/radius(1)^2; invy = 1/radius(2)^2;
        a = nx^2*invx + ny^2*invy;
        b = 2*(px*nx*invx + py*ny*invy);
        c = px.^2*invx + py.^2*invy - 1;
        disc = b.^2 - 4*a*c;
        root = sqrt(max(disc, 0));
        exitT = (-b + root)/(2*a); entryT = (-b - root)/(2*a);
        chord = max(exitT - max(entryT, sourceN), 0);
        valid = valid & (disc > 1e-12);
    end
    valid = valid & (chord > 1e-6);
    plane = voxelVolume ./ max(chord, 1e-6) .* ones(1, cols);
    plane(~valid) = 0;
    weights(:,:,kk) = cast(plane, options.cType);
end
end
