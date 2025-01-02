function [x, z] = get_coordinates_SPECT(options)
    % Error check
    n1 = numel(options.angles);
    n2 = numel(options.radiusPerProj);
    n3 = numel(options.swivelAngles);
    assert((n1 == n2) && (n2 == n3), 'The amount of angles, radii and swivel angles have to be equal.')

    % Variable initialization
    nProjections = n1;
    x = zeros(6, nProjections); % Normal vectors of each projection
    z = zeros(2, nProjections); % Orthogonal to x, |z|=detector pixel size

    % Handle each projection
    for ii = 1:nProjections
        r1 = options.radiusPerProj(ii);
        r2 = options.CORtoDetectorSurface;

        % The (gantry angle + home angle) should be relative to positive x axis with positive rotation direction counterclockwise
        % The swivel angle should also be relative to positive x axis with positive rotation direction counterclockwise
        alpha1 = options.angles(ii); 
        alpha2 = options.homeAngles(ii) - 360/(2*pi)*options.offangle;
        alpha3 = options.swivelAngles(ii) - 360/(2*pi)*options.offangle; 

        if ~options.flip_image
            alpha2 = alpha2 + 180;
            alpha3 = alpha3 + 180;
        end
        
        x(4, ii) = r1 * cosd(alpha1+alpha2) + r2 * cosd(alpha3);
        x(5, ii) = r1 * sind(alpha1+alpha2) + r2 * sind(alpha3);
        x(6, ii) = 0;
        x(1, ii) = x(4, ii) + options.colL * cosd(alpha3);
        x(2, ii) = x(5, ii) + options.colL * sind(alpha3);
        x(3, ii) = 0;
        
        z(1, ii) = options.crXY * cosd(alpha3-90);
        z(2, ii) = options.crXY * sind(alpha3-90);
    end
    x = cast(x, options.cType);
    z = cast(z, options.cType);
end