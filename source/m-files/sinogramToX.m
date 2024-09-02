function returnList = sinogramToX(angles, radii, nx, ny, crxy)
%SINOGRAMTOX A function for converting sinogram pixel position and orientation information to a list of coordinates. Each column of the output has 6 elements, two pairs of xyz points. The line spanned by the points corresponds to the detector pixel normal vector.
%   Utility function for OMEGA

    if (numel(angles) ~= numel(radii))
        error('Different amount of angles and radii')
    end

    nIter = numel(angles);
    returnList = zeros(6, nIter * ny * nx);

    panelXmin = -crxy * (ny - 1) / 2;
    panelXmax = -panelXmin;
    panelYmin = -crxy * (nx - 1) / 2;
    panelYmax = -panelYmin;

    % Detector x points
    x = zeros(nx, ny);

    % Detector y points
    y = repmat(linspace(panelYmin, panelYmax, nx)', [1, ny]);

    % Detector z points
    z = repmat(linspace(panelXmin, panelXmax, ny), [nx, 1]);

    % Rotate and move
    idxCounter = 1;
    for nn = 1:nIter
        ang = angles(nn);
        rr = radii(nn);

        nVec = rr*[cosd(ang); sind(ang)]; % Panel normal
        R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)]; % Rotation matrix
        for ii = 1:ny
            for jj = 1:nx
                detXY = R * [x(jj, ii); y(jj, ii)] + nVec;
                detZ = z(jj, ii);

                returnList(:, idxCounter) = [detXY + nVec; detZ; detXY; detZ];

                idxCounter = idxCounter + 1;
            end
        end
    end
end
