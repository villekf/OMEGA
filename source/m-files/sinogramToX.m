function returnList = sinogramToX(angles, radii, nx, ny, crxy)
    arguments (Input)
        angles  (:, 1) double % Column vector of panel angles
        radii   (:, 1) double % Column vector of panel radii
        nx      (1, 1) double % Detector column count
        ny      (1, 1) double % Detector row count
        crxy    (1, 1) double % Detector crystal pitch
    end

    if (numel(angles) ~= numel(radii))
        error('Different amount of angles and radii')
    end

    nIter = numel(angles);
    returnList = zeros(6, nIter * nx * ny);
    
    panelXmin = -crxy * (nx - 1) / 2;
    panelXmax = -panelXmin;
    panelYmin = -crxy * (ny - 1) / 2;
    panelYmax = -panelYmin;
    
    % Detector x points
    x = zeros(ny, nx);
    
    % Detector y points
    y = repmat(linspace(panelYmin, panelYmax, ny)', [1, nx]);
    
    % Detector z points
    z = repmat(linspace(panelXmin, panelXmax, nx), [ny, 1]);
    
    % Rotate and move
    idxCounter = 1;
    for nn = 1:nIter
        ang = angles(nn);
        rr = radii(nn);
        
        nVec = rr*[cosd(ang); sind(ang)]; % Panel normal
        R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)]; % Rotation matrix
        for ii = 1:nx
            for jj = 1:ny
                detXY = R * [x(jj, ii); y(jj, ii)] + nVec;
                detZ = z(jj, ii);
    
                returnList(:, idxCounter) = [detXY + nVec; detZ; detXY; detZ];
                
                idxCounter = idxCounter + 1;
            end
        end
    end
end