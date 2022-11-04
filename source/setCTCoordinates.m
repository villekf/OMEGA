function options = setCTCoordinates(options)
%SETCTCOORDINATES Computes the CT source-detector coordinates
%   Utility function
if options.subsets > 1 && options.subset_type < 8
    [options.x,options.y,options.z] = CTDetectorCoordinatesFull(options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
        options.ySize,options.horizontalOffset,options.verticalOffset,options.bedOffset);
end

if ~isfield(options,'x') && ~isfield(options,'y') && ~isfield(options,'z') && ~isfield(options,'z_det')
    [options.x,options.y,options.z] = CTDetSource(options.angles,options.sourceToDetector,options.sourceToCRot,...
        options.horizontalOffset,options.verticalOffset,options.bedOffset, options.uCenter, options.vCenter);
    if numel(options.z)/2 > numel(options.angles)
        if size(options.angles,1) == 1
            options.angles = reshape(options.angles, [],1);
        end
        options.angles = repmat(options.angles,numel(options.z)/2/numel(options.angles),1);
    end
end
if ~isfield(options, 'uV') && numel(options.x) ~= options.xSize * options.ySize * options.nProjections
    options.uV = CTDetectorCoordinates(options.angles,options.pitchRoll);
    if options.implementation == 3 || options.implementation == 2
        options.uV = single(options.uV);
    end
end
if isfield(options,'z_det')
    options.z = options.z_det;
    options = rmfield(options,'z_det');
end
if isfield(options, 'uV')
    options.x = [options.x(:,1) options.y(:,1) options.z(:,1) options.x(:,2) options.y(:,2) options.z(:,2)]';
    options.z = options.uV;
end
if (isfield(options,'pitchRoll') && ~isempty(options.pitchRoll)) && isfield(options, 'uV')
    options.PITCH = true;
    options.z(1,:) = options.z(1,:) * options.dPitchX;
    options.z(2,:) = options.z(2,:) * options.dPitchX;
    options.z(4,:) = options.z(4,:) * options.dPitchX;
    options.z(5,:) = options.z(5,:) * options.dPitchX;
    options.z(3,:) = options.z(3,:) * options.dPitchY;
    options.z(6,:) = options.z(6,:) * options.dPitchY;
elseif isfield(options, 'uV')
    options.z(1,:) = options.z(1,:) * options.dPitchX;
    options.z(2,:) = options.z(2,:) * options.dPitchX;
end
end