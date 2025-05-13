function options = setCTCoordinates(options)
%SETCTCOORDINATES Computes the CT source-detector coordinates
%   Utility function
if options.subsets > 1 && options.subset_type < 8 && options.subset_type > 0
    [options.x,options.y,options.z] = CTDetectorCoordinatesFull(options.angles,options.nProjections,options.nColsD,options.nRowsD,options.sourceToDetector,options.sourceToCRot,options.dPitch,...
        options.sourceOffsetRow,options.sourceOffsetCol,options.bedOffset, options.detOffsetRow, options.detOffsetCol);
    if abs(options.offangle) > 0
        R = [cos(options.offangle) -sin(options.offangle); sin(options.offangle) cos(options.offangle)];
        det = [options.x(:,1) options.y(:,1)]';
        sXY = (R * det)';
        det = [options.x(:,2) options.y(:,2)]';
        dXY = (R * det)';
        options.x = [sXY(:,1) dXY(:,1)];
        options.y = [sXY(:,2) dXY(:,2)];
    end
end

if ~isfield(options,'x') && ~isfield(options,'y') && ~isfield(options,'z') && ~isfield(options,'z_det')
    [options.x,options.y,options.z] = CTDetSource(options.angles,options.nProjections, options.sourceToDetector,options.sourceToCRot,...
        options.sourceOffsetRow,options.sourceOffsetCol,options.bedOffset, options.detOffsetRow, options.detOffsetCol);
    if numel(options.z)/2 > numel(options.angles)
        if size(options.angles,1) == 1
            options.angles = reshape(options.angles, [],1);
        end
        options.angles = repmat(options.angles,numel(options.z)/2/numel(options.angles),1);
    end
elseif abs(options.offangle) > 0
    R = [cos(options.offangle) -sin(options.offangle); sin(options.offangle) cos(options.offangle)];
    det = [options.x(:,1) options.y(:,1)]';
    sXY = (R * det)';
    det = [options.x(:,2) options.y(:,2)]';
    dXY = (R * det)';
    options.x = [sXY(:,1) dXY(:,1)];
    options.y = [sXY(:,2) dXY(:,2)];
end
if (~isfield(options, 'uV') || (isfield(options,'uV') && isempty(options.uV))) && numel(options.x) ~= options.nColsD * options.nRowsD * options.nProjections
    options.uV = CTDetectorCoordinates(options.angles,options.pitchRoll);
    if options.implementation == 3 || options.implementation == 2 || options.implementation == 5 || options.useSingles
        options.uV = single(options.uV);
    end
end
if options.flip_image && isfield(options,'x')
    options.y = -options.y;
    if isfield(options, 'uV')
        options.uV(2,:) = -options.uV(2,:);
        if (isfield(options,'pitchRoll') && ~isempty(options.pitchRoll))
            options.uV(5,:) = -options.uV(5,:);
        end
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
if ((isfield(options,'pitchRoll') && ~isempty(options.pitchRoll)) && isfield(options, 'uV')) || (isfield(options, 'uV') && size(options.uV, 1) == 6)
    options.pitch = true;
    options.z(1,:) = options.z(1,:) * options.dPitchX;
    options.z(2,:) = options.z(2,:) * options.dPitchX;
    options.z(3,:) = options.z(3,:) * options.dPitchY;
    %     if options.projector_type == 5
    %         options.z(4,:) = 0;
    %         options.z(5,:) = 0;
    %         options.z(3,:) = 0;
    %     else
    options.z(4,:) = options.z(4,:) * options.dPitchX;
    options.z(5,:) = options.z(5,:) * options.dPitchX;
    options.z(6,:) = options.z(6,:) * options.dPitchY;
%     end
elseif isfield(options, 'uV')
    options.z(1,:) = options.z(1,:) * options.dPitchX;
    options.z(2,:) = options.z(2,:) * options.dPitchX;
end
end