function varargout = type6NativeOperator(action, varargin)
%TYPE6NATIVEOPERATOR Shared numerical primitives for the MATLAB type-6 path.
%   This file deliberately contains one public dispatcher so the primitives
%   can be reused by forwardProjection and backwardProjection.
%
%   type6NativeOperator('gather', image, angleDeg, tx, ty[, geom]) applies a
%   cropped bilinear rotation/translation. Optional geom fields are inputSize,
%   outputSize, inputSpacing, outputSpacing, inputOffset, outputOffset, and
%   translation (physical [x,y] mm). The adjoint action takes the exact same
%   arguments and scatters to the input grid.
%   type6NativeOperator('type6-rotate',image,angleDeg,offX,offY,geom) matches
%   the native rotateCustomAF grid mapping; 'type6-rotate-adjoint' scatters its
%   exact transpose. Offsets are in native voxel coordinates.
%   type6NativeOperator('conv-adjoint', image, kernel) applies the transpose
%   of conv2(image,kernel,'same'), including MATLAB's even-kernel alignment.
%   'grid-to-detector' gathers at the native affine detector map;
%   'detector-to-grid-adjoint' scatters the identical tent weights.

switch lower(char(action))
    case {'gather', 'adjoint'}
        if numel(varargin) < 4 || numel(varargin) > 5
            error('type6NativeOperator:AffineArguments', ...
                'Affine actions require image, angleDeg, tx, ty, and optional geometry.');
        end
        image = varargin{1};
        angleDeg = double(varargin{2});
        tx = double(varargin{3});
        ty = double(varargin{4});
        geom = struct();
        if numel(varargin) == 5, geom = varargin{5}; end
        if strcmpi(action, 'gather')
            varargout{1} = bilinearGather(image, angleDeg, tx, ty, geom);
        else
            varargout{1} = bilinearAdjoint(image, angleDeg, tx, ty, geom);
        end
    case {'rotation-parameters', 'rotation_parameters'}
        if numel(varargin) ~= 7, error('type6NativeOperator:RotationParameterArguments', 'rotation-parameters requires options, volume, view, gridSize, gridNz, padLo and padHi.'); end
        [varargout{1},varargout{2},varargout{3}] = nativeRotationParameters(varargin{:});
    case {'type6-rotate', 'type6_rotate'}
        if numel(varargin) ~= 5, error('type6NativeOperator:NativeRotationArguments', 'type6-rotate requires image, angleDeg, offX, offY, and geometry.'); end
        varargout{1} = nativeRotationGather(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    case {'type6-rotate-adjoint', 'type6_rotate_adjoint'}
        if numel(varargin) ~= 5, error('type6NativeOperator:NativeRotationArguments', 'type6-rotate-adjoint requires image, angleDeg, offX, offY, and geometry.'); end
        varargout{1} = nativeRotationAdjoint(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    case {'attenuation-factors', 'attenuation_factors'}
        if ~ismember(numel(varargin), [3, 6]), error('type6NativeOperator:AttenuationArguments', 'attenuation-factors requires options, volume, view, and optional grid size, grid Nz, and panel shift.'); end
        varargout{1} = attenuationFactors(varargin{1}, varargin{2}, varargin{3}, varargin(4:end));
    case {'cdrf', 'cdrf-kernel'}
        if numel(varargin) ~= 3, error('type6NativeOperator:CDRFArguments', 'cdrf requires options, one-based volume, and one-based view.'); end
        varargout{1} = cdrfKernel(varargin{1}, varargin{2}, varargin{3});
    case {'axial-geometry', 'axial_geometry'}
        if numel(varargin) ~= 2, error('type6NativeOperator:AxialArguments', 'axial-geometry requires options and a one-based volume index.'); end
        [varargout{1}, varargout{2}, varargout{3}] = axialGeometry(varargin{1}, varargin{2});
    case {'volume-geometry', 'volume_geometry'}
        if numel(varargin) ~= 2, error('type6NativeOperator:GridArguments', 'volume-geometry requires options and a one-based volume index.'); end
        [varargout{1}, varargout{2}] = volumeGeometry(varargin{1}, varargin{2});
    case {'grid-to-detector', 'grid_to_detector', 'grid-to-detector-adjoint', 'grid_to_detector_adjoint'}
        if numel(varargin) ~= 4, error('type6NativeOperator:ResampleArguments', 'grid-to-detector requires image, detectorRows, trA, and trB.'); end
        varargout{1} = gridToDetector(varargin{1}, varargin{2}, varargin{3}, varargin{4});
    case {'detector-to-grid-adjoint', 'detector_to_grid_adjoint', 'detector-to-grid', 'detector_to_grid'}
        if numel(varargin) ~= 4, error('type6NativeOperator:ResampleArguments', 'detector-to-grid-adjoint requires image, gridRows, trA, and trB.'); end
        varargout{1} = detectorToGridAdjoint(varargin{1}, varargin{2}, varargin{3}, varargin{4});
    case {'conv-adjoint', 'conv_adjoint'}
        if numel(varargin) ~= 2
            error('type6NativeOperator:ConvolutionArguments', ...
                'conv-adjoint requires image and kernel.');
        end
        varargout{1} = convolutionAdjoint(varargin{1}, varargin{2});
    otherwise
        error('type6NativeOperator:UnknownAction', 'Unknown action "%s".', char(action));
end
end

function [offX,offY,resample] = nativeRotationParameters(options,volume,view,gridSize,gridNz,padLo,padHi)
% offsetsFPType6/gridType6 offsets from precomp.h, in native voxel coordinates.
volume=double(volume); gridNx=double(gridSize(1)); gridNy=double(gridSize(2)); gridNz=double(gridNz);
centerX=double(options.bx(volume))+0.5*double(options.Nx(volume))*double(options.dx(volume));
centerY=double(options.by(volume))+0.5*double(options.Ny(volume))*double(options.dy(volume));
centerZ=double(options.bz(volume))+0.5*double(options.Nz(volume))*double(options.dz(volume));
mainCentered=abs(double(options.bx(1))+0.5*double(options.Nx(1))*double(options.dx(1)))<1e-3*double(options.dx(1)) && ...
    abs(double(options.by(1))+0.5*double(options.Ny(1))*double(options.dy(1)))<1e-3*double(options.dy(1)) && ...
    abs(double(options.bz(1))+0.5*double(options.Nz(1))*double(options.dz(1)))<1e-3*double(options.dz(1));
resample=volume>1 || ~mainCentered || padLo~=padHi || double(options.nRowsD)>gridNy || double(options.nColsD)>gridNz;
panel=double(options.blurPlanes2F(view))*double(options.dx(1))/double(options.dx(volume));
if ~resample && mod(gridNy-double(options.nRowsD),2)~=0, panel=panel+0.5; end
offX=-sind(double(options.swivelAngles(view)))*panel-centerX/double(options.dx(volume));
offY= cosd(double(options.swivelAngles(view)))*panel-centerY/double(options.dy(volume));
end

function output = nativeRotationGather(input, angleDeg, offX, offY, geom)
% Exact host analogue of C++ rotateCustomAF/auxKernels.cl rotate coordinates.
[srcNx,srcNy,nz]=size3(input); outNx=double(geom.outputSize(1)); outNy=double(geom.outputSize(2));
[rowSrc,colSrc,valid]=nativeRotationCoordinates(srcNx,srcNy,outNx,outNy,angleDeg,offX,offY);
[r0,c0,wr,wc]=zeroPaddedStencil(rowSrc,colSrc,valid,srcNx,srcNy);
output=zeros(outNx,outNy,nz,'like',input);
for z=1:nz
    output(:,:,z)=addGather(zeros(outNx,outNy,'like',input),input(:,:,z),r0,c0,wr,wc,srcNx,srcNy);
end
end

function output = nativeRotationAdjoint(input, angleDeg, offX, offY, geom)
% Scatter transpose of nativeRotationGather. The same sampling coordinates,
% zero-outside mask and edge clamping are used by the forward gather.
outNx=double(geom.outputSize(1)); outNy=double(geom.outputSize(2));
srcNx=double(geom.inputSize(1)); srcNy=double(geom.inputSize(2));
[~,~,nPlanes]=size3(input);
[rowSrc,colSrc,valid]=nativeRotationCoordinates(srcNx,srcNy,outNx,outNy,angleDeg,offX,offY);
[r0,c0,wr,wc]=zeroPaddedStencil(rowSrc,colSrc,valid,srcNx,srcNy);
output=zeros(srcNx,srcNy,nPlanes,'like',input);
for z=1:nPlanes
    output(:,:,z)=addScatter(zeros(srcNx,srcNy,'like',input),input(:,:,z),r0,c0,wr,wc,srcNx,srcNy);
end
end

function [rowSrc,colSrc,valid]=nativeRotationCoordinates(srcNx,srcNy,outNx,outNy,angleDeg,offX,offY)
[rowOut,colOut]=ndgrid(1:outNx,1:outNy);
qRow=rowOut-.5-outNx/2; qCol=colOut-.5-outNy/2;
c=cosd(-double(angleDeg)); s=sind(-double(angleDeg));
rowSrc=c.*qRow-s.*qCol+double(offX)+srcNx/2+.5;
colSrc=s.*qRow+c.*qCol+double(offY)+srcNy/2+.5;
% The OpenCL/CUDA kernel admits fractional edge coordinates whose floor or
% ceil neighbor is in range, then contributes zero for each out-of-range
% neighbor. In one-based coordinates that is the open interval (0,N+1).
valid=rowSrc>0 & rowSrc<srcNx+1 & colSrc>0 & colSrc<srcNy+1;
end

function [rows,cols,wr,wc]=zeroPaddedStencil(rowSrc,colSrc,valid,nRows,nCols)
r0=floor(rowSrc); c0=floor(colSrc); fr=rowSrc-r0; fc=colSrc-c0;
rows={r0, r0+1, r0, r0+1};
cols={c0, c0, c0+1, c0+1};
wr={1-fr,fr,1-fr,fr}; wc={1-fc,1-fc,fc,fc};
for k=1:4, wr{k}=wr{k}.*valid; end
end

function output = bilinearGather(input, angleDeg, tx, ty, geom)
% Gather from the inverse affine map. Coordinates are physical in x/y and
% output values outside the finite source grid are zero, with no renormalizing.
[nRows, nCols, nPlanes] = size3(input);
g = affineGeometry(geom, [nRows, nCols], [nRows, nCols], tx, ty);
[rowSrc, colSrc] = sourceCoordinates(nRows, nCols, g, angleDeg);
[r0, c0, wr, wc] = bilinearStencil(rowSrc, colSrc, nRows, nCols);
output = zeros([g.outputSize, nPlanes], 'like', input);
for z = 1:nPlanes
    output(:, :, z) = addGather(zeros(g.outputSize, 'like', input), input(:, :, z), ...
        r0, c0, wr, wc, nRows, nCols);
end
end

function output = bilinearAdjoint(input, angleDeg, tx, ty, geom)
% Exact transpose of bilinearGather; input is on the gather output lattice.
[nRowsOut, nColsOut, nPlanes] = size3(input);
g = affineGeometry(geom, [nRowsOut, nColsOut], [nRowsOut, nColsOut], tx, ty);
nRowsIn = g.inputSize(1);
nColsIn = g.inputSize(2);
[rowSrc, colSrc] = sourceCoordinates(nRowsIn, nColsIn, g, angleDeg);
[r0, c0, wr, wc] = bilinearStencil(rowSrc, colSrc, nRowsIn, nColsIn);
output = zeros([nRowsIn, nColsIn, nPlanes], 'like', input);
for z = 1:nPlanes
    output(:, :, z) = addScatter(zeros(nRowsIn, nColsIn, 'like', input), input(:, :, z), ...
        r0, c0, wr, wc, nRowsIn, nColsIn);
end
end

function g = affineGeometry(geom, defaultInputSize, defaultOutputSize, tx, ty)
g.inputSize = defaultInputSize;
g.outputSize = defaultOutputSize;
g.inputSpacing = [1, 1]; % [row, column] mm per pixel
g.outputSpacing = [1, 1];
g.inputOffset = [0, 0]; % physical [x,y] center of the input lattice
g.outputOffset = [0, 0];
g.translation = [tx, ty]; % legacy translations are in output pixels
if isempty(geom), return, end
fields = {'inputSize','outputSize','inputSpacing','outputSpacing','inputOffset','outputOffset','translation'};
for k = 1:numel(fields)
    if isfield(geom, fields{k}) && ~isempty(geom.(fields{k}))
        value = geom.(fields{k});
        g.(fields{k}) = double(value(:).');
    end
end
if ~isfield(geom, 'translation') || isempty(geom.translation)
    g.translation = [tx * g.outputSpacing(2), ty * g.outputSpacing(1)];
end
g.inputSize = round(g.inputSize(1:2));
g.outputSize = round(g.outputSize(1:2));
end

function [rowSrc, colSrc] = sourceCoordinates(nRowsIn, nColsIn, g, angleDeg)
[rowOut, colOut] = ndgrid(1:g.outputSize(1), 1:g.outputSize(2));
xOut = (colOut - (g.outputSize(2) + 1) / 2) .* g.outputSpacing(2) + g.outputOffset(1) - g.translation(1);
yOut = (rowOut - (g.outputSize(1) + 1) / 2) .* g.outputSpacing(1) + g.outputOffset(2) - g.translation(2);
c = cosd(angleDeg);
s = sind(angleDeg);
% Inverse of forward image-coordinate rotation [c s; -s c].
xSrc = c .* xOut - s .* yOut - g.inputOffset(1);
ySrc = s .* xOut + c .* yOut - g.inputOffset(2);
rowSrc = ySrc ./ g.inputSpacing(1) + (nRowsIn + 1) / 2;
colSrc = xSrc ./ g.inputSpacing(2) + (nColsIn + 1) / 2;
end

function [r0, c0, wr, wc] = bilinearStencil(rowSrc, colSrc, nRows, nCols)
r0 = floor(rowSrc);
c0 = floor(colSrc);
wr = rowSrc - r0;
wc = colSrc - c0;
% Neighbor offsets and weights are handled by the caller through arrays.
% Keep dimensions as matrices to avoid a large four-dimensional stencil.
r0 = {r0, r0 + 1, r0, r0 + 1};
c0 = {c0, c0, c0 + 1, c0 + 1};
wr = {1 - wr, wr, 1 - wr, wr};
wc = {1 - wc, 1 - wc, wc, wc};
end

function output = addGather(output, input, rows, cols, wr, wc, nRows, nCols)
for k = 1:4
    valid = rows{k} >= 1 & rows{k} <= nRows & cols{k} >= 1 & cols{k} <= nCols;
    if any(valid(:))
        dstIdx = find(valid);
        srcIdx = sub2ind([nRows, nCols], rows{k}(valid), cols{k}(valid));
        output(dstIdx) = output(dstIdx) + input(srcIdx) .* cast(wr{k}(valid) .* wc{k}(valid), 'like', input);
    end
end
end

function output = addScatter(output, input, rows, cols, wr, wc, nRows, nCols)
for k = 1:4
    valid = rows{k} >= 1 & rows{k} <= nRows & cols{k} >= 1 & cols{k} <= nCols;
    if any(valid(:))
        dstIdx = sub2ind([nRows, nCols], rows{k}(valid), cols{k}(valid));
        vals = input(valid) .* cast(wr{k}(valid) .* wc{k}(valid), 'like', input);
        contrib = accumarray(double(dstIdx), double(vals), [nRows*nCols, 1], @sum, 0);
        output = output + reshape(cast(contrib, 'like', input), nRows, nCols);
    end
end
end

function factors = attenuationFactors(options, volume, view, extra)
% Native type-6 per-view attenuation factors from the union image-domain map.
volume=double(volume); view=double(view);
[~,gridXY]=volumeGeometry(options,volume); [gridNz,padLo,~]=axialGeometry(options,volume);
if nargin>=4 && numel(extra)>=2
    gridXY=double(extra{1}(:).'); gridNz=double(extra{2});
end
if ~isfield(options,'vaimennus') || isempty(options.vaimennus)
    factors=ones([gridXY,gridNz],options.cType); return
end
nVol=double(options.nMultiVolumes)+1;
Nx=double(options.Nx(1:nVol)); Ny=double(options.Ny(1:nVol)); Nz=double(options.Nz(1:nVol));
dx=double(options.dx(1:nVol)); dy=double(options.dy(1:nVol)); dz=double(options.dz(1:nVol));
bx=double(options.bx(1:nVol)); by=double(options.by(1:nVol)); bz=double(options.bz(1:nVol));
lo=[min(bx),min(by),min(bz)]; hi=[max(bx+Nx.*dx),max(by+Ny.*dy),max(bz+Nz.*dz)];
map=options.vaimennus;
if iscell(map)
    if numel(map)~=1, error('type6NativeOperator:DynamicAttenuationUnsupported','Per-frame type-6 image attenuation is not supported by the MATLAB helper.'); end
    map=map{1};
end
if isfield(options,'multiResUnionDims') && ~isempty(options.multiResUnionDims)
    mapDims=double(options.multiResUnionDims(:).');
else
    mapDims=double([size(map,1),size(map,2),size(map,3)]);
end
if numel(mapDims)~=3 || any(mapDims<1) || numel(map)~=prod(mapDims)
    error('type6NativeOperator:AttenuationDimensions','The type-6 attenuation map must cover its union grid.');
end
map=reshape(single(map),mapDims(1),mapDims(2),mapDims(3));
mapSpacing=[dx(1),dy(1),dz(1)]; mapOrigin=lo;
% In the native public map the first voxel is the low corner of the volume union.
if any(abs((hi-lo)./mapSpacing-mapDims)>1.1) && isfield(options,'multiResPadLow')
    mapOrigin=[bx(1),by(1),bz(1)]-double(options.multiResPadLow(:).').*mapSpacing;
end
runion=0;
for jj=1:nVol
    runion=max(runion,sqrt(max(bx(jj)^2,(bx(jj)+Nx(jj)*dx(jj))^2)+max(by(jj)^2,(by(jj)+Ny(jj)*dy(jj))^2)));
end
gridNx=gridXY(1); gridNy=gridXY(2);
minNx=ceil(2*runion/dx(volume)); minNy=ceil(2*runion/dy(volume));
attNx=max(1,minNx+mod(minNx-gridNx,2)); srcNy=max(1,minNy+mod(minNy-gridNy,2));
srcNz=max(1,ceil((hi(3)-lo(3))/dz(volume)-1e-10));
sourceDims=[attNx,srcNy,srcNz]; sourceSpacing=[dx(volume),dy(volume),dz(volume)];
sourceOrigin=[-.5*attNx*sourceSpacing(1),-.5*srcNy*sourceSpacing(2),lo(3)];
source=resamplePhysical3D(map,mapOrigin,mapSpacing,sourceDims,sourceOrigin,sourceSpacing);
% Rotate onto the emission grid's axial-centred row lattice. This source offset
% is the native offsetsFPType6 panel shift expressed in physical [column,row].
swivel=double(options.swivelAngles(view));
panel=double(options.blurPlanes2F(view))*dx(1)/dx(volume);
if nargin>=4 && numel(extra)>=3, panel=double(extra{3}); end
offX=-sind(swivel)*panel; offY=cosd(swivel)*panel;
geom=struct('inputSize',[attNx,srcNy],'outputSize',[attNx,gridNy]);
rotated=nativeRotationGather(source,-swivel,offX,offY,geom);
cumulative=cumsum(rotated,1); previous=zeros(size(cumulative),'like',cumulative); previous(2:end,:,:)=cumulative(1:end-1,:,:);
mode=2; envMode=getenv('TYPE6_ATTCONV'); if ~isempty(envMode), parsed=str2double(envMode(1)); if isfinite(parsed)&&ismember(parsed,[0,1,2]), mode=parsed; end, end
if mode==0
    factorsGrid=exp(-dx(volume).*cumulative);
elseif mode==1
    factorsGrid=exp(-dx(volume).*(cumulative+previous).*0.5);
else
    interval=dx(volume).*(cumulative-previous); average=ones(size(interval),'like',interval);
    nonzero=abs(interval)>1e-6; average(nonzero)=-expm1(-interval(nonzero))./interval(nonzero);
    small=~nonzero; average(small)=1-interval(small).*0.5+interval(small).^2/6;
    factorsGrid=exp(-dx(volume).*previous).*average;
end
rowStart=floor((attNx-gridNx)/2)+1; rows=rowStart:(rowStart+gridNx-1);
factors=ones(gridNx,gridNy,gridNz,'single');
for iz=1:gridNz
    zPhysical=bz(volume)+(iz-padLo-0.5)*dz(volume);
    zMap=(zPhysical-sourceOrigin(3))/sourceSpacing(3)+0.5;
    plane=sampleVolumePlane(factorsGrid,zMap);
    if rows(1)>=1 && rows(end)<=attNx, factors(:,:,iz)=plane(rows,:); end
end
factors=cast(factors,options.cType);
end

function output = resamplePhysical3D(input,inputOrigin,inputSpacing,outputDims,outputOrigin,outputSpacing)
% Mean-preserving box resampling via a trilinearly sampled summed-area table.
[nx,ny,nz]=size(input);
ex= min(max((outputOrigin(1)+(0:outputDims(1))*outputSpacing(1)-inputOrigin(1))/inputSpacing(1),0),nx);
ey= min(max((outputOrigin(2)+(0:outputDims(2))*outputSpacing(2)-inputOrigin(2))/inputSpacing(2),0),ny);
ez= min(max((outputOrigin(3)+(0:outputDims(3))*outputSpacing(3)-inputOrigin(3))/inputSpacing(3),0),nz);
sat=zeros(nx+1,ny+1,nz+1,'double');
sat(2:end,2:end,2:end)=cumsum(cumsum(cumsum(double(input),1),2),3);
[X,Y,Z]=ndgrid(ex+1,ey+1,ez+1); corners=trilinearSample(sat,X,Y,Z);
boxSum=corners(2:end,2:end,2:end)-corners(1:end-1,2:end,2:end)-corners(2:end,1:end-1,2:end)-corners(2:end,2:end,1:end-1) ...
    +corners(1:end-1,1:end-1,2:end)+corners(1:end-1,2:end,1:end-1)+corners(2:end,1:end-1,1:end-1)-corners(1:end-1,1:end-1,1:end-1);
nominalVolume=prod(outputSpacing./inputSpacing);
output=cast(boxSum./max(nominalVolume,1e-6),'like',input);
end

function output = sampleVolumePlane(input, zIndex)
lo=floor(zIndex); f=zIndex-lo; output=zeros(size(input,1),size(input,2),'like',input);
if lo>=1 && lo<=size(input,3), output=output+input(:,:,lo).*(1-f); end
hi=lo+1; if hi>=1 && hi<=size(input,3), output=output+input(:,:,hi).*f; end
end

function output = trilinearSample(input, x, y, z)
[nx,ny,nz]=size(input); x0=floor(x); y0=floor(y); z0=floor(z); fx=x-x0; fy=y-y0; fz=z-z0;
output=zeros(size(x),'like',input);
for oz=0:1, for oy=0:1, for ox=0:1
    ix=x0+ox; iy=y0+oy; iz=z0+oz;
    valid=ix>=1 & ix<=nx & iy>=1 & iy<=ny & iz>=1 & iz<=nz;
    if any(valid(:))
        idx=sub2ind([nx,ny,nz],ix(valid),iy(valid),iz(valid));
        wx=fx; if ox==0, wx=1-fx; end
        wy=fy; if oy==0, wy=1-fy; end
        wz=fz; if oz==0, wz=1-fz; end
        output(valid)=output(valid)+input(idx).*cast(wx(valid).*wy(valid).*wz(valid),'like',input);
    end
end, end, end
end

function kernel = cdrfKernel(options, volume, view)
% Per-view Gaussian CDRF on the expanded rotation grid (Python reference).
volume=double(volume); view=double(view);
if isfield(options,'cdrfModel6') && double(options.cdrfModel6) ~= 0
    error('type6NativeOperator:UnsupportedCDRF','The native MATLAB helper currently supports only the Gaussian on-the-fly CDRF model.');
end
if isfield(options,'CDRFOnTheFly') && ~options.CDRFOnTheFly
    error('type6NativeOperator:UnsupportedCDRF','Manual type-6 PSF banks are not supported by the native MATLAB helper.');
end
[~,gridXY]=volumeGeometry(options,volume); [gridNz,~,~]=axialGeometry(options,volume);
gridNx=gridXY(1); gridNy=gridXY(2); gridNz=double(gridNz);
dx=double(options.dx(volume)); dy=double(options.dy(volume)); dz=double(options.dz(volume)); dx0=double(options.dx(1));
distFloat=double(options.blurPlanesF(view)); collimator=(distFloat-.5*double(options.Nx(1)))*dx0;
distance=(collimator+.5*gridNx*dx)/dx;
firstPlane=max(0,min(gridNx-1,ceil(distance-.5))); firstPlane=floor(firstPlane/16)*16;
coeff=double(options.CDRFCoeffs6(:)); Az=coeff(1); Bz=coeff(2); Axy=coeff(3); Bxy=coeff(4); intrinsic=coeff(5);
sigmaCoeff=1/(2*sqrt(2*log(2)));
planes=(0:gridNx-1).'; depth=(planes+.5-distance)*dx;
rgZ=max(Az.*depth+Bz,0); rgXY=max(Axy.*depth+Bxy,0);
sigmaZ=sqrt(rgZ.^2+intrinsic^2)./dz.*sigmaCoeff;
fwhmXY=sqrt(rgXY.^2+intrinsic^2)./dy;
if ~isfield(options,'type6RotationComp') || options.type6RotationComp
    comp=rotationVariance(options,volume,view,gridXY);
else
    comp=sigmaCoeff^2;
end
sigmaXY=sqrt(max((fwhmXY.*sigmaCoeff).^2-comp,1e-16*sigmaCoeff^2));
extent=0;
for zmm=[(firstPlane+.5-distance)*dx,(gridNx-.5-distance)*dx]
    rz=max(Az*zmm+Bz,0); rxy=max(Axy*zmm+Bxy,0);
    extent=max(extent,5.2565*max(sqrt(rz^2+intrinsic^2)/dz,sqrt(rxy^2+intrinsic^2)/dy)*sigmaCoeff);
end
radius=min(ceil(extent/8)*8,max(gridNy,gridNz));
taps=-radius:radius;
kernel=zeros(gridNx,2*radius+1,2*radius+1,options.cType);
for ix=1:gridNx
    if ix-1<firstPlane || depth(ix)<=0, continue, end
    sy=max(sigmaXY(ix),1e-8); sz=max(sigmaZ(ix),1e-8);
    trans=exp(-(taps.^2)/(2*sy^2)); axial=exp(-(taps.^2)/(2*sz^2));
    trans(abs(taps)>ceil(5.2565*sy))=0;
    axial(abs(taps)>ceil(5.2565*sz))=0;
    if sum(trans)>0, trans=trans/sum(trans); end
    if sum(axial)>0, axial=axial/sum(axial); end
    kernel(ix,:,:)=cast(trans(:)*axial(:).','like',kernel);
end
end

function variance=rotationVariance(options,volume,view,grid)
% Port of Python _type6_rotation_variance; exact tent-weighted variance.
nx=double(options.Nx(volume)); ny=double(options.Ny(volume)); outNx=grid(1); outNy=grid(2);
theta=deg2rad(double(options.swivelAngles(view))); c=cos(theta); s=sin(theta);
panel=double(options.blurPlanes2F(view))*double(options.dx(1))/double(options.dx(volume));
axisX=-(double(options.bx(volume))+.5*nx*double(options.dx(volume)))/double(options.dx(volume));
axisY=-(double(options.by(volume))+.5*ny*double(options.dy(volume)))/double(options.dy(volume));
offX=-s*panel+axisX; offY=c*panel+axisY; hh=abs(c)+abs(s);
csx=.5*nx-.5; csy=.5*ny-.5; cox=.5*outNx-.5; coy=.5*outNy-.5;
if abs(s)<1e-6 || abs(c)<1e-6
    samples=[floor(nx/2),floor(ny/2)];
elseif nx*ny<=1024
    [ii,jj]=ndgrid(0:nx-1,0:ny-1); samples=[ii(:),jj(:)];
else
    % Deterministic Halton sampling, matching the Python fallback.
    samples=zeros(1024,2);
    for k=1:1024, samples(k,:)=[min(floor(halton(k,2)*nx),nx-1),min(floor(halton(k,3)*ny),ny-1)]; end
end
num=0; den=0;
for k=1:size(samples,1)
    ux=samples(k,1)-csx; uy=samples(k,2)-csy;
    qx=c*(ux-offX)+s*(uy-offY); qy=-s*(ux-offX)+c*(uy-offY);
    jx0=ceil(qx+cox-hh); jx1=floor(qx+cox+hh); jy0=ceil(qy+coy-hh); jy1=floor(qy+coy+hh);
    for jy=jy0:jy1
        ay=jy-coy;
        for jx=jx0:jx1
            ax=jx-cox; sx=c*ax-s*ay+offX; sy=s*ax+c*ay+offY;
            weight=max(1-abs(sx-ux),0)*max(1-abs(sy-uy),0);
            num=num+weight*(ay-qy)^2; den=den+weight;
        end
    end
end
if den>0, variance=num/den; else, variance=1/6; end
if volume>1
    nDet=double(options.nRowsD); a=double(options.dy(1))/double(options.dy(volume)); b=.5*outNy+(.5-.5*nDet)*a-.5;
    num=0; den=0;
    for row=0:nDet-1
        t=a*row+b; low=floor(t);
        for off=0:1, d=t-(low+off); w=max(1-abs(d),0); num=num+w*d^2; den=den+w; end
    end
    if den>0, variance=variance+num/den; end
end
end

function value=halton(index,base)
value=0; factor=1/base;
while index>0, value=value+factor*mod(index,base); index=floor(index/base); factor=factor/base; end
end

function [gridNz, padLo, padHi] = axialGeometry(options, volume)
% Axial blur spill padding from computePadZType6 in precomp.h.
volume = double(volume);
Nz = double(options.Nz(volume)); dz = double(options.dz(volume));
padLo = 0; padHi = 0; padZ = 0;
needsPad = numel(options.Nx) > 1 || ...
    abs(double(options.bx(1)) + 0.5 * double(options.Nx(1)) * double(options.dx(1))) > 1e-3 * double(options.dx(1)) || ...
    abs(double(options.by(1)) + 0.5 * double(options.Ny(1)) * double(options.dy(1))) > 1e-3 * double(options.dy(1)) || ...
    abs(double(options.bz(1)) + 0.5 * double(options.Nz(1)) * double(options.dz(1))) > 1e-3 * double(options.dz(1));
if needsPad && isfield(options, 'CDRFCoeffs6') && numel(options.CDRFCoeffs6) >= 5 && ...
        isfield(options, 'blurPlanesF') && ~isempty(options.blurPlanesF)
    coeff = double(options.CDRFCoeffs6(:));
    minDist = min(double(options.blurPlanesF(:)));
    geom = type6NativeOperator('volume-geometry', options, volume);
    gridNx = geom.outputSize(1);
    collimator = (minDist - 0.5 * double(options.Nx(1))) * double(options.dx(1));
    dist = (collimator + 0.5 * gridNx * double(options.dx(volume))) / double(options.dx(volume));
    zMax = (gridNx - 0.5 - dist) * double(options.dx(volume));
    sigmaCoeff = 1 / (2 * sqrt(2 * log(2)));
    extentMax = 0;
    model = 0;
    if isfield(options, 'cdrfModel6'), model = double(options.cdrfModel6); end
    for z = [0, max(zMax, 0)]
        Rg = max(coeff(1) * z + coeff(2), 0);
        C = coeff(5);
        if model == 1
            extent = Rg / dz + 5.2565 * (C / dz * sigmaCoeff);
        else
            extent = 5.2565 * (sqrt(Rg^2 + C^2) / dz * sigmaCoeff);
        end
        extentMax = max(extentMax, extent);
    end
    radius = ceil(extentMax);
    padZ = ceil(radius / 8) * 8;
    zColLo = -0.5 * double(options.nColsD) * double(options.dz(1)) + 0.5 * double(options.dz(1));
    zColHi = 0.5 * double(options.nColsD) * double(options.dz(1)) - 0.5 * double(options.dz(1));
    posLo = (zColLo - double(options.bz(volume))) / dz - 0.5;
    posHi = (zColHi - double(options.bz(volume))) / dz - 0.5;
    needLo = -floor(posLo + 1e-4);
    needHi = ceil(posHi - 1e-4) - (Nz - 1);
    padLo = max(min(padZ, needLo), 0);
    padHi = max(min(padZ, needHi), 0);
end
gridNz = Nz + padLo + padHi;
end

function [geom, gridSize] = volumeGeometry(options, volume)
% Physical XY geometry of the smallest axis-centered square grid that contains
% one reconstruction volume, following rotGridType6 in precomp.h.
volume = double(volume);
Nx = double(options.Nx(volume)); Ny = double(options.Ny(volume));
if volume == 1 && numel(options.Nx) == 1 && ...
        abs(double(options.bx(1)) + 0.5*Nx*double(options.dx(1))) < 1e-3*double(options.dx(1)) && ...
        abs(double(options.by(1)) + 0.5*Ny*double(options.dy(1))) < 1e-3*double(options.dy(1))
    geom=struct('inputSize',[Nx,Ny],'outputSize',[Nx,Ny], ...
        'inputSpacing',[double(options.dx(1)),double(options.dy(1))], ...
        'outputSpacing',[double(options.dx(1)),double(options.dy(1))], ...
        'inputOffset',[0,0],'outputOffset',[0,0]);
    gridSize=[Nx,Ny]; return
end
dx = double(options.dx(volume)); dy = double(options.dy(volume));
bx = double(options.bx(volume)); by = double(options.by(volume));
x0 = bx; x1 = bx + Nx * dx; y0 = by; y1 = by + Ny * dy;
R = sqrt(max(x0.^2, x1.^2) + max(y0.^2, y1.^2));
Rcol = R;
if 0.5 * double(options.nRowsD) * double(options.dy(1)) > R
    shifts = abs(double(options.radiusPerProj(:)) .* sind(double(options.swivelAngles(:)) - double(options.angles(:)) + 180));
    if ~isempty(shifts), Rcol = R + max(shifts); end
end
nxGrid = fitGridDimension(R, dx, double(options.Nx(1)), volume == 1);
nyGrid = fitGridDimension(Rcol, dy, double(options.Ny(1)), volume == 1);
centerX = bx + 0.5 * Nx * dx;
centerY = by + 0.5 * Ny * dy;
geom = struct('inputSize', [Nx, Ny], 'outputSize', [nxGrid, nyGrid], ...
    'inputSpacing', [dx, dy], 'outputSpacing', [dx, dy], ...
    'inputOffset', [centerY, centerX], 'outputOffset', [0, 0]);
gridSize = [nxGrid, nyGrid];
end

function n = fitGridDimension(radius, spacing, mainSize, isMain)
n = ceil(2 * radius / spacing);
if isMain, parity = mod(round(mainSize), 2); else, parity = 0; end
if mod(n, 2) ~= parity, n = n + 1; end
n = max(n, 1);
end

function output = gridToDetector(input, detectorRows, trA, trB)
% Native forward resampling: sample the grid at t=trA*r+trB for each detector row r.
detectorRows=double(detectorRows); trA=double(trA); trB=double(trB);
if detectorRows<1 || trA<=0, error('type6NativeOperator:ResampleGeometry','Detector rows and scale must be positive.'); end
sz=size(input); nGrid=sz(1); rest=prod(sz(2:end)); input=reshape(input,nGrid,rest);
output=zeros(detectorRows,rest,'like',input);
for r=0:detectorRows-1
    t=trA*r+trB; lo=floor(t); f=t-lo;
    if lo>=0 && lo<nGrid, output(r+1,:)=output(r+1,:)+input(lo+1,:).*cast(1-f,'like',input); end
    hi=lo+1;
    if hi>=0 && hi<nGrid, output(r+1,:)=output(r+1,:)+input(hi+1,:).*cast(f,'like',input); end
end
output=reshape(output,[detectorRows,sz(2:end)]);
end

function output = detectorToGridAdjoint(input, gridRows, trA, trB)
% Exact transpose: scatter each detector row to floor/ceil(trA*r+trB).
gridRows=double(gridRows); trA=double(trA); trB=double(trB);
if gridRows<1 || trA<=0, error('type6NativeOperator:ResampleGeometry','Grid rows and scale must be positive.'); end
sz=size(input); nDet=sz(1); rest=prod(sz(2:end)); input=reshape(input,nDet,rest);
output=zeros(gridRows,rest,'like',input);
for r=0:nDet-1
    t=trA*r+trB; lo=floor(t); f=t-lo;
    if lo>=0 && lo<gridRows, output(lo+1,:)=output(lo+1,:)+input(r+1,:).*cast(1-f,'like',input); end
    hi=lo+1;
    if hi>=0 && hi<gridRows, output(hi+1,:)=output(hi+1,:)+input(r+1,:).*cast(f,'like',input); end
end
output=reshape(output,[gridRows,sz(2:end)]);
end

function output = convolutionAdjoint(input, kernel)
% Exact transpose of MATLAB conv2(input,kernel,'same'). The selected full-
% convolution crop starts at floor(size(kernel)/2)+1 on each axis.
if ndims(input) > 2 || ndims(kernel) > 2
    error('type6NativeOperator:ConvolutionDimensions', 'Convolution adjoint expects 2-D arrays.');
end
[nRows, nCols] = size(input);
[kRows, kCols] = size(kernel);
r0 = floor(kRows / 2) + 1;
c0 = floor(kCols / 2) + 1;
output = zeros(size(input), 'like', input);
for kr = 1:kRows
    dr = kr - r0;
    if dr >= 0
        srcRows = (1 + dr):nRows;
        dstRows = 1:(nRows - dr);
    else
        srcRows = 1:(nRows + dr);
        dstRows = (1 - dr):nRows;
    end
    if isempty(srcRows), continue, end
    for kc = 1:kCols
        dc = kc - c0;
        if dc >= 0
            srcCols = (1 + dc):nCols;
            dstCols = 1:(nCols - dc);
        else
            srcCols = 1:(nCols + dc);
            dstCols = (1 - dc):nCols;
        end
        if ~isempty(srcCols)
            output(dstRows, dstCols) = output(dstRows, dstCols) + ...
                input(srcRows, srcCols) .* cast(kernel(kr, kc), 'like', input);
        end
    end
end
end

function [nRows, nCols, nPlanes] = size3(input)
s = size(input);
nRows = s(1);
nCols = s(2);
if numel(s) < 3
    nPlanes = 1;
else
    nPlanes = prod(s(3:end));
end
input = reshape(input, nRows, nCols, nPlanes); %#ok<NASGU>
end
