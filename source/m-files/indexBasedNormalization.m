function options = indexBasedNormalization(options, forSensitivity, isCorrectionFactor)
%INDEXBASEDNORMALIZATION Converts the normalization of index-based data
%
% Converts the normalization of index-based list-mode data into the form
% used by the reconstruction: one value for each event
% (options.normalization). The normalization can be given either as a
% sinogram (Ndist * Nang * NSinos elements, or Ndist * Nang * TotSinos
% elements of which the first NSinos planes are used) or as one value for
% each event. The sinogram value of an event is found with the lookup
% tables of sinogramPairLookup, i.e. the event is placed to the same
% sinogram element as the loaders would place it. Events that are outside
% the sinogram (radial FOV, ring difference) get the value zero.
%
% A sinogram bin is the sum of all the LORs (detector pair and ring pair)
% that are placed to it (span > 1, mashing), but an index-based event is a
% single LOR. A sinogram normalization is therefore divided by the number
% of LORs in its bin (or multiplied by it, if it is a correction factor)
% before it is used. A normalization with one value for each event is used
% as it is, so it has to be given for a single LOR.
%
% With the sinogram form and forSensitivity set, the sinogram and the lookup
% tables are also stored in options.normSens, options.sensXY, options.sensZ
% and options.sensBins, which are used to apply the normalization in the
% sensitivity image (all detector pairs, not only the measured events).
%
% If a pair normalization table was input (options.normSens without
% options.sensXY, see OMEGA_error_check), the normalization must already be
% given as one value for each event. It is then used as it is and the table
% is not touched (apart from its conversion to single precision).
%
% INPUTS:
%   options = Struct containing the reconstruction options. The
%   normalization has to be already inverted, if it was loaded from a file.
%   forSensitivity = If true, the normalization is also used in the
%   sensitivity image if it is computed (options.compute_sensitivity_image)
%   isCorrectionFactor = (optional) If true, a sinogram normalization is a
%   correction factor that multiplies the data (the normalization in the
%   precorrection), so it is multiplied by the number of LORs in the bin
%   instead of dividing by it. The default is false (the normalization is
%   a multiplier in the reconstruction).
%
% OUTPUTS:
%   options = The same struct with the normalization as single precision
%   values for each event
%
% See also sinogramPairLookup, loadCorrections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    isCorrectionFactor = false;
end
% The sinogram dimensions are collapsed for index-based data, use the saved ones
if isfield(options, 'indexSinogramDims')
    dims = double(options.indexSinogramDims);
else
    dims = double([options.Ndist, options.Nang, options.NSinos, options.TotSinos]);
end
Ndist = dims(1);
Nang = dims(2);
NSinos = dims(3);
TotSinos = dims(4);
nSino = Ndist * Nang * NSinos;
nTot = Ndist * Nang * TotSinos;
if iscell(options.trIndex)
    if numel(options.trIndex) > 1
        error('Normalization correction is not supported for dynamic index-based data!')
    end
    trIndex = options.trIndex{1};
    axIndex = options.axIndex{1};
else
    trIndex = options.trIndex;
    axIndex = options.axIndex;
end
nEvents = numel(trIndex) / 2;
normalization = options.normalization(:);
if isfield(options, 'normSens') && ~isempty(options.normSens) && ~(isfield(options, 'sensXY') && ~isempty(options.sensXY))
    % A pair table (options.normSens, see OMEGA_error_check) is used in the sensitivity image, the normalization has to be given
    % for the events and is used as it is: a sinogram is not converted, since that would replace the pair table
    if numel(normalization) ~= nEvents
        error(['With a pair normalization table (options.normSens) options.normalization has to contain one value for each event (' ...
            num2str(nEvents) ' elements), but it has ' num2str(numel(normalization)) ' elements! The table is used only in the sensitivity ' ...
            'image, a sinogram normalization is not converted to the events.'])
    end
    options.normSens = single(options.normSens(:));
    options.normalization = single(normalization);
    return
end
if numel(normalization) == nSino || numel(normalization) == nTot
    if options.nLayers > 1
        error(['Index-based reconstruction supports a sinogram normalization only with a single layer! Input the normalization ' ...
            'as one value for each event (' num2str(nEvents) ' elements) instead.'])
    end
    if options.det_w_pseudo ~= options.det_per_ring || (~isempty(options.pseudot) && any(options.pseudot(:) > 0))
        error(['Index-based reconstruction supports a sinogram normalization only without pseudo detectors and pseudo rings! ' ...
            'Input the normalization as one value for each event (' num2str(nEvents) ' elements) instead.'])
    end
    detPerRing = double(options.det_per_ring);
    rings = double(options.rings);
    if mod(floor(detPerRing / 2), Nang) ~= 0
        error(['The number of sinogram angles (' num2str(Nang) ') has to divide the number of detectors per ring divided by two (' ...
            num2str(floor(detPerRing / 2)) ') in order to use a sinogram normalization!'])
    end
    tr = double(reshape(trIndex, 2, []));
    ax = double(reshape(axIndex, 2, []));
    if ~isempty(tr) && max(tr(:)) >= detPerRing
        error(['options.trIndex contains detector indices (maximum ' num2str(max(tr(:))) ') that are beyond the number of detectors per ring (' num2str(detPerRing) ')!'])
    end
    if ~isempty(ax) && max(ax(:)) >= rings
        error(['options.axIndex contains ring indices (maximum ' num2str(max(ax(:))) ') that are beyond the number of rings (' num2str(rings) ')!'])
    end
    sinogram = single(normalization(1 : nSino));
    [xyLookup, zLookup] = sinogramPairLookup(detPerRing, rings, Ndist, Nang, options.ring_difference, options.span, ...
        options.segment_table, NSinos, options.ndist_side);
    % A sinogram bin is the sum of its LORs (span > 1, mashing) but each event is a single LOR, so the number of
    % LORs in the bin (detector pairs a < b in the radial/angular bin times ordered ring pairs in the plane) is used
    outside = intmax('uint32');
    nz = accumarray(double(zLookup(zLookup ~= outside)) + 1, 1, [NSinos, 1]);
    xyPairs = xyLookup(triu(true(detPerRing), 1));
    nt = accumarray(double(bitand(xyPairs(xyPairs ~= outside), uint32(2^31 - 1))) + 1, 1, [Ndist * Nang, 1]);
    clear xyPairs
    sinogram = reshape(sinogram, Ndist * Nang, NSinos);
    for k = 1 : NSinos
        count = single(nt * nz(k));
        reachable = count > 0;
        if isCorrectionFactor
            sinogram(reachable, k) = sinogram(reachable, k) .* count(reachable);
        else
            sinogram(reachable, k) = sinogram(reachable, k) ./ count(reachable);
        end
    end
    sinogram = sinogram(:);
    xy = xyLookup(tr(1, :).' + tr(2, :).' * detPerRing + 1);
    % The rings are exchanged for the events whose detectors are swapped in the sinogram
    flip = bitand(xy, uint32(2^31)) ~= 0;
    ring1 = ax(1, :).';
    ring2 = ax(2, :).';
    ring1(flip) = ax(2, flip).';
    ring2(flip) = ax(1, flip).';
    plane = zLookup(ring1 + ring2 * rings + 1);
    inside = xy ~= intmax('uint32') & plane ~= intmax('uint32');
    index = double(bitand(xy, uint32(2^31 - 1))) + double(plane) * (Ndist * Nang) + 1;
    eventNormalization = zeros(nEvents, 1, 'single');
    eventNormalization(inside) = sinogram(index(inside));
    options.normalization = eventNormalization;
    if options.compute_sensitivity_image && forSensitivity
        options.normSens = sinogram;
        options.sensXY = xyLookup;
        options.sensZ = zLookup;
        options.sensBins = uint32(Ndist * Nang);
    end
elseif numel(normalization) == nEvents
    % The normalization has been converted by an earlier call with the same options (the sinogram is stored in
    % options.normSens), it is used as it is
    if options.compute_sensitivity_image && forSensitivity && ~(isfield(options, 'normSens') && ~isempty(options.normSens))
        sizeText = [num2str(nSino) ' elements'];
        if nTot ~= nSino
            sizeText = [sizeText ', or ' num2str(nTot) ' elements of which the first ' num2str(nSino) ' are used'];
        end
        error(['The sensitivity image of index-based data needs the normalization as a sinogram (' sizeText ') in order to apply it ' ...
            'to all detector pairs, but one value for each event (' num2str(nEvents) ' elements) was input!'])
    end
    options.normalization = single(normalization);
else
    sizeText = ['Ndist * Nang * NSinos = ' num2str(nSino) ' elements'];
    if nTot ~= nSino
        sizeText = [sizeText ', or Ndist * Nang * TotSinos = ' num2str(nTot) ' elements of which the first ' num2str(nSino) ' are used'];
    end
    error(['The size of the normalization (' num2str(numel(normalization)) ' elements) does not match index-based data! It has to be either a sinogram ' ...
        '(' sizeText ') or one value for each event (' num2str(nEvents) ' elements).'])
end
end
