function varargout = volume3Dviewer(volume, varargin)
%VOLUME3DVIEWER Visualizes a stack of slices from any 3D (grayscale) image
%volume
%   Default setting uses the minimum and maximum values of the current
%   slice as the respective minimum and maximum values of the figure.
%   Default setting also uses transverse stack. The user can adjust the
%   minimum maximum values as well as the orientation.
%   Examples:
%       volume3Dviewer(volume)
%       volume3Dviewer(volume, [0 1])
%       volume3Dviewer(volume, 'fit')
%       volume3Dviewer(volume, [0 1], [1 0 0])
%       volume3Dviewer(volume, [], [0 1 0])
%       volume3Dviewer(volume, 0)
%   Inputs:
%       volume = The input 3D image volume
%       scale = The minimum and maximum values of the color range, i.e.
%       [minimum maximum]. If you input 'fit' instead, the minimum will be
%       minimum value in the volume and maximum the maximum value in the 
%       volume. If only scalar is input, the minimum is set as the scalar
%       and maximum as the image maximum value. Default uses the minimum
%       and maximum values for the current slice. 
%       orientation = The orientation. [0 0 1] is the default transverse
%       orientation. Other possibilities are [1 0 0] and [0 1 0].
%
% See also sliceViewer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2024 Nargiza Djurabekova, Ville-Veikko Wettenhovi
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
    % Create a figure
    varargout{1} = figure;
    if nargin >= 2 && ~isempty(varargin{1})
        if strcmp('fit',varargin{1})
            scale = [min(volume(:)) max(volume(:))];
        elseif numel(varargin{1}) == 1
            scale = [0 max(volume(:))];
        else
            scale = varargin{1};
        end
        if numel(scale) ~= 2
            error('Scale value needs to have a minimum and a maximum value!')
        end
    else
        scale = [];
    end
    orientation = [0 0];
    if nargin >= 3 && ~isempty(varargin{2})
        if varargin{2}(1) == 1
            orientation = [1 0];
        %     volume = flipud(permute(volume, [1 3 2]));
        elseif varargin{2}(2) == 1
            orientation = [0 1];
        %     volume = flipud(permute(volume, [3 2 1]));
        end
    end



    % Create a handle for the axes
    ax = axes('Parent', varargout{1}, 'Position', [0.1, 0.15, 0.8, 0.8]);



    % Create a slider
    slider = uicontrol('Parent', varargout{1}, 'Style', 'slider', 'Units', 'normalized', 'Position', [.05, .05, .9, .05], ...
        'Min', 1, 'Max', size(volume, 3), 'Value', 1, 'SliderStep', [1, 5]./(size(volume, 3)-1));



    % Add a callback function to the slider
    if exist('OCTAVE_VERSION','builtin') == 0 && ~verLessThan('matlab', '8.3')
        addlistener(slider, 'ContinuousValueChange', @(src, event) updateSlice(src, ax, volume, scale, orientation));
    else
        set(slider, 'Callback', @(src, event) updateSlice(src, ax, volume, scale, orientation));
    end



    % Add a callback for mouse hover over the image
    set(varargout{1}, 'WindowButtonMotionFcn', @(src, event) displayPixelValue(src, ax, volume, slider));

    % Initial display of the first slice
    updateSlice(slider, ax, volume, scale, orientation);
end



% Function to update the displayed slice
function updateSlice(src, ax, volume, scale, orientation)
    % Get the current slider value (z-slice number)
    slice_num = round(get(src, 'Value'));



    % Display the selected slice using imagesc
    if isempty(scale)
        if orientation(2) == 1
            imagesc(ax, squeeze(volume(slice_num, :, :)));
        elseif orientation(1) == 1
            imagesc(ax, squeeze(volume(:, slice_num, :)));
        else
            imagesc(ax, volume(:, :, slice_num));
        end
    else
        if orientation(2) == 1
            imagesc(ax, squeeze(volume(slice_num, :, :)), scale);
        elseif orientation(1) == 1
            imagesc(ax, squeeze(volume(:, slice_num, :)), scale);
        else
            imagesc(ax, volume(:, :, slice_num), scale);
        end
    end


    % Adjust the axis properties as needed (colormap, title, etc.)
    colormap(ax, gray); % Adjust the colormap as needed



    % Update the axis limits if necessary
    axis(ax, 'off', 'equal'); % Turn off axis labels and ensure equal aspect ratio

    title(ax, ['Slice ' num2str(slice_num) '/' num2str(size(volume,3))]);


    % Refresh the figure
    drawnow;
end



% Function to display pixel value on mouse hover
function displayPixelValue(~, ax, volume, slider)
    cp = get(ax, 'CurrentPoint');
    x = round(cp(1, 1));
    y = round(cp(1, 2));
    slice_num = round(get(slider, 'Value'));
    img = volume(:, :, slice_num);
    img_v = NaN;



    if (x > 0 && y > 0 && x <= size(img, 2) && y <= size(img, 1))
        img_v = img(y, x);
    end
    if (~isnan(img_v))
        title(ax, sprintf("Slice %i/%i : (%i, %i) = %i", slice_num, size(volume,3), x, y, img_v));
    else
        title(ax, sprintf("Slice %i/%i", slice_num, size(volume,3)));
    end
end
