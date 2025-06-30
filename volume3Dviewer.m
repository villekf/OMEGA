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
%       volume3Dviewer(volume, [], [], volume2)
%       volume3Dviewer(volume, 'fit', [], volume2)
%   Inputs:
%       volume = The input 3D image volume
%       scale = The minimum and maximum values of the color range, i.e.
%       [minimum maximum]. If you input 'fit' instead, the minimum will be
%       minimum value in the volume and maximum the maximum value in the
%       volume. If only scalar is input, the minimum is set as the scalar
%       and maximum as the image maximum value. Default uses the minimum
%       and maximum values for the current slice.
%       orientation = The orientation. [0 0 1] is the default transverse
%       orientation used with empty input. Other possibilities include
%       [1 0 0] and [0 1 0].
%       volume2 = An optional CT volume for PET/SPECT visualisation. Has to
%       match the first volume in dimensions. Only for MATLAB as Octave has
%       no support for transparency with imagesc.
% See also sliceViewer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2024-2025 Nargiza Djurabekova, Ville-Veikko Wettenhovi,
% Niilo Saarlemo
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
    %% Check environment
    USING_MATLAB = 0;
    if exist('OCTAVE_VERSION','builtin') == 0 && ~verLessThan('matlab', '8.3')
        USING_MATLAB = 1;
    end

    %% Figure creation
    varargout{1} = figure; % Create figure

    % Create a handle for the axes
    ax = axes('Parent', varargout{1}, 'Position', [0.1, 0.15, 0.8, 0.8]);

    %% Variable arguments
    scale = []; % Scale
    if nargin >= 2 && ~isempty(varargin{1})
        if strcmp('fit',varargin{1})
            scale = [min(volume(:)) max(volume(:))];
        elseif isscalar(varargin{1})
            scale = [varargin{1} max(volume(:))];
        else
            scale = varargin{1};
        end
        if numel(scale) ~= 2
            error('Scale value needs to have a minimum and a maximum value!')
        end
    end

    if nargin >= 3 && ~isempty(varargin{2})
        if varargin{2}(1) == 1
            volume = flipud(permute(volume, [1 3 2]));
        elseif varargin{2}(2) == 1
            volume = flipud(permute(volume, [3 2 1]));
        end
    end

    volume2 = []; % Volume 2 (CT image in PET/SPECT)
    if nargin >= 4 && ~isempty(varargin{3})
        if size(varargin{3}) == size(volume)
            volume2 = varargin{3};
        else
            error('The size of input volumes do not match')
        end
    end

    useColorbar = false;
    if nargin >= 5 && ~isempty(varargin{4})
        useColorbar = varargin{4};
    end

    %% Volume 1 UI controls
    % Create a slider for the slice
    sliderSlice = uicontrol('Parent', varargout{1}, 'Style', 'slider', 'Units', 'normalized', 'Position', [.05, .05, .9, .05], ...
        'Min', 1, 'Max', size(volume, 3), 'Value', 1, 'SliderStep', [1, 5]./(size(volume, 3)-1));

    % Add a callback function to the slider
    if USING_MATLAB
        addlistener(sliderSlice, 'ContinuousValueChange', @(src, event) updateSlice());
    else
        set(sliderSlice, 'Callback', @(src, event) updateSlice());
        disp("test")
    end

    % Add a callback for mouse hover over the image
    set(varargout{1}, 'WindowButtonMotionFcn', @(src, event) displayPixelValue(src, ax, volume, sliderSlice));

    %% Volume 2 scale, axes and UI controls
    if ~isempty(volume2) && USING_MATLAB % Octave does not support alpha
        ax2 = copyobj(ax, varargout{1}); % Create axes for other volume
        scale2 = [];
        if ~isempty(varargin{1}) && strcmp('fit',varargin{1})
            scale2 = [min(volume2(:)) max(volume2(:))];
        end

        sliderWeight = uicontrol('Parent', varargout{1}, 'Style', 'slider', 'Units', 'normalized', 'Position', [.05, .15, .05, .8], ...
        'Min', 0, 'Max', 1, 'Value', 0.5, 'SliderStep', [.05, .25]);
        if USING_MATLAB
            addlistener(sliderWeight, 'ContinuousValueChange', @(src, event) updateSlice());
        else
            set(sliderWeight, 'Callback', @(src, event) updateSlice());
        end
    end

    %% Initial display of the first slice
    updateSlice();

    %% Functions
    % Function to draw slice using imagesc
    function drawSlice(ax, volume, scale, slice_num)
        if isempty(scale) % Display the selected slice using imagesc
            imagesc(ax, volume(:, :, slice_num));
        else
            imagesc(ax, volume(:, :, slice_num), scale);
        end
    end

    % Function to update the displayed slice
    function updateSlice()
        % Get the current slider value (displayed z-slice)
        slice_num = round(get(sliderSlice, 'Value'));

        drawSlice(ax, volume, scale, slice_num) % Draw volume 1

        if ~isempty(volume2) && USING_MATLAB % Octave does not support alpha
            drawSlice(ax2, volume2, scale2, slice_num) % Draw volume 2

            colormap(ax, 'hot')
            colormap(ax2, 'gray')
            w = get(sliderWeight, 'Value');
            alpha(ax, w)
            alpha(ax2, 1-w)
            axis(ax2, "off", "equal")
            ax2.UserData = linkprop([ax,ax2],...
                {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
                'ydir','xdir','xlim','ylim'});
            if useColorbar
                colorbar(ax2, 'Position', [.9, .15, .02, .8])
            end
        else
            colormap(ax, 'gray');
        end
        if (USING_MATLAB && useColorbar); colorbar(ax, 'Position', [.85, .15, .02, .8]); end
        axis(ax, "off", "equal")
        title(ax, ['Slice ' num2str(slice_num) '/' num2str(size(volume,3))]);

        drawnow; % Refresh the figure
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
end
