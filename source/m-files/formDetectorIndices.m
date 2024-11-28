function L = formDetectorIndices(det_w_pseudo, varargin)
%FORMDETECTORINDICES Outputs the detector index for the specified detector
%numbers
% This is a helper function
% if ~isempty(varargin) && ~isempty(varargin{1})
%     nLayers = varargin{1};
% else
%     nLayers = 1;
% end
% if ~isempty(varargin) && ~isempty(varargin{2})
%     crystN = varargin{2};
% else
%     crystN = 0;
% end
L = zeros(sum(1:det_w_pseudo),2,'int32');
jh = int32(1);
for kk = int32(1) : (det_w_pseudo)
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repeat_elem((kk), det_w_pseudo-(kk-1)), ((kk):det_w_pseudo)'];
    elseif exist('OCTAVE_VERSION','builtin') == 5 && verLessThan('Octave','7')
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1)), ((kk):det_w_pseudo)'];
    else
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1))', ((kk):det_w_pseudo)'];
    end
    jh = jh + (det_w_pseudo) -kk + 1;
end
L(L(:,1) == 0,:) = [];
% if nLayers > 1
%     temp = (crystN:crystN:det_w_pseudo)';
%     L(ismember(L(:,1),temp),:) = [];
%     L(ismember(L(:,2),temp),:) = [];
% end
end