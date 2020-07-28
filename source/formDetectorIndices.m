function L = formDetectorIndices(det_w_pseudo)
%FORMDETECTORINDICES Outputs the detector index for the specified detector
%numbers
% This is a helper function
L = zeros(sum(1:det_w_pseudo),2,'int32');
jh = int32(1);
for kk = int32(1) : (det_w_pseudo)
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repeat_elem((kk), det_w_pseudo-(kk-1)), ((kk):det_w_pseudo)'];
    elseif exist('OCTAVE_VERSION','builtin') == 5
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1)), ((kk):det_w_pseudo)'];
    else
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1))', ((kk):det_w_pseudo)'];
    end
    jh = jh + (det_w_pseudo) -kk + 1;
end
L(L(:,1) == 0,:) = [];
end