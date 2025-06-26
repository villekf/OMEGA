function [useCUDA] = checkCUDA(varargin)
%CHECKCUDA Checks whether the device is a NVIDIA device or not
msg = ArrayFire_OpenCL_device_info();
device = 0;
if ~isempty(varargin)
    device = varargin{1};
end
k1 = strfind(msg, ['[' num2str(device) ']']);
if isempty(k1)
    k1 = strfind(msg, ['-' num2str(device) '-']);
end
k3 = strfind(msg(k1:end), ':');
vendor = msg(k1:k1+k3(1)-2);
if ~isempty(strfind(vendor,'NVIDIA'))
    useCUDA = true;
else
    useCUDA = false;
end
end