function install_mex(varargin)
%install_mex Build the necessary MEX-files
% This file can be used to build all the necessary MEX-files. ROOT and/or
% OpenCL MEX-files are optional and require additional dependencies. ROOT
% support requires that ROOT is installed and OpenCL MEX-files require
% OpenCL headers and library and (for some features) ArrayFire library.
% OpenC headers and libraries as well as ArrayFire path are automatically
% checked from standard folders, but they can also be inputted manually.
% Second input is the OpenCL include folder, third OpenCL library folder
% (x64) and fourth is the ArrayFire folder (one containing include and
% lib-folders, e.g. with version 3 on Windows it is by default C:\Program
% Files\ArrayFire\v3.  For ROOT support the user should either ensure that
% the main ROOT folder is either on path or specify it in the fifth input.
% The installation of all features will be done automatically if the
% relevant headers and libraries are found. If not, a warning is displayed.
%
% Example:
%   install_mex
%   install_mex(1)
%   install_mex(1, OPENCL_INCLUDE_PATH, OPENCL_LIB_PATH, AF_PATH, ROOT_PATH, USE_CUDA, CUDA_PATH)
%   install_mex(1, '/path/to/cl.h', '/path/to/libOpenCL.so', '/path/to/Arrayfire', '/path/to/ROOT', true, '/path/to/cuda')
% INPUTS:
%   Input 1: If the value is set to 1, then the (possible) compiler errors
%   are displayed. Otherwise only warnings are displayed if the
%   compilations were unsuccessful. This input is optional. On Octave
%   compilation errors are shown regardless of the choice.
%   OPENCL_INCLUDE_PATH: The path to the OpenCL SDK include folder. Default
%   folder on Windows platform is taken from the environment variables. On
%   Linux, checks are made for Intel SDK, CUDA toolkit and AMD drivers.
%   When using Octave, installation on non-standard location with spaces on
%   path can cause problems.
%   OPENCL_LIB_PATH: The path to the OpenCL SDK library folder (folder
%   containing opencl.lib). Default folder on Windows platform is taken
%   from the environment variables. On Linux, checks are made for Intel
%   SDK, CUDA toolkit and AMD drivers. When using Octave, installation on
%   non-standard location with spaces on path can cause problems.
%   AF_PATH: The path to the ArrayFire installation folder. Default folder
%   on Windows platforms is taken from the environment variables. On Linux
%   /opt/arrayfire/ is used. When using Octave, installation on
%   non-standard location with spaces on path can cause problems.
%   ROOT_PATH: The path to the ROOT installation folder. Default folder on
%   Windows platforms is C:\Program Files\root, on others /opt/root
%   USE_CUDA: Whether the CUDA code is compiled. Default is TRUE.
%   CUDA_PATH: The path to the CUDA installation folder. Default folder on
%   Windows platforms is taken from the environment value (CUDA_PATH),
%   elsewhere it is /usr/local/cuda/. Only needed if USE_CUDA = true.
%
% OUTPUTS:
%   Various MEX-files needed by OMEGA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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



numvarargs = length(varargin);
if numvarargs > 8
    error('Requires at most eight input parameters')
end



if nargin > 0
    verbose = varargin{1};
    if nargin > 1 && ~isempty(varargin{2})
        opencl_include_path = varargin{2};
    else
        opencl_include_path = '';
    end
    if nargin > 2 && ~isempty(varargin{3})
        opencl_lib_path = varargin{3};
    else
        opencl_lib_path = '';
    end
    if nargin > 3 && ~isempty(varargin{4})
        af_path = varargin{4};
    else
        af_path = '';
    end
    if nargin > 4 && ~isempty(varargin{5})
        root_path = varargin{5};
    else
        root_path = '';
    end
    if nargin > 5 && ~isempty(varargin{6})
        use_CUDA = varargin{6};
    else
        use_CUDA = true;
    end
    if nargin > 6 && ~isempty(varargin{7})
        cuda_path = varargin{7};
    else
        cuda_path = '';
    end
    if nargin > 7 && ~isempty(varargin{8})
        implementation = varargin{8};
    else
        implementation = 0;
    end
else
    verbose = false;
    opencl_include_path = '';
    opencl_lib_path = '';
    af_path = '';
    root_path = '';
    use_CUDA = true;
    cuda_path = '';
    implementation = 0;
end

af_path_include = [af_path '/include'];

if ispc
    [~,homedir] = system('echo %userprofile%');
    if exist('OCTAVE_VERSION','builtin') == 5
        proggis1_orig = 'Program Files';
        proggis2_orig = 'Program Files (x86)';
        proggis1 = 'PROGRA~1';
        proggis2 = 'PROGRA~2';

        nvidia = 'NVIDIA~2';
        nvidia_orig = 'NVIDIA GPU Computing Toolkit';
        amd = 'AMDAPP~1';
        amd_orig = 'AMD APP SDK';
        if isempty(root_path)
            root_path = 'C:/PROGRA~2/root/';
            for yy = 6:12
                if exist(root_path,'dir') == 7
                    break;
                end
                for xx = 0:40
                    for zz = 0:40
                        if xx < 10
                            xc = ['0' num2str(xx)];
                        else
                            xc = num2str(xx);
                        end
                        if zz < 10
                            zc = ['0' num2str(zz)];
                        else
                            zc = num2str(zz);
                        end
                        dir = ['C:/root_v' num2str(yy) '.' xc '.' zc];
                        if exist(dir,'dir') == 7
                            root_path = dir;
                            break;
                        end
                        dir2 = ['C:/PROGRA~2/root_v' num2str(yy) '.' xc '.' zc];
                        if exist(dir2,'dir') == 7
                            root_path = dir2;
                            break;
                        end
                    end
                end
            end
%             root_path = 'C:/PROGRA~2/root/';
        end
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('INTELOCLSDKROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('CUDA_PATH')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('CUDA_PATH');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_include_path = strrep(opencl_include_path, nvidia_orig, nvidia);
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('OCL_ROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('OCL_ROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('AMDAPPSDKROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('AMDAPPSDKROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_include_path = strrep(opencl_include_path, amd_orig, amd);
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif exist([homedir '/Documents/OCL_SDK_Light'], 'dir') == 7
            opencl_include_path = [homedir '/Documents/OCL_SDK_Light/include'];
            opencl_lib_path = [homedir '/Documents/OCL_SDK_Light/lib/x86_64'];
        elseif exist([homedir '/Downloads/OCL_SDK_Light'], 'dir') == 7
            opencl_include_path = [homedir '/Downloads/OCL_SDK_Light/include'];
            opencl_lib_path = [homedir '/Downloads/OCL_SDK_Light/lib/x86_64'];
        elseif isempty(opencl_lib_path)
            opencl_include_path = '';
            opencl_lib_path = '';
            warning('No OpenCL SDK detected. If one is installed, insert the paths manually by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'')')
        end
        if isempty(af_path)
            af_path = getenv('AF_PATH');
            af_path = strrep(af_path, proggis2_orig, proggis2);
            af_path = strrep(af_path, proggis1_orig, proggis1);
            af_path_include = [af_path '/include'];
        end
    else
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('INTELOCLSDKROOT');
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('CUDA_PATH')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('CUDA_PATH');
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('OCL_ROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('OCL_ROOT');
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('AMDAPPSDKROOT')) && isempty(opencl_lib_path)
            opencl_include_path = getenv('AMDAPPSDKROOT');
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif exist([homedir '/Documents/OCL_SDK_Light'], 'dir') == 7
            opencl_include_path = [homedir '/Documents/OCL_SDK_Light/include'];
            opencl_lib_path = [homedir '/Documents/OCL_SDK_Light/lib/x86_64'];
        elseif exist([homedir '/Downloads/OCL_SDK_Light'], 'dir') == 7
            opencl_include_path = [homedir '/Downloads/OCL_SDK_Light/include'];
            opencl_lib_path = [homedir '/Downloads/OCL_SDK_Light/lib/x86_64'];
        elseif isempty(opencl_lib_path)
            opencl_include_path = '';
            opencl_lib_path = '';
        end
        if isempty(af_path)
            af_path = getenv('AF_PATH');
        end
        if isempty(root_path)
            [res, ldir] = system('root-config --libdir');
            if res == 0
                root_path = ldir(2:end-4);
            else
                root_path = 'C:/Program Files/root/';
                for yy = 6:12
                    if exist(root_path,'dir') == 7
                        break;
                    end
                    for xx = 0:60
                        if exist(root_path,'dir') == 7
                            break;
                        end
                        for zz = 0:60
                            if xx < 10
                                xc = ['0' num2str(xx)];
                            else
                                xc = num2str(xx);
                            end
                            if zz < 10
                                zc = ['0' num2str(zz)];
                            else
                                zc = num2str(zz);
                            end
                            dir = ['C:/root_v' num2str(yy) '.' xc '.' zc];
                            if exist(dir,'dir') == 7
                                root_path = dir;
                                break;
                            end
                            dir2 = ['C:/Program Files/root_v' num2str(yy) '.' xc '.' zc];
                            if exist(dir2,'dir') == 7
                                root_path = dir2;
                                break;
                            end
                        end
                    end
                end
            end
%             root_path = 'C:/PROGRA~2/root/';
        end
%         if isempty(root_path)
%             root_path = 'C:/Program Files/root/';
%         end
    end
    if isempty(cuda_path)
        cuda_path = getenv('CUDA_PATH');
    end
elseif ~ismac % Linux
    if isempty(cuda_path)
        cuda_path = '/usr/local/cuda';
    end
    [~,ldirA] = system('echo $HOME');
    ldir = [ldirA(1:end-1) '/arrayfire/'];
    sdir = fileparts(which('install_mex.m'));
    sdir = [sdir(1:end-20) '/arrayfire/'];
    sdirC = [sdir(1:end-20) '/Arrayfire/'];
    sdirC2 = [sdir(1:end-20) '/ArrayFire/'];
    if exist('/opt/arrayfire/','dir') == 7 && isempty(af_path)
        af_path = '/opt/arrayfire';
        af_path_include = '/opt/arrayfire/include/';
    elseif exist('/usr/local/include/af/','dir') == 7 && isempty(af_path)
        af_path = '/usr/local';
        af_path_include = '/usr/local/include/af';
    elseif exist('/usr/local/arrayfire/','dir') == 7 && isempty(af_path)
        af_path = '/usr/local/arrayfire';
        af_path_include = '/usr/local/arrayfire/include/';
    elseif exist(ldir,'dir') == 7  && isempty(af_path)
        af_path = ldir;
        af_path_include = [ldir 'include/'];
    elseif exist(sdir,'dir') == 7  && isempty(af_path)
        af_path = sdir;
        af_path_include = [sdir 'include/'];
    elseif exist(sdirC,'dir') == 7  && isempty(af_path)
        af_path = sdirC;
        af_path_include = [sdirC 'include/'];
    elseif exist(sdirC2,'dir') == 7  && isempty(af_path)
        af_path = sdirC2;
        af_path_include = [sdirC2 'include/'];
    elseif isempty(af_path)
        basePath = [ldirA(1:end-1) '/ArrayFire-'];
        AFound = false;
        for kk = 3:10
            for ll = 0:20
                for ii = 0:10
                    if exist([basePath num2str(kk) '.' num2str(ll) '.' num2str(ii) '-Linux/'],'dir') == 7
                        af_path = [basePath num2str(kk) '.' num2str(ll) '.' num2str(ii) '-Linux/'];
                        af_path_include = [[basePath num2str(kk) '.' num2str(ll) '.' num2str(ii) '-Linux/'] 'include/'];
                        AFound = true;
                        break;
                    end
                end
            end
        end
        if AFound == 0
            warning('ArrayFire not found. Please specify AF_PATH, if you wish to use implementation 2.')
            af_path = '';
            af_path_include = '';
        end
    end
    if exist(cuda_path,'dir') == 0 && isempty(cuda_path)
        breikki = false;
        for kk = 25 : -1 : 7
            for ll = 5 : -1 : 0
                cuda_path = ['/usr/local/cuda-' num2str(kk) '.' num2str(ll)];
                if exist(cuda_path,'dir') == 7
                    breikki = true;
                    break;
                end
            end
            if breikki
                break;
            end
        end
    end
    if exist('/opt/intel/opencl/include/','dir') == 7 && isempty(opencl_lib_path)
        opencl_include_path = '/opt/intel/opencl/include/';
        opencl_lib_path = '/opt/intel/opencl/lib/x64/';
    elseif exist('/usr/local/cuda/targets/x86_64-linux/include/','dir') == 7 && isempty(opencl_lib_path)
        opencl_include_path = '/usr/local/cuda/targets/x86_64-linux/include/';
        opencl_lib_path = '/usr/local/cuda/lib64/';
    elseif isempty(opencl_lib_path)
        opencl_include_path = '/usr/local/include/';
        opencl_lib_path = '/usr/lib/x86_64-linux-gnu/';
    end
else % Apple silicon

end

if ispc && (isempty(opencl_include_path) || isempty(opencl_lib_path))
    warning('No OpenCL SDK detected. If one is installed, insert the paths manually by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'')')
end

folder = fileparts(which('install_mex.m'));
folder = [folder(1:end-7), 'cpp'];
folder = strrep(folder, '\','/');

if ~isempty(opencl_include_path)
    opencl_include_path = strrep(opencl_include_path, '\','/');
end
if ~isempty(opencl_lib_path)
    opencl_lib_path = strrep(opencl_lib_path, '\','/');
end
if ~isempty(af_path)
    af_path = strrep(af_path, '\','/');
end
if ~isempty(root_path)
    root_path = strrep(root_path, '\','/');
end
if ~isempty(root_path) && strcmp(root_path(end),'/')
    root_path = [root_path 'bin/root-config'];
else
    if isempty(root_path)
        root_path = [root_path '/opt/root/bin/root-config'];
    else
        root_path = [root_path '/bin/root-config'];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATLAB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist('OCTAVE_VERSION','builtin') == 0) && ~ismac
    cc = mex.getCompilerConfigurations('C++','Selected');
    if isempty(cc)
        error('No C++ compiler selected! Use mex -setup C++ to select a C++ compiler')
    end
    if ispc
        OMPPath = ['"' matlabroot '/bin/win64"'];
        if strcmp(cc.Manufacturer, 'GNU')
            OMPLib = '';
        else
            OMPLib = '-liomp5md';
        end
        LPLib = '';
        ldflags = 'LDFLAGS="$LDFLAGS -fopenmp"';
        compiler = '';
    elseif isunix
        OMPPath = [matlabroot '/sys/os/glnxa64'];
        OMPLib = '-liomp5';
        LPLib = '-lpthread';
        ldflags = 'LDFLAGS="\$LDFLAGS -z noexecstack"';
        compiler = '';
        if strcmp(cc.Manufacturer, 'GNU')
            gccV = str2double(cc.Version);
            if verLessThan('matlab', '9.4')
                if gccV >= 5 || isnan(gccV)
                    if exist('/usr/bin/g++-4.9','file') == 2
                        compiler = 'GCC=/usr/bin/g++-4.9';
%                     elseif exist('/usr/bin/g++-6','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-6';
%                     elseif exist('/usr/bin/g++-7/','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-7';
%                     elseif exist('/usr/bin/g++-8/','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-8';
                    else
                        compiler = '';
                    end
                else
                    compiler = '';
                end
            elseif verLessThan('matlab', '9.9')
                if gccV >= 7 || isnan(gccV)
                    if exist('/usr/bin/g++-6','file') == 2
                        compiler = 'GCC=/usr/bin/g++-6';
%                     elseif exist('/usr/bin/g++-7/','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-7';
%                     elseif exist('/usr/bin/g++-8/','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-8';
%                     elseif exist('/usr/bin/g++-9/','file') == 2
%                         compiler = 'GCC=/usr/bin/g++-9';
                    else
                        compiler = '';
                    end
                else
                    compiler = '';
                end
            elseif verLessThan('matlab', '9.12')
                if gccV >= 10 || isnan(gccV)
                    if exist('/usr/bin/g++-9','file') == 2
                        compiler = 'GCC=/usr/bin/g++-9';
                    elseif exist('/usr/bin/g++-8','file') == 2
                        compiler = 'GCC=/usr/bin/g++-8';
                    elseif exist('/usr/bin/g++-7','file') == 2
                        compiler = 'GCC=/usr/bin/g++-7';
                    elseif exist('/usr/bin/g++-10','file') == 2
                        compiler = 'GCC=/usr/bin/g++-10';
                    else
                        compiler = '';
                    end
                else
                    compiler = '';
                end
            elseif verLessThan('matlab', '9.13')
                if gccV >= 11 || isnan(gccV)
                    if exist('/usr/bin/g++-10','file') == 2
                        compiler = 'GCC=/usr/bin/g++-10';
                    elseif exist('/usr/bin/g++-9','file') == 2
                        compiler = 'GCC=/usr/bin/g++-9';
                    elseif exist('/usr/bin/g++-8','file') == 2
                        compiler = 'GCC=/usr/bin/g++-8';
                    elseif exist('/usr/bin/g++-7','file') == 2
                        compiler = 'GCC=/usr/bin/g++-7';
                    else
                        compiler = '';
                    end
                else
                    compiler = '';
                end
            else
                compiler = '';
            end
        else
            compiler = '';
        end
    end
    if ~verLessThan('matlab','9.4')
        complexFlag = '-R2018a';
        %         complexFlag = '-R2017b';
    else
        complexFlag = '-largeArrayDims';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementations 1 & 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OMPh = '';
    if strcmp(cc.Manufacturer, 'Microsoft')
        compflags = 'COMPFLAGS="$COMPFLAGS /openmp"';
        cxxflags = 'CXXFLAGS="$CXXFLAGS"';
    elseif strcmp(cc.Manufacturer, 'Intel')
        if ispc
            compflags = 'COMPFLAGS="$COMPFLAGS /Qopenmp"';
        else
            compflags = 'COMPFLAGS="$COMPFLAGS -qopenmp"';
        end
        cxxflags = 'CXXFLAGS="$CXXFLAGS"';
    else
        compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
        cxxflags = 'CXXFLAGS="$CXXFLAGS -fopenmp"';
    end
    if implementation == 0 || implementation == 1 || implementation == 4
        try
            mex(compiler, complexFlag, '-outdir', folder, ['-L' OMPPath], ['-I ' folder], OMPh, OMPLib, LPLib, compflags, cxxflags, '-DMATLAB',...
                ldflags, [folder '/projector_mex.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'])
            mex(compiler, complexFlag, '-outdir', folder, ['-L' OMPPath], ['-I ' folder], OMPh, OMPLib, LPLib, compflags, cxxflags, '-DMATLAB',...
                ldflags, [folder '/projector_mexSingle.cpp'])
            disp('Implementations 1 & 4 built with OpenMP (parallel) support')
        catch ME
            mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11"', '-DMATLAB', ...
                [folder '/projector_mex.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'])
            mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11"', '-DMATLAB',...
                ldflags, [folder '/projector_mexSingle.cpp'])
            if verbose
                warning('Implementations 1 & 4 built WITHOUT OpenMP (parallel) support. Compiler error: ')
                disp(ME.message);
            else
                warning('Implementations 1 & 4 built WITHOUT OpenMP (parallel) support, Use install_mex(1) to see compiler error')
            end
        end
        try
            mex(compiler, complexFlag, '-outdir', folder, compflags, cxxflags, ['-I ' folder], ['-L' OMPPath], OMPh, OMPLib, LPLib, ldflags, ...
                [folder '/NLM_func.cpp'])
            mex(compiler, complexFlag, '-outdir', folder, compflags, cxxflags, ['-I ' folder], ['-L' OMPPath], OMPh, OMPLib, LPLib, ldflags, ...
                [folder '/NLM_funcSingle.cpp'])
        catch ME
            mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], [folder '/NLM_func.cpp'])
            mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], [folder '/NLM_funcSingle.cpp'])
            if verbose
                warning('NLM support for implementations 1 and 4 built WITHOUT OpenMP (parallel) support. Compiler error: ')
                disp(ME.message);
            else
                warning('NLM support for implementations 1 and 4 built WITHOUT OpenMP (parallel) support. Use install_mex(1) to see compiler error.')
            end
        end
        try
            if verLessThan('matlab','9.4')
                mex(compiler, complexFlag, '-outdir', folder, compflags, cxxflags, ['-L' OMPPath], OMPh, OMPLib, LPLib, ['-I ' folder], ldflags, ...
                    [folder '/createSinogramASCII.cpp'])
            else
                mex(compiler, '-largeArrayDims', '-outdir', folder, compflags, cxxflags, ['-L' OMPPath], OMPh, OMPLib, LPLib, ['-I ' folder], ldflags, ...
                    [folder '/createSinogramASCIICPP.cpp'])
            end
        catch ME
            if verLessThan('matlab','9.4')
                mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], [folder '/createSinogramASCII.cpp'])
            else
                mex(compiler, '-largeArrayDims', '-outdir', folder, ['-I ' folder], [folder '/createSinogramASCIICPP.cpp'])
            end
            if verbose
                warning('ASCII sinogram creation built WITHOUT OpenMP (parallel) support. Compiler error: ')
                disp(ME.message);
            else
                warning('ASCII sinogram creation built WITHOUT OpenMP (parallel) support. Use install_mex(1) to see compiler error.')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMF support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % try
        %     mex(compiler, '-largeArrayDims', '-outdir', folder, [folder '/gate_lmf_matlab.cpp'])
        %     disp('LMF support enabled')
        % catch ME
        %     if verbose
        %         warning('LMF support not enabled. Compiler error: ')
        %         disp(ME.message);
        %     else
        %         warning('LMF support not enabled. Use install_mex(1) to see compiler error.')
        %     end
        % end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inveon support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            mex(compiler, '-largeArrayDims', '-outdir', folder, ldflags, ['-I ' folder], [folder '/inveon_list2matlab.cpp'])
            disp('Inveon support enabled')
        catch ME
            if verbose
                warning('Inveon support not enabled. Compiler error: ')
                disp(ME.message);
            else
                warning('Inveon support not enabled. Use install_mex(1) to see compiler error.')
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROOT support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ispc
            if ~strcmp(cc.Manufacturer, 'Microsoft')
                warning('ROOT support is only supported when compiling with Visual Studio!')
            end
            try
                compflags = 'COMPFLAGS="$COMPFLAGS /std:c++17 -MD -EHsc -GR"';
                r_path = root_path;
                if strcmp(r_path,'root-config')
                    r_path = 'C:\root\';
                else
                    r_path = root_path(1:end-15);
                end
                mex(compiler, '-largeArrayDims', '-outdir', folder, compflags, '-DMATLABC', '-lCore', '-lTree', '-lRIO', ['-L"' r_path '/lib"'],...
                    ['-I"' r_path '/include"'], [folder '/GATE_root_matlab_C.cpp'])
                if verLessThan('matlab','9.6')
                else
                    mex(compiler, '-largeArrayDims', '-outdir', folder, compflags, '-DMATLABCPP', '-lCore', '-lTree', '-lRIO', ['-L"' r_path '/lib"'],...
                        ['-I"' r_path '/include"'], [folder '/GATE_root_matlab.cpp'])
                    % mex(compiler, '-largeArrayDims', '-outdir', folder, compflags, '-lCore', '-lTree', ['-L"' r_path '/lib"'],...
                    %     ['-I"' r_path '/include"'], [folder '/GATE_root_matlab_uint16.cpp'])
                    mex(compiler, '-largeArrayDims', '-outdir', folder, compflags, '-DMATLABC', '-lCore', '-lTree', '-lRIO', ['-L"' r_path '/lib"'],...
                        ['-I"' r_path '/include"'], [folder '/GATE_root_matlab_C.cpp'])
                end
                disp('ROOT support enabled')
            catch ME
                if verbose
                    warning('Unable to build ROOT support! Compiler error: ')
                    disp(ME.message);
                else
                    warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise, make sure that ROOT was compiled with the same compiler you are using to compile this' ...
                        ' mex-file and provide the path to ROOT install folder with install_mex(0, [], [], [], ''C:/path/to/ROOT''). Compiler error is shown with install_mex(1)']);
                end
            end
        else
            try
                if verLessThan('matlab','9.6')
                    mex(compiler, '-largeArrayDims', '-outdir', folder, ['CXXFLAGS="$CXXFLAGS $(' root_path ' --cflags)"'], '-DMATLABC', '-lCore', '-lTree', '-lRIO', '-ldl', ['LDFLAGS="$LDFLAGS $(' root_path ' --libs)"'], ...
                        [folder '/GATE_root_matlab_C.cpp'])
                    warning('Importing ROOT files will cause MATLAB to crash if you are not using R2019a or newer. Use MATLAB with `matlab -nojvm` to circumvent this.')
                else
                    mex(compiler, '-largeArrayDims', '-outdir', folder, ['CXXFLAGS="$CXXFLAGS $(' root_path ' --cflags)"'], '-DMATLABCPP', '-lCore', '-lTree', '-lRIO', '-ldl', '-lpthread', ...
                        ['LDFLAGS="$LDFLAGS $(' root_path ' --libs) -z noexecstack"'], ['-L' matlabroot '/sys/os/glnxa64'], ...
                        [folder '/GATE_root_matlab.cpp'])
                    % mex(compiler, '-largeArrayDims', '-outdir', folder, ['CXXFLAGS="$CXXFLAGS $(' root_path ' --cflags)"'], '-lCore', '-lTree', '-ldl', '-lpthread', ...
                    %     ['LDFLAGS="$LDFLAGS $(' root_path ' --libs)"'], ['-L' matlabroot '/sys/os/glnxa64'], ...
                    %     [folder '/GATE_root_matlab_uint16.cpp'])
                    mex(compiler, '-largeArrayDims', '-outdir', folder, ['CXXFLAGS="$CXXFLAGS $(' root_path ' --cflags)"'], '-DMATLABC', '-lCore', '-lTree', '-lRIO', '-ldl', ['LDFLAGS="$LDFLAGS $(' root_path ' --libs) -z noexecstack"'], ...
                        [folder '/GATE_root_matlab_C.cpp'])
                end
                disp('ROOT support enabled')
            catch ME
                if verbose
                    warning(['Unable to build ROOT support! Make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and that root-config '...
                        'is working properly, i.e. you have sourced the ROOT binaries with thisroot.*. You can manually input the ROOT installation folder with ' ...
                        'install_mex(0, [], [], [], ''/path/to/ROOT''). Compiler error: '])
                    disp(ME.message);
                else
                    warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
                        'the same compiler you are using to compile this mex-file and that root-config is working properly, i.e. you have sourced the ROOT binaries with thisroot.*. '...
                        'You can manually input the ROOT installation folder with install_mex(0, [], [], [], ''/path/to/ROOT''). Compiler error is shown with install_mex(1)']);
                end
            end
        end
    end

    % If you get any compilation errors then change the paths if
    % necessary. The examples below are for CUDA, AMD APP SDK (untested)
    % OCL SDK, and Intel SDK.

    %%%%%%%%%%%%%%%%%%%%%%%%% Implementations 2, 3 & 5 %%%%%%%%%%%%%%%%%%%%%%%%
    if (ispc)
        %%% Windows %%%
        if isempty(opencl_include_path)
            warning('No OpenCL SDK found. Implementations 2 and 3 will not be built. Use install_mex(1, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'') to set OpenCL include and library paths.')
        else
            if isempty(af_path)
                warning('ArrayFire not found on path. Implementation 2 will not be built. Use install_mex(1, [], [], ''C:/PATH/TO/ARRAYFIRE/'') to set ArrayFire path.')
            else
                if implementation == 0 || implementation == 2
                    if use_CUDA
                        if strcmp(cc.Manufacturer, 'Microsoft')
                        elseif strcmp(cc.Manufacturer, 'Intel')
                        else
                            compflags = 'COMPFLAGS="$COMPFLAGS -std=c++17"';
                            cxxflags = 'CXXFLAGS="$CXXFLAGS -Wp"';
                        end
                        try
                            disp('Attemping to build CUDA code.')
                            mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree', compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                                '-lnvrtc', ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                            mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree_uint16', compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE', ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                                '-lnvrtc', ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                            mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree_uint8', compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE2', ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                                '-lnvrtc', ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                            disp('CUDA support enabled')
                        catch
                            if strcmp(cc.Manufacturer, 'GNU')
                                warning('CUDA support is not available with MinGW. Please use Visual Studio instead.')
                            else
                                warning('CUDA support not enabled')
                            end
                        end
                    end
                end
                if strcmp(cc.Manufacturer, 'Microsoft')
                    compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                elseif strcmp(cc.Manufacturer, 'Intel')
                    %                     if ispc
                    compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    %                     else
                    %                         compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    %                     end
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                else
                    compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
                    cxxflags = 'CXXFLAGS="$CXXFLAGS -DOPENCL -Wno-ignored-attributes"';
                end
                if implementation == 0 || implementation == 2
                    try
                        %%%%%%%%%%%%%%%%%%%%%% Implementation 2 %%%%%%%%%%%%%%%%%%%%%%%
                        mex(compiler, complexFlag, '-outdir', folder, ['-I ' folder], ['-I"' opencl_include_path '"'], compflags, cxxflags, '-DMATLAB', '-DAF', '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib"'],...
                            ['-L"' opencl_lib_path '"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                        mex(compiler, complexFlag, '-outdir', folder, '-output', 'OpenCL_matrixfree_uint16', ['-I ' folder], ['-I"' opencl_include_path '"'], compflags, cxxflags, '-DMATLAB', '-DAF', '-DMTYPE', '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib"'],...
                            ['-L"' opencl_lib_path '"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                        mex(compiler, complexFlag, '-outdir', folder, '-output', 'OpenCL_matrixfree_uint8', ['-I ' folder], ['-I"' opencl_include_path '"'], compflags, cxxflags, '-DMATLAB', '-DAF', '-DMTYPE2', '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib"'],...
                            ['-L"' opencl_lib_path '"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                        mex(compiler, '-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib"'],['-L"' opencl_lib_path '"'], ...
                            ['-I ' folder], ['-I"' opencl_include_path '"'], ['-I"' af_path '/include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])

                        disp('Implementation 2 built')
                    catch ME
                        if verbose
                            warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                                'If OpenCL SDK and\or ArrayFire has been installed in a non-standard path, they can be added manually '...
                                'by using install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'', ''C:/PATH/TO/ARRAYFIRE/V3''). Compiler error:']);
                            disp(ME.message)
                        else
                            if strcmp(cc.Manufacturer, 'GNU')
                                warning(['Implementation 2 support with MinGW compiler requires ArrayFire to be compiled from source with MinGW. Compiler error is shown with install_mex(1).'...
                                    'If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. If you compiled ArrayFire with MinGW, add the path manually by using '...
                                    'install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'', ''/PATH/TO/ARRAYFIRE/V3'')'])
                            else
                                warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                                    'Compiler error is shown with install_mex(1). If OpenCL SDK and\or ArrayFire has been installed in a non-standard path, they can be added manually '...
                                    'by using install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'', ''/PATH/TO/ARRAYFIRE/V3'')']);
                            end
                        end
                    end
                end
                if implementation == 0 || implementation == 2
                    if strcmp(cc.Manufacturer, 'Microsoft')
                        compflags = 'COMPFLAGS="$COMPFLAGS /openmp -DCPU"';
                        cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                    elseif strcmp(cc.Manufacturer, 'Intel')
                        if ispc
                            compflags = 'COMPFLAGS="$COMPFLAGS /Qopenmp -DCPU"';
                        else
                            compflags = 'COMPFLAGS="$COMPFLAGS -qopenmp -DCPU"';
                        end
                        cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                    else
                        compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
                        cxxflags = 'CXXFLAGS="$CXXFLAGS -fopenmp -DCPU -Wno-ignored-attributes"';
                    end
                    try
                        %%%%%%%%%%%%%%%%%%%%%% Implementation 2 %%%%%%%%%%%%%%%%%%%%%%%
                        mex(compiler, complexFlag, '-outdir', folder, '-output', 'CPU_matrixfree', ['-L' OMPPath], OMPh, OMPLib, LPLib, ['-I ' folder], compflags, cxxflags, '-DMATLAB', '-DAF', '-lafcpu', ['-L"' af_path '/lib"'],...
                            ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'])

                        disp('Implementation 2 built with CPU support')
                    catch ME
                        if verbose
                            warning('Unable to build CPU support for implementation 2! Compiler error:');
                            disp(ME.message)
                        else
                            if strcmp(cc.Manufacturer, 'GNU')
                                warning(['Implementation 2 support with MinGW compiler requires ArrayFire to be compiled from source with MinGW. Compiler error is shown with install_mex(1).'...
                                    'If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. If you compiled ArrayFire with MinGW, add the path manually by using '...
                                    'install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'', ''C:/PATH/TO/ARRAYFIRE/V3'')'])
                            else
                                warning(['Unable to build CPU support for implementation 2! Compiler error is shown with install_mex(1). If ArrayFire has been installed in a non-standard path, it can be added manually '...
                                    'by using install_mex(0, '''', '''', ''C:/PATH/TO/ARRAYFIRE/V3'')']);
                            end
                        end
                    end
                end
            end
            if implementation == 0 || implementation == 3
                if strcmp(cc.Manufacturer, 'Microsoft')
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                elseif strcmp(cc.Manufacturer, 'Intel')
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                else
                    cxxflags = 'CXXFLAGS="$CXXFLAGS -Wno-ignored-attributes"';
                end
                try
                    %%%%%%%%%%%%%%%%%%%%%%%%% Implementation 3 %%%%%%%%%%%%%%%%%%%%%%%%
                    mex(compiler, '-largeArrayDims', '-outdir', folder, cxxflags, '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                        ['-I"' opencl_include_path '"'], [folder '/OpenCL_device_info.cpp'])

                    mex(compiler, complexFlag, '-outdir', folder, cxxflags, '-DMATLAB', '-DOPENCL', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                        ['-I"' opencl_include_path '"'], [folder '/OpenCL_matrixfree_multi_gpu.cpp'])


                    disp('Implementation 3 built')
                catch ME
                    if verbose
                        warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                            'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                            'install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64''). Compiler error:']);
                        disp(ME.message)
                    else
                        warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                            'Compiler error is shown with install_mex(1). If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                            'install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'')']);
                    end
                end
            end
        end
    else
        %%% Linux/MAC %%%
        if implementation == 0 || implementation == 2
            if use_CUDA
                compflags = 'COMPFLAGS="$COMPFLAGS -std=c++17"';
                cxxflags = 'CXXFLAGS="$CXXFLAGS -w"';
                try
                    disp('Attempting to build CUDA code.')
                    mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree', ldflags, compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', ['-I' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', '-lnvrtc', ...
                        ['-L"' af_path '/lib64"'], ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib64"'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);

                    mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree_uint16', ldflags, compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE', ['-I' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', '-lnvrtc', ...
                        ['-L"' af_path '/lib64"'], ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib64"'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);

                    mex(compiler, complexFlag, '-outdir', folder, '-output', 'CUDA_matrixfree_uint8', ldflags, compflags, cxxflags, '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE2', ['-I' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', '-lnvrtc', ...
                        ['-L"' af_path '/lib64"'], ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib64"'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);

                    disp('CUDA support enabled')
                catch ME
                    warning('CUDA support not enabled')
                    if verbose
                        disp(ME.message)
                    end
               end
            end
        end
        if strcmp(cc.Manufacturer, 'Microsoft')
            compflags = 'COMPFLAGS="$COMPFLAGS"';
            cxxflags = 'CXXFLAGS="$CXXFLAGS"';
        elseif strcmp(cc.Manufacturer, 'Intel')
            %             if ispc
            compflags = 'COMPFLAGS="$COMPFLAGS"';
            %             else
            %                 compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
            %             end
            cxxflags = 'CXXFLAGS="$CXXFLAGS"';
        else
            compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
            cxxflags = 'CXXFLAGS="$CXXFLAGS -Wno-ignored-attributes -Wno-narrowing"';
        end
        cxxlib = '-lOpenCL';
        if implementation == 0 || implementation == 2
            try
                mex(compiler, '-largeArrayDims', '-outdir', folder, ldflags, '-lafopencl', cxxlib, ['-L' af_path '/lib64'], ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib64"'], ...
                    ['-L"' opencl_lib_path '"'], '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64' ,['-I ' folder], ['-I' af_path_include], ...
                    ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/ArrayFire_OpenCL_device_info.cpp'])

                mex(compiler, complexFlag, '-outdir', folder, ldflags, compflags, cxxflags, '-DMATLAB', '-DOPENCL', '-DAF', '-lafopencl', cxxlib, ['-L' af_path '/lib64'], ['-L"' af_path '/lib"'], ...
                    ['-L"' cuda_path '/lib64"'], ['-L"' opencl_lib_path '"'], '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ...
                    ['-I' af_path_include], ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', ...
                    [folder '/OpenCL_matrixfree.cpp'])

                mex(compiler, complexFlag, '-outdir', folder, '-output', 'OpenCL_matrixfree_uint16', ldflags, compflags, cxxflags, '-DMATLAB', '-DOPENCL', '-DAF', '-DMTYPE', '-lafopencl', cxxlib, ['-L' af_path '/lib64'], ['-L"' af_path '/lib"'], ...
                    ['-L"' cuda_path '/lib64"'], ['-L"' opencl_lib_path '"'], '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ...
                    ['-I' af_path_include], ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', ...
                    [folder '/OpenCL_matrixfree.cpp'])

                mex(compiler, complexFlag, '-outdir', folder, '-output', 'OpenCL_matrixfree_uint8', ldflags, compflags, cxxflags, '-DMATLAB', '-DOPENCL', '-DAF', '-DMTYPE', '-lafopencl', cxxlib, ['-L' af_path '/lib64'], ['-L"' af_path '/lib"'], ...
                    ['-L"' cuda_path '/lib64"'], ['-L"' opencl_lib_path '"'], '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ...
                    ['-I' af_path_include], ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', ...
                    [folder '/OpenCL_matrixfree.cpp'])
                disp('Implementation 2 built')
            catch ME
                if verbose
                    warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                        'If OpenCL SDK and/or ArrayFire has been installed in a non-standard path, they can be added manually '...
                        'by using install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'', ''/PATH/TO/ARRAYFIRE''). Compiler error:']);
                    disp(ME.message)
                else
                    warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                        'Compiler error is shown with install_mex(1). If OpenCL SDK and/or ArrayFire has been installed in a non-standard path, they can be added manually '...
                        'by using install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'', ''/PATH/TO/ARRAYFIRE'')']);
                end
            end
            if strcmp(cc.Manufacturer, 'Microsoft')
                compflags = 'COMPFLAGS="$COMPFLAGS /openmp -DCPU"';
                cxxflags = 'CXXFLAGS="$CXXFLAGS"';
            elseif strcmp(cc.Manufacturer, 'Intel')
                if ispc
                    compflags = 'COMPFLAGS="$COMPFLAGS /Qopenmp -DCPU"';
                else
                    compflags = 'COMPFLAGS="$COMPFLAGS -qopenmp -DCPU"';
                end
                cxxflags = 'CXXFLAGS="$CXXFLAGS"';
            else
                compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
                cxxflags = 'CXXFLAGS="$CXXFLAGS -fopenmp -DCPU -Wno-ignored-attributes"';
            end
            try
                mex(compiler, complexFlag, '-outdir', folder, '-output', 'CPU_matrixfree', ['-L' OMPPath], OMPh, OMPLib, LPLib, ['-I ' folder], ldflags, compflags, cxxflags, '-DMATLAB', '-DCPU', '-DAF', '-lafcpu', ['-L' af_path '/lib64'], ['-L"' af_path '/lib"'], ...
                    ['-I ' folder], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp'])
                disp('Implementation 2 built with CPU support')
            catch ME
                if verbose
                    warning(['Unable to build CPU support for implementation 2! If rrayFire has been installed in a non-standard path, it can be added manually '...
                        'by using install_mex(0, '''', '''', ''/PATH/TO/ARRAYFIRE''). Compiler error:']);
                    disp(ME.message)
                else
                    warning(['Unable to build CPU support for implementation 2! Compiler error is shown with install_mex(1). If ArrayFire has been installed in a non-standard path, it can be added manually '...
                        'by using install_mex(0, '''', '''', ''/PATH/TO/ARRAYFIRE'')']);
                end
            end
        end
        if implementation == 0 || implementation == 3
            if strcmp(cc.Manufacturer, 'Microsoft')
                cxxflags = 'CXXFLAGS="$CXXFLAGS"';
            elseif strcmp(cc.Manufacturer, 'Intel')
                cxxflags = 'CXXFLAGS="$CXXFLAGS"';
            else
                cxxflags = 'CXXFLAGS="$CXXFLAGS -Wno-ignored-attributes"';
            end
            try
                mex(compiler, '-largeArrayDims', '-outdir', folder, ldflags, cxxflags, cxxlib, ['-L"' cuda_path '/lib64"'], ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                    ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_device_info.cpp'])

                mex(compiler, complexFlag, '-outdir', folder, ldflags, cxxflags, '-DMATLAB', '-DOPENCL', cxxlib, ['-L"' cuda_path '/lib64"'], ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                    ['-I"' cuda_path '/include"'], ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree_multi_gpu.cpp'])
                disp('Implementation 3 and 5 built')
            catch ME
                if verbose
                    warning(['Unable to build OpenCL files for implementation 3 and 5! If you do not need implementation 3 or 5 (matrix free OpenCL) ignore this warning. '...
                        'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                        'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64''). Compiler error:']);
                    disp(ME.message)
                else
                    warning(['Unable to build OpenCL files for implementation 3 and 5! If you do not need implementation 3 or 5 (matrix free OpenCL) ignore this warning. '...
                        'Compiler error is shown with install_mex(1). If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                        'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'')']);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCTAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~ismac
    joku = getenv('CXXFLAGS');
    cxxflags = '-std=c++11 -fopenmp -fPIC';
    OMPlib = '-lgomp';
    if ~any(strfind(joku,'-fopenmp'))
        cxxflags = [cxxflags ' ', joku];
        setenv('CXXFLAGS',cxxflags);
    end
    jokuL = getenv('LDFLAGS');
    setenv('LDFLAGS','-fopenmp');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementations 1 & 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, sys] = mkoctfile('--mex','-DMATLAB', ['-I ' folder], OMPlib, ...
        [folder '/projector_mex.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp']);
    if sys == 0
        movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
    end
    [~, sys2] = mkoctfile('--mex','-DMATLAB', ['-I ' folder], OMPlib, [folder '/projector_mexSingle.cpp']);
    if sys2 == 0 && sys == 0
        movefile('projector_mexSingle.mex', [folder '/projector_mexSingle.mex'],'f');
        disp('Implementations 1 & 4 built with OpenMP (parallel) support')
    end
    setenv('CXXFLAGS',joku);
    setenv('LDFLAGS',jokuL);
    if sys ~= 0 || sys2 ~= 0
        if verbose
            mkoctfile('--mex','-DMATLAB', ['-I ' folder], ...
                [folder '/projector_mex.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'])
            movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
            warning('Implementations 1 & 4 built WITHOUT OpenMP (parallel) support.')
        else
            mkoctfile('--mex','-DMATLAB', ['-I ' folder], '"--std=c++11"', ...
                [folder '/projector_mex.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'])
            movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
            warning('Implementations 1 & 4 built WITHOUT OpenMP (parallel) support, Use install_mex(1) to see compiler error')
        end
    end
    if ~any(strfind(joku,'-fopenmp'))
        cxxflags = [cxxflags ' ', joku];
        setenv('CXXFLAGS',cxxflags);
    end
    setenv('LDFLAGS','-fopenmp');
    [~, sys] = mkoctfile('--mex',['-I ' folder], OMPlib, [folder '/NLM_func.cpp']);
    [~, sys2] = mkoctfile('--mex',['-I ' folder], OMPlib, [folder '/NLM_funcSingle.cpp']);
    setenv('CXXFLAGS',joku);
    setenv('LDFLAGS',jokuL);
    if sys == 0 && sys2 == 0
        movefile('NLM_func.mex', [folder '/NLM_func.mex'],'f');
        movefile('NLM_funcSingle.mex', [folder '/NLM_funcSingle.mex'],'f');
    else
        [~, sys] = mkoctfile('--mex',['-I ' folder], [folder '/NLM_func.cpp']);
        [~, sys2] = mkoctfile('--mex',['-I ' folder], [folder '/NLM_funcSingle.cpp']);
        if sys == 0 && sys2 == 0
            movefile('NLM_func.mex', [folder '/NLM_func.mex'],'f');
            movefile('NLM_funcSingle.mex', [folder '/NLM_funcSingle.mex'],'f');
            warning('NLM support built WITHOUT OpenMP (parallel) support.')
        else
            if verbose
                warning('NLM support for implementations 1 and 4 not enabled. Compiler error: ')
            else
                warning('NLM support for implementations 1 and 4 not enabled. Use install_mex(1) to see compiler error.')
            end
        end
    end
    if ~any(strfind(joku,'-fopenmp'))
        cxxflags = [cxxflags ' ', joku];
        setenv('CXXFLAGS',cxxflags);
    end
    setenv('LDFLAGS','-fopenmp');
    [~, sys] = mkoctfile(['-I' folder], OMPlib, [folder '/createSinogramASCIIOct.cpp']);
    setenv('CXXFLAGS',joku);
    setenv('LDFLAGS',jokuL);
    if sys == 0
        movefile('createSinogramASCIIOct.oct', [folder '/createSinogramASCIIOct.oct'],'f');
    else
        [~, sys] = mkoctfile(['-I' folder], [folder '/createSinogramASCIIOct.cpp']);
        if sys == 0
            movefile('createSinogramASCIIOct.oct', [folder '/createSinogramASCIIOct.oct'],'f');
            warning('ASCII sinogram creation built WITHOUT OpenMP (parallel) support.')
        else
            if verbose
                warning('ASCII sinogram creation not enabled. Compiler error: ')
            else
                warning('ASCII sinogram creation not enabled. Use install_mex(1) to see compiler error.')
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMF support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mkoctfile('--mex', [folder '/gate_lmf_matlab.cpp'])
    % movefile('gate_lmf_matlab.mex', [folder '/gate_lmf_matlab.mex'],'f');
    % disp('LMF support enabled')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inveon support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkoctfile('--mex', ['-I ' folder], [folder '/inveon_list2matlab.cpp'])
    movefile('inveon_list2matlab.mex', [folder '/inveon_list2matlab.mex'],'f');
    disp('Inveon support enabled')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROOT support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc
      warning('ROOT data load is not supported on Octave on Windows!')
%        try
%            r_path = root_path;
%            if strcmp(r_path,'root-config')
%                r_path = 'C:\root\';
%            else
%                r_path = root_path(1:end-15);
%            end
%            mkoctfile('-v', '-MD', '-lCore', '-lTree', ['-I' r_path 'include'], ['-L' r_path '/lib'],...
%                [folder '/GATE_root_matlab_oct.cpp'])
%            movefile('GATE_root_matlab_oct.oct', [folder '/GATE_root_matlab_oct.oct'],'f');
%            disp('ROOT support enabled')
%        catch ME
%            if verbose
%                warning('Unable to build ROOT support! Compiler error: ')
%                disp(ME.message);
%            else
%                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise, make sure that ROOT was compiled with the same compiler you are using to compile '...
%                    'this mex-file and provide the path to ROOT install folder with install_mex(0, [], [], [], ''C:\path\to\ROOT''). Compiler error is shown with install_mex(1)']);
%            end
%        end
    else
        try
            mkoctfile(['"$(' root_path '  --cflags)"'], '-DOCTAVE', ['"-Wl,-lCore -lTree -lRIO -ldl $(' root_path ' --libs)"'], ...
                [folder '/GATE_root_matlab_oct.cpp'])
            movefile('GATE_root_matlab_oct.oct', [folder '/GATE_root_matlab_oct.oct'],'f');
            disp('ROOT support enabled')
        catch ME
            if verbose
                warning(['Unable to build ROOT support! Make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and that root-config '...
                    'is working properly, i.e. you have sourced the ROOT binaries with thisroot.*. You can manually input the ROOT installation folder with ' ...
                    'install_mex(0, [], [], [], ''/path/to/ROOT''). Compiler error: '])
                disp(ME.message);
            else
                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
                    'the same compiler you are using to compile this mex-file and that root-config is working properly, i.e. you have sourced the ROOT binaries with thisroot.*. '...
                    'You can manually input the ROOT installation folder with install_mex(0, [], [], [], ''/path/to/ROOT''). Compiler error is shown with install_mex(1)']);
            end
        end
    end

    % If you get any compilation errors then change the paths if
    % necessary. The examples below are for CUDA, AMD APP SDK (untested)
    % OCL SDK, and Intel SDK.

    %%%%%%%%%%%%%%%%%%%%%%%%% Implementations 2, 3 & 5 %%%%%%%%%%%%%%%%%%%%%%%%
    if (ispc)
        %%% Windows %%%
        warning('ArrayFire needs to be built with MinGW in order for implementation 2 to work. You cannot use the prebuilt binaries.')
        if isempty(opencl_include_path)
            warning(['No OpenCL SDK found. Implementations 2 and 3 will not be built. Use install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'') '...
                'to set OpenCL include and library paths.'])
        else
            warning('Implemention 2 is ONLY supported with Octave on Windows when ArrayFire has been built from source with MinGW-w64.')
            if isempty(af_path)
                warning('ArrayFire not found on path. Implementation 2 will not be built. Use install_mex(1, [], [], ''C:/PATH/TO/ARRAYFIRE/'') to set ArrayFire path.')
            else
                if use_CUDA
                    %                    disp('Attemping to build CUDA code.')
                    %                    [~, sys] = mkoctfile('--mex', '-w', ['-I"' folder '"'], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                    %                        '-lnvrtc', ['-L"' af_path '/lib64"'], ['-L"' af_path '/lib"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path_include '"'], [folder '/CUDA_matrixfree.cpp'],...
                    %                        [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                    %                        [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter_CUDA.cpp'], [folder '/mexFunktio.cpp']);
                    syst = 1;
                    %                    if sys ~= 0
                    %                        charArray = mkoctfile('--mex', '-v', '-w', ['-I' folder], ['-I' cuda_path '/include'], '-lafcuda', '-lcuda', ...
                    %                            '-lnvrtc', ['-L' af_path '/lib64'], ['-L' af_path '/lib'], ['-L' cuda_path '/lib/x64'], ['-I' af_path_include], [folder '/CUDA_matrixfree.cpp'],...
                    %                            [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                    %                            [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter_CUDA.cpp'], [folder '/mexFunktio.cpp']);
                    %                        sys = makeOCT(charArray);
                    %                    end
                    %                    syst = syst + sys;
                    if syst == 0
                        disp('CUDA support enabled')
                    else
                        warning('CUDA support not supported on Octave on Windows')
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%% Implementation 2 %%%%%%%%%%%%%%%%%%%%%%%
                [~, sys] = mkoctfile('--mex', '-DMATLAB', '-DAF', '-DOPENCL', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ...
                    ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path ''], ['-I' folder ''], ['-I' opencl_include_path ''], ['-I' af_path_include], ...
                    [folder '/OpenCL_matrixfree.cpp']);
                syst = 0;
                if sys ~= 0
                    charArray = mkoctfile('--mex', '-v', '-DMATLAB', '-DAF', '-DOPENCL', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path], ...
                        ['-I' folder], ['-I' opencl_include_path], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);
                    sys = makeOCT(charArray);
                end
                syst = syst + sys;

                [~, sys] = mkoctfile('--mex', '-DOPENCL', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path], ...
                    ['-I ' folder], ['-I' opencl_include_path], ['-I' af_path_include], [folder '/ArrayFire_OpenCL_device_info.cpp']);
                if sys ~= 0
                    charArray = mkoctfile('--mex', '-v', '-DOPENCL', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path], ...
                        ['-I ' folder], ['-I' opencl_include_path], ['-I' af_path_include], [folder '/ArrayFire_OpenCL_device_info.cpp']);
                    sys = makeOCT(charArray);
                end
                syst = syst + sys;
                if syst == 0
                    movefile('ArrayFire_OpenCL_device_info.mex', [folder '/ArrayFire_OpenCL_device_info.mex'],'f');
                    movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree.mex'],'f');

                    [~, sys] = mkoctfile('--mex', '-DMATLAB', '-DAF', '-DOPENCL', '-DMTYPE', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ...
                        ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path ''], ['-I' folder ''], ['-I' opencl_include_path ''], ['-I' af_path_include], ...
                        [folder '/OpenCL_matrixfree.cpp']);
                    syst = 0;
                    if sys ~= 0
                        charArray = mkoctfile('--mex', '-v', '-DMATLAB', '-DAF', '-DOPENCL', '-DMTYPE', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path], ...
                            ['-I' folder], ['-I' opencl_include_path], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);
                        sys = makeOCT(charArray);
                    end
                    movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree_uint16.mex'],'f');
                    syst = syst + sys;

                    [~, sys] = mkoctfile('--mex', '-DMATLAB', '-DAF', '-DOPENCL', '-DMTYPE2', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ...
                        ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path ''], ['-I' folder ''], ['-I' opencl_include_path ''], ['-I' af_path_include], ...
                        [folder '/OpenCL_matrixfree.cpp']);
                    if sys ~= 0
                        charArray = mkoctfile('--mex', '-v', '-DMATLAB', '-DAF', '-DOPENCL', '-DMTYPE2', '-Wno-ignored-attributes', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib64'], ['-L' af_path '\lib'],['-L' opencl_lib_path], ...
                            ['-I' folder], ['-I' opencl_include_path], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);
                        sys = makeOCT(charArray);
                    end
                    movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree_uint8.mex'],'f');
                    syst = syst + sys;
                    if syst == 0
                        disp('Implementation 2 built')
                    else
                        disp('Implementation 2 built without high-resolution support!')
                    end
                else
                    warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                        'If OpenCL SDK and/or ArrayFire has been installed in a non-standard path, they can be added manually '...
                        'by using install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'', ''C:/PATH/TO/ARRAYFIRE'')']);
                end


                setenv('CXXFLAGS',cxxflags);
                setenv('LDFLAGS','-fopenmp');
                [~, sys] = mkoctfile('--mex', '-DMATLAB', '-DAF', '-DCPU', '-Wno-ignored-attributes', '-lafcpu', ...
                    ['-L' af_path '\lib64'], ['-L' af_path '\lib'], ['-I' folder ''], ['-I' af_path_include], ...
                    [folder '/OpenCL_matrixfree.cpp']);
                syst = 0;
                if sys ~= 0
                    charArray = mkoctfile('--mex', '-v', '-DMATLAB', '-DAF', '-DCPU', '-Wno-ignored-attributes', '-lafcpu', ['-L' af_path '\lib64'], ['-L' af_path '\lib'], ...
                        ['-I' folder], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp']);
                    sys = makeOCT(charArray);
                end
                syst = syst + sys;
                if syst == 0
                    movefile('OpenCL_matrixfree.mex', [folder '/CPU_matrixfree.mex'],'f');
                    disp('Implementation 2 built with CPU support')
                else
                    warning(['Unable to build CPU support for implementation 2! If ArrayFire has been installed in a non-standard path, it can be added manually '...
                        'by using install_mex(0, '''', '''', ''C:/PATH/TO/ARRAYFIRE'')']);
                end
                setenv('CXXFLAGS',joku);
                setenv('LDFLAGS',jokuL);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% Implementation 3 %%%%%%%%%%%%%%%%%%%%%%%%
            [~, sys] = mkoctfile('--mex', '-DOPENCL', '-Wno-ignored-attributes', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I' folder], ...
                ['-I' opencl_include_path], [folder '/OpenCL_device_info.cpp']);
            syst = 0;
            if sys ~= 0
                [charArray, ~] = mkoctfile('--mex', '-v', '-DOPENCL', '-Wno-ignored-attributes', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I' folder], ...
                    ['-I' opencl_include_path], [folder '/OpenCL_device_info.cpp']);
                sys = makeOCT(charArray);
            end
            syst = syst + sys;

            [~, sys] = mkoctfile('--mex', '-DMATLAB', '-DOPENCL', '-Wno-ignored-attributes', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I' folder], ...
                ['-I' opencl_include_path], [folder '/OpenCL_matrixfree_multi_gpu.cpp']);
            if sys ~= 0
                [charArray, ~] = mkoctfile('--mex', '-v', '-DMATLAB', '-DOPENCL', '-Wno-ignored-attributes', '-lOpenCL', ['-L' opencl_lib_path], ['-I' folder], ...
                    ['-I' opencl_include_path], [folder '/OpenCL_matrixfree_multi_gpu.cpp']);
                sys = makeOCT(charArray);
            end
            syst = syst + sys;

            if syst == 0
                movefile('OpenCL_device_info.mex', [folder '/OpenCL_device_info.mex'],'f');
                movefile('OpenCL_matrixfree_multi_gpu.mex', [folder '/OpenCL_matrixfree_multi_gpu.mex'],'f');
                disp('Implementation 3 and 5 built')
            else
                warning(['Unable to build OpenCL files for implementation 3 and 5! If you do not need implementation 3 and 5 (matrix free OpenCL) ignore this warning. ' ...
                    'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''C:/PATH/TO/OPENCL/INCLUDE'', ''C:/PATH/TO/OPENCL/LIB/X64'').)']);
            end
        end
    else
        cxxflags = '-Wno-ignored-attributes';
        %%% Linux/MAC %%%
        if use_CUDA
            try
                disp('Attempting to build CUDA code.')
                mkoctfile('--mex', '-DMATLAB', '-DCUDA', '-DAF', '-lafcuda', '-lnvrtc', ['-L' af_path '/lib64'], ['-L' af_path '/lib'],...
                    ['-L' cuda_path '/lib64'], ['-I ' folder], ['-I' cuda_path '/include'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp'])
                movefile('OpenCL_matrixfree.mex', [folder '/CUDA_matrixfree.mex'],'f');
                mkoctfile('--mex', '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE', '-lafcuda', '-lnvrtc', ['-L' af_path '/lib64'], ['-L' af_path '/lib'],...
                    ['-L' cuda_path '/lib64'], ['-I ' folder], ['-I' cuda_path '/include'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp'])
                movefile('OpenCL_matrixfree.mex', [folder '/CUDA_matrixfree_uint16.mex'],'f');
                mkoctfile('--mex', '-DMATLAB', '-DCUDA', '-DAF', '-DMTYPE2', '-lafcuda', '-lnvrtc', ['-L' af_path '/lib64'], ['-L' af_path '/lib'],...
                    ['-L' cuda_path '/lib64'], ['-I ' folder], ['-I' cuda_path '/include'], ['-I' af_path_include], [folder '/OpenCL_matrixfree.cpp'])
                movefile('OpenCL_matrixfree.mex', [folder '/CUDA_matrixfree_uint8.mex'],'f');
                disp('CUDA support enabled.')
            catch
                warning('CUDA support not enabled')
            end
        end

        setenv('CXXFLAGS',cxxflags);
        setenv('LDFLAGS','-fopenmp');
        try
            mkoctfile('--mex', '-DMATLAB', '-DAF', '-DCPU', cxxflags, '-lafcpu', ['-L' af_path '/lib64'], ['-L' af_path '/lib'], ['-I ' folder], ['-I' af_path_include], ...
                [folder '/OpenCL_matrixfree.cpp']);
            movefile('OpenCL_matrixfree.mex', [folder '/CPU_matrixfree.mex'],'f');
            disp('Implementation 2 built with CPU support')
        catch ME
            if verbose
                warning(['Unable to build CPU support for implementation 2! If ArrayFire has been installed in a non-standard path, it can be added manually '...
                    'by using install_mex(0, '''', '''', ''/PATH/TO/ARRAYFIRE''). Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build CPU support for implementation 2! Compiler error is shown with install_mex(1). If ArrayFire has been installed in a non-standard path, it can be added manually '...
                    'by using install_mex(0, '''', '''', ''/PATH/TO/ARRAYFIRE'')']);
            end
        end
        setenv('CXXFLAGS',joku);
        setenv('LDFLAGS',jokuL);
        try
            mkoctfile('--mex', cxxflags, '-lafopencl', ['-L' af_path '/lib64'], ['-L' af_path '/lib'], ['-L' cuda_path '/lib64'], ['-L' opencl_lib_path], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64' ,['-I ' folder], ['-I' af_path_include], ['-I' cuda_path '/include'], ['-I' opencl_include_path], ...
                '-I/opt/AMDAPPSDK-3.0/include', [folder '/ArrayFire_OpenCL_device_info.cpp'])
            movefile('ArrayFire_OpenCL_device_info.mex', [folder '/ArrayFire_OpenCL_device_info.mex'],'f');

            mkoctfile('--mex', '-DMATLAB', '-DAF', '-DOPENCL', '-DMTYPE', cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], ['-L' af_path '/lib'], ['-L' cuda_path '/lib64'], ['-L' opencl_lib_path], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ['-I' af_path_include], ['-I' cuda_path '/include'], ['-I' opencl_include_path], ...
                '-I/opt/AMDAPPSDK-3.0/include', '-DOPENCL', [folder '/OpenCL_matrixfree.cpp']);
            movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree_uint16.mex'],'f');

            mkoctfile('--mex', '-DMATLAB', '-DAF', '-DOPENCL', cxxflags, '-lOpenCL', '-lafopencl', ['-L' af_path '/lib64'], ['-L' af_path '/lib'], ['-L' cuda_path '/lib64'], ['-L' opencl_lib_path], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ['-I' af_path_include], ['-I' cuda_path '/include'], ['-I' opencl_include_path], ...
                '-I/opt/AMDAPPSDK-3.0/include', '-DOPENCL', [folder '/OpenCL_matrixfree.cpp']);

            movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree.mex'],'f');
            disp('Implementation 2 built')
        catch ME
            if verbose
                warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                    'Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                    'Compiler error is shown with install_mex(1). If OpenCL SDK and/or ArrayFire has been installed in a non-standard path, they can be added manually '...
                    'by using install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'', ''/PATH/TO/ARRAYFIRE'')']);
            end
        end
        try
            mkoctfile('--mex', '-DOPENCL', cxxflags, '-lOpenCL', ['-L' cuda_path '/lib64'], ['-L' opencl_lib_path], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                ['-I' cuda_path '/include'], ['-I' opencl_include_path ''], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_device_info.cpp'])

            mkoctfile('--mex', '-DMATLAB', '-DOPENCL', cxxflags, '-lOpenCL', ['-L' cuda_path '/lib64'], ['-L' opencl_lib_path], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                ['-I' cuda_path '/include'], ['-I' opencl_include_path], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree_multi_gpu.cpp'])
            movefile('OpenCL_device_info.mex', [folder '/OpenCL_device_info.mex'],'f');
            movefile('OpenCL_matrixfree_multi_gpu.mex', [folder '/OpenCL_matrixfree_multi_gpu.mex'],'f');
            disp('Implementation 3 and 5 built')
        catch ME
            if verbose
                warning(['Unable to build OpenCL files for implementation 3 and 5! If you do not need implementation 3 and 5 (matrix free OpenCL) ignore this warning. '...
                    'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64''). Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build OpenCL files for implementation 3 and 5! If you do not need implementation 3 and 5 (matrix free OpenCL) ignore this warning. '...
                    'Compiler error is shown with install_mex(1). If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'')']);
            end
        end
    end
else % Apple silicon MATLAB/Octave
    folderMetal = strcat(folder, '/metal');
    compiler = '';
    complexFlag = '';
    ldflags = 'LDFLAGS="\$LDFLAGS -framework Metal -framework Foundation"';
    cxxflags = 'CXXFLAGS="\$CXXFLAGS -std=c++17 -fobjc-arc"';
    mex(compiler, complexFlag, '-outdir', folderMetal, ldflags, cxxflags, '-DMATLAB', ['-I' folder], ['-I' folderMetal], [folderMetal '/Metal_matrixfree.mm'])
    disp('Implementation 2 built')
end
