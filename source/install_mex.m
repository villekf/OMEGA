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
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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
if numvarargs > 6
    error('Requires at most six input parameters')
end
if ispc
    if exist('OCTAVE_VERSION','builtin') == 5
        proggis1_orig = 'Program Files';
        proggis2_orig = 'Program Files (x86)';
        proggis1 = 'PROGRA~1';
        proggis2 = 'PROGRA~2';
        
        nvidia = 'NVIDIA~2';
        nvidia_orig = 'NVIDIA GPU Computing Toolkit';
        amd = 'AMDAPP~1';
        amd_orig = 'AMD APP SDK';
        root_path = 'C:/PROGRA~2/root/';
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
            opencl_include_path = getenv('INTELOCLSDKROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('CUDA_PATH'))
            opencl_include_path = getenv('CUDA_PATH');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_include_path = strrep(opencl_include_path, nvidia_orig, nvidia);
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('OCL_ROOT'))
            opencl_include_path = getenv('OCL_ROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('AMDAPPSDKROOT'))
            opencl_include_path = getenv('AMDAPPSDKROOT');
            opencl_include_path = strrep(opencl_include_path, proggis2_orig, proggis2);
            opencl_include_path = strrep(opencl_include_path, proggis1_orig, proggis1);
            opencl_include_path = strrep(opencl_include_path, amd_orig, amd);
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        else
            opencl_include_path = '';
            opencl_lib_path = '';
            warning('No OpenCL SDK detected. If one is installed, insert the paths manually by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'')')
        end
        af_path = getenv('AF_PATH');
        af_path = strrep(af_path, proggis2_orig, proggis2);
        af_path = strrep(af_path, proggis1_orig, proggis1);
    else
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
            opencl_include_path = getenv('INTELOCLSDKROOT');
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('CUDA_PATH'))
            opencl_include_path = getenv('CUDA_PATH');
            opencl_lib_path = [opencl_include_path '/lib/x64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('OCL_ROOT'))
            opencl_include_path = getenv('OCL_ROOT');
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        elseif ~isempty(getenv('AMDAPPSDKROOT'))
            opencl_include_path = getenv('AMDAPPSDKROOT');
            opencl_lib_path = [opencl_include_path '/lib/x86_64'];
            opencl_include_path = [opencl_include_path '/include'];
        else
            opencl_include_path = '';
            opencl_lib_path = '';
        end
        af_path = getenv('AF_PATH');
        root_path = 'C:/Program Files/root/';
    end
    cuda_path = getenv('CUDA_PATH');
    optargs = {false, opencl_include_path, opencl_lib_path, af_path, root_path, true, cuda_path};
else
    cuda_path = '/usr/local/cuda';
    optargs = {false, '/opt/intel/opencl/include/','/opt/intel/opencl/lib/x64/','/opt/arrayfire/', '/opt/root/', true, cuda_path};
end

if nargin > 0
    optargs{1} = varargin{1};
    if nargin > 1 && ~isempty(varargin{2})
        optargs{2} = varargin{2};
    end
    if nargin > 2 && ~isempty(varargin{3})
        optargs{3} = varargin{3};
    end
    if nargin > 3 && ~isempty(varargin{4})
        optargs{4} = varargin{4};
    end
    if nargin > 4 && ~isempty(varargin{5})
        optargs{5} = varargin{5};
    end
    if nargin > 5 && ~isempty(varargin{6})
        optargs{6} = varargin{6};
    end
    if nargin > 6 && ~isempty(varargin{7})
        optargs{7} = varargin{7};
    end
end

% use_root = optargs{1};
% use_opencl = optargs{2};
verbose = optargs{1};
opencl_include_path = optargs{2};
opencl_lib_path = optargs{3};
af_path = optargs{4};
root_path = optargs{5};
use_CUDA = optargs{6};
cuda_path = optargs{7};

if isempty(opencl_include_path) || isempty(opencl_lib_path)
    warning('No OpenCL SDK detected. If one is installed, insert the paths manually by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'')')
end

folder = fileparts(which('install_mex.m'));
folder = strrep(folder, '\','/');

opencl_include_path = strrep(opencl_include_path, '\','/');
opencl_lib_path = strrep(opencl_lib_path, '\','/');
af_path = strrep(af_path, '\','/');
root_path = strrep(root_path, '\','/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATLAB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('OCTAVE_VERSION','builtin') == 0
    cc = mex.getCompilerConfigurations('C++','Selected');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementations 1 & 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     try
    %         mex('-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11 /openmp /Qopenmp -qopenmp"', 'CXXFLAGS="$CXXFLAGS -fopenmp"', ...
    %             'LDFLAGS="$LDFLAGS -fopenmp"', [folder '/improved_Siddon_algorithm_discard.cpp'], [folder '/projector_functions.cpp'])
    %         disp('Precomputation for implementations 1 & 4 built with OpenMP support')
    %     catch
    %         try
    %             mex('-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11"', ...
    %             [folder '/improved_Siddon_algorithm_discard.cpp'], ...
    %                 [folder '/projector_functions.cpp'])
    %             disp('Precomputation for implementations 1 & 4 built without OpenMP support')
    %         catch ME
    %             if verbose
    %                 disp(ME.message);
    %             else
    %                 warning('Build failed, make sure you have a C++11 compatible compiler. Use install_mex(1) to see compiler error')
    %             end
    %         end
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementations 1 & 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    try
        mex('-largeArrayDims', '-outdir', folder, ['-I ' folder], compflags, cxxflags, ...
            'LDFLAGS="$LDFLAGS -fopenmp"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], [folder '/improved_siddon_precomputed.cpp'], ...
            [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], ...
            [folder '/improved_siddon_no_precompute.cpp'], [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'],...
            [folder '/volume_projector_functions.cpp'], [folder '/vol_siddon_precomputed.cpp'])
        disp('Implementations 1 & 4 built with OpenMP (parallel) support')
    catch ME
        if verbose
            mex('-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
                [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
                [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], ...
                [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'], [folder '/volume_projector_functions.cpp'], ...
                [folder '/vol_siddon_precomputed.cpp'])
            warning('Implementations 1 & 4 built without OpenMP (parallel) support. Compiler error: ')
            disp(ME.message);
        else
            mex('-largeArrayDims', '-outdir', folder, ['-I ' folder], 'COMPFLAGS="$COMPFLAGS -std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
                [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
                [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], ...
                [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'], [folder '/volume_projector_functions.cpp'], ...
                [folder '/vol_siddon_precomputed.cpp'])
            warning('Implementations 1 & 4 built without OpenMP (parallel) support, Use install_mex(1) to see compiler error')
        end
    end
    try
        mex('-largeArrayDims', '-outdir', folder, compflags, cxxflags, 'LDFLAGS="$LDFLAGS -fopenmp"', [folder '/NLM_func.cpp'])
    catch ME
        if verbose
            warning('NLM support for implementations 1 and 4 not enabled. Compiler error: ')
            disp(ME.message);
        else
            warning('NLM support for implementations 1 and 4 not enabled. Use install_mex(1) to see compiler error.')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMF support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mex('-largeArrayDims', '-outdir', folder, [folder '/gate_lmf_matlab.cpp'])
    disp('LMF support enabled')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inveon support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mex('-largeArrayDims', '-outdir', folder, [folder '/inveon_list2matlab.cpp'])
    disp('Inveon support enabled')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROOT support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc
        try
            if verLessThan('matlab','9.6')
                mex('-largeArrayDims', '-outdir', folder, 'COMPFLAGS="$COMPFLAGS -MD -EHsc -GR"', '-lCore', '-lTree', ['-L"' root_path '/lib"'],...
                    ['-I"' root_path '/include"'], [folder '/GATE_root_matlab_C.cpp'])
            else
                mex('-largeArrayDims', '-outdir', folder, 'COMPFLAGS="$COMPFLAGS -MD -EHsc -GR"', '-lCore', '-lTree', ['-L"' root_path '/lib"'],...
                    ['-I"' root_path '/include"'], [folder '/GATE_root_matlab.cpp'])
            end
            disp('ROOT support enabled')
        catch ME
            if verbose
                warning('Unable to build ROOT support! On Windows ROOT is supported only on 32-bit MATLAB (untested). Compiler error: ')
                disp(ME.message);
            else
                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. On Windows ROOT supports only 32-bit MATLAB (untested). '...
                    'If you are using 32-bit MATLAB make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and provide the path to ROOT '...
                    'install folder with install_mex(0, [], [], [], ''C:/path/to/ROOT''). Compiler error is shown with install_mex(1)']);
            end
        end
    else
        try
            if verLessThan('matlab','9.6')
                mex('-largeArrayDims', '-outdir', folder, 'CXXFLAGS="$CXXFLAGS $(root-config --cflags)"', '-lCore', '-lTree', '-ldl', 'LDFLAGS="$LDFLAGS $(root-config --libs)"', ...
                    [folder '/GATE_root_matlab_C.cpp'])
                warning('Importing ROOT files will cause MATLAB to crash if you are not using R2019a or newer. Use MATLAB with `matlab -nojvm` to circumvent this.')
            else
                mex('-largeArrayDims', '-outdir', folder, 'CXXFLAGS="$CXXFLAGS $(root-config --cflags)"', '-lCore', '-lTree', '-ldl', 'LDFLAGS="$LDFLAGS $(root-config --libs)"', ...
                    [folder '/GATE_root_matlab.cpp'])
            end
            disp('ROOT support enabled')
        catch ME
            if verbose
                warning(['Unable to build ROOT support! Make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and provide the path '...
                    'to ROOT install folder with install_mex(1, [], [], [], ''/path/to/ROOT''). Compiler error: '])
                disp(ME.message);
            else
                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
                    'the same compiler you are using to compile this mex-file and provide the path to ROOT install folder with install_mex(0, [], [], [], ''/path/to/ROOT''). '...
                    'Compiler error is shown with install_mex(1)']);
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
            warning('No OpenCL SDK found. Implementations 2 and 3 will not be built. Use install_mex(1, ''C:/PATH/TO/OPENCL/'') to set OpenCL SDK path.')
        else
            if isempty(af_path)
                warning('ArrayFire not found on path. Implementation 2 will not be built. Use install_mex(1, [],''C:/PATH/TO/ARRAYFIRE/'') to set ArrayFire path.')
            else
                if use_CUDA
                    try
                        disp('Attemping to build CUDA code.')
                        mex('-largeArrayDims','-outdir', folder, ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                            '-lnvrtc', ['-L"' af_path '/lib64"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path '/include"'], [folder '/CUDA_matrixfree.cpp'],...
                        [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                        [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter_CUDA.cpp'])
                    
                        disp('CUDA support enabled')
                    catch
                        warning('CUDA support not enabled')
                    end
                end
                if strcmp(cc.Manufacturer, 'Microsoft')
                    compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                elseif strcmp(cc.Manufacturer, 'Intel')
                    if ispc
                        compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    else
                        compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
                    end
                    cxxflags = 'CXXFLAGS="$CXXFLAGS"';
                else
                    compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
                    cxxflags = 'CXXFLAGS="$CXXFLAGS -DOPENCL"';
                end
                try
                    %%%%%%%%%%%%%%%%%%%%%% Implementations 2 & 5 %%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%% Implementation 5 %%%%%%%%%%%%%%%%%%%%%%%
                    mex('-largeArrayDims', '-outdir', folder, compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib64"'],['-L"' opencl_lib_path '"'], ...
                        ['-I ' folder], ['-I"' opencl_include_path '"'], ['-I"' af_path '/include"'], [folder '/improved_Siddon_openCL.cpp'], ...
                        [folder '/functions.cpp'], [folder '/AF_opencl_functions.cpp'],[folder '/opencl_error.cpp'])
                end
                try
                    %%%%%%%%%%%%%%%%%%%%%% Implementation 2 %%%%%%%%%%%%%%%%%%%%%%%
                    mex('-largeArrayDims','-outdir', folder, ['-I ' folder], ['-I"' opencl_include_path '"'], compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib"'],...
                        ['-L"' opencl_lib_path '"'], ['-I"' af_path '/include"'], [folder '/OpenCL_matrixfree.cpp'],...
                        [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/precomp.cpp'], [folder '/AF_opencl_functions.cpp'], ...
                        [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], ...
                        [folder '/find_lors.cpp'], [folder '/compute_OS_estimates_subiter.cpp'])
                    
                    mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' af_path '/lib64"'],['-L"' opencl_lib_path '"'], ...
                        ['-I ' folder], ['-I"' opencl_include_path '"'], ['-I"' af_path '/include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    
                    disp('Implementation 2 built')
                catch ME
                    if verbose
                        warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                            'If OpenCL SDK and\or ArrayFire has been installed in a non-standard path, they can be added manually '...
                            'by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'', ''C:\PATH\TO\ARRAYFIRE\V3''). Compiler error:']);
                        disp(ME.message)
                    else
                        warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                            'Compiler error is shown with install_mex(1). If OpenCL SDK and\or ArrayFire has been installed in a non-standard path, they can be added manually '...
                            'by using install_mex(0, ''C:\PATH\TO\OPENCL\INCLUDE'', ''C:\PATH\TO\OPENCL\LIB\X64'', ''C:\PATH\TO\ARRAYFIRE\V3'')']);
                    end
                end
            end
            try
                %%%%%%%%%%%%%%%%%%%%%%%%% Implementation 3 %%%%%%%%%%%%%%%%%%%%%%%%
                mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                    ['-I"' opencl_include_path '"'], [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
                
                mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                    ['-I"' opencl_include_path '"'], [folder '/OpenCL_matrixfree_multi_gpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], ...
                    [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'],[folder '/precomp.cpp'],[folder '/forward_backward_projections.cpp'], ...
                    [folder '/multigpu_OSEM.cpp'])
                
                
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
    else
        %%% Linux/MAC %%%
        if use_CUDA
            try
                disp('Attemping to build CUDA code.')
                mex('-largeArrayDims','-outdir', folder, 'CXXFLAGS="$CXXFLAGS -w"', ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', '-lnvrtc', ['-L"' af_path '/lib64"'],...
                    ['-L"' cuda_path '/lib64"'], ['-I"' af_path '/include"'], [folder '/CUDA_matrixfree.cpp'],...
                    [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                    [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], ...
                    [folder '/compute_OS_estimates_subiter_CUDA.cpp']);
                
                disp('CUDA support enabled')
            catch
                warning('CUDA support not enabled')
            end
        end
        if strcmp(cc.Manufacturer, 'Microsoft')
            compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
            cxxflags = 'CXXFLAGS="$CXXFLAGS"';
        elseif strcmp(cc.Manufacturer, 'Intel')
            if ispc
                compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
            else
                compflags = 'COMPFLAGS="$COMPFLAGS -DOPENCL"';
            end
            cxxflags = 'CXXFLAGS="$CXXFLAGS"';
        else
            compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
            cxxflags = 'CXXFLAGS="$CXXFLAGS -DOPENCL"';
        end
        try
            mex('-largeArrayDims', '-outdir', folder, compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/AMDAPPSDK-3.0/lib/x86_64' , '-L/opt/amdgpu-pro/lib64',['-I ' folder], ['-I' af_path '/include'], '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', [folder '/improved_Siddon_openCL.cpp'], [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/precomp.cpp'],...
                [folder '/AF_opencl_functions.cpp'])
        end
        try
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64' ,['-I ' folder], '-I/opt/arrayfire/include', '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', [folder '/ArrayFire_OpenCL_device_info.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ['-I' af_path '/include'], '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree.cpp'], [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], ...
                [folder '/precomp.cpp'], [folder '/opencl_error.cpp'], [folder '/find_lors.cpp'], [folder '/AF_opencl_functions.cpp'], [folder '/compute_ML_estimates.cpp'], ...
                [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter.cpp'])
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
        try
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree_multi_gpu.cpp'], ...
                [folder '/functions_multigpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], [folder '/precomp.cpp'], [folder '/opencl_error.cpp'],[folder '/forward_backward_projections.cpp'], ...
                    [folder '/multigpu_OSEM.cpp'])
            disp('Implementation 3 built')
        catch ME
            if verbose
                warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                    'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64''). Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                    'Compiler error is shown with install_mex(1). If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'')']);
            end
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCTAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementations 1 & 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        mkoctfile('--mex', ['-I ' folder], '"--std=c++11"', '"-Wl,-fopenmp"', ...
            [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], ...
            [folder '/sequential_improved_siddon_openmp.cpp'], [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], ...
            [folder '/improved_siddon_no_precompute.cpp'], [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'], ...
            [folder '/volume_projector_functions.cpp'], [folder '/vol_siddon_precomputed.cpp'])
        movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
        disp('Implementations 1 & 4 built with OpenMP (parallel) support')
    catch ME
        if verbose
            mkoctfile('--mex', ['-I ' folder], '"--std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
                [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
                [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], ...
                [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'], ...
                [folder '/volume_projector_functions.cpp'], [folder '/vol_siddon_precomputed.cpp'])
            movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
            warning('Implementations 1 & 4 built without OpenMP (parallel) support. Compiler error: ')
            disp(ME.message);
        else
            mkoctfile('--mex', ['-I ' folder], '"--std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
                [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
                [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], ...
                [folder '/original_siddon_function.cpp'], [folder '/improved_Siddon_algorithm_discard.cpp'], ...
                [folder '/volume_projector_functions.cpp'], [folder '/vol_siddon_precomputed.cpp'])
            movefile('projector_mex.mex', [folder '/projector_mex.mex'],'f');
            warning('Implementations 1 & 4 built without OpenMP (parallel) support, Use install_mex(1) to see compiler error')
        end
    end
    try
        mkoctfile('--mex', ['-I ' folder], '"--std=c++11"', '"-Wl,-fopenmp"', [folder '/NLM_func.cpp'])
        movefile('NLM_func.mex', [folder '/NLM_func.mex'],'f');
    catch ME
        if verbose
            warning('NLM support for implementations 1 and 4 not enabled. Compiler error: ')
            disp(ME.message);
        else
            warning('NLM support for implementations 1 and 4 not enabled. Use install_mex(1) to see compiler error.')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMF support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkoctfile('--mex', [folder '/gate_lmf_matlab.cpp'])
    movefile('gate_lmf_matlab.mex', [folder '/gate_lmf_matlab.mex'],'f');
    disp('LMF support enabled')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inveon support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkoctfile('--mex', [folder '/inveon_list2matlab.cpp'])
    movefile('inveon_list2matlab.mex', [folder '/inveon_list2matlab.mex'],'f');
    disp('Inveon support enabled')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROOT support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc
        try
            mkoctfile('--mex', '-v', '-MD -EHsc -GR"', '"-Wl,-lCore -lTree"', ['-L ' root_path '/lib'],...
                ['-I ' root_path '/include'], [folder '/GATE_root_matlab_oct.cpp'])
            movefile('GATE_root_matlab_oct.mex', [folder '/GATE_root_matlab_oct.mex'],'f');
            disp('ROOT support enabled')
        catch ME
            if verbose
                warning('Unable to build ROOT support! On Windows ROOT is supported only on 32-bit OCTAVE (untested). Compiler error: ')
                disp(ME.message);
            else
                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. On Windows ROOT supports only 32-bit OCTAVE (untested). '...
                    'If you are using 32-bit OCTAVE make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and provide the path to ROOT '...
                    'install folder with install_mex(0, [], [], [], ''C:\path\to\ROOT''). Compiler error is shown with install_mex(1)']);
            end
        end
    else
        try
            mkoctfile('"$(root-config --cflags)"', '"-Wl,-lCore -lTree -ldl $(root-config --libs)"', ...
                [folder '/GATE_root_matlab_oct.cpp'])
            movefile('GATE_root_matlab_oct.oct', [folder '/GATE_root_matlab_oct.oct'],'f');
            disp('ROOT support enabled')
        catch ME
            if verbose
                warning(['Unable to build ROOT support! Make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and provide the path '...
                    'to ROOT install folder with install_mex(1, [], [], [] ''/path/to/ROOT''). Compiler error: '])
                disp(ME.message);
            else
                warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
                    'the same compiler you are using to compile this mex-file and provide the path to ROOT install folder with install_mex(0, [], [], [] ''/path/to/ROOT''). '...
                    'Compiler error is shown with install_mex(1)']);
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
            warning('No OpenCL SDK found. Implementations 2 and 3 will not be built. Use install_mex(1, ''C:/PATH/TO/OPENCL/'') to set OpenCL SDK path.')
        else
            warning('Implemention 2 is ONLY supported with Octave on Windows when ArrayFire has been built from source with MinGW-w64.')
            if isempty(af_path)
                warning('ArrayFire not found on path. Implementation 2 will not be built. Use install_mex(1, [],''C:/PATH/TO/ARRAYFIRE/'') to set ArrayFire path.')
            else
                if use_CUDA
                    try
                        disp('Attemping to build CUDA code.')
                        mkoctfile('--mex', ['-I ' folder], ['-I"' cuda_path '/include"'], '-lafcuda', '-lcuda', ...
                            '-lnvrtc', ['-L"' af_path '/lib64"'], ['-L"' cuda_path '/lib/x64"'], ['-I"' af_path '/include"'], [folder '/CUDA_matrixfree.cpp'],...
                        [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                        [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter_CUDA.cpp'])
                    
                        disp('CUDA support enabled')
                    catch
                        warning('CUDA support not enabled')
                    end
                end
                compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
                cxxflags = 'CXXFLAGS="$CXXFLAGS -DOPENCL"';
                try
                    %%%%%%%%%%%%%%%%%%%%%% Implementations 2 & 5 %%%%%%%%%%%%%%%%%%%%%%
                    if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
                        %%%%%%%%%%%%%%%%%%%%%% Implementation 5 %%%%%%%%%%%%%%%%%%%%%%%
                        try
                            mkoctfile('--mex', compflags, cxxflags, ['-I' af_path '/include/'], ['-I' opencl_path '/include/'], ['-L' af_path '/lib64/'], ['-L' opencl_path 'lib/x64/'], '-lafopencl', '-lOpenCL', ...
                                [folder '/improved_Siddon_openCL.cpp'])
                            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%% Implementation 2 %%%%%%%%%%%%%%%%%%%%%%%
                        mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'],['-L' opencl_path '\lib\x64'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/OpenCL_matrixfree.cpp'],...
                            [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/precomp.cpp'], [folder '/find_lors.cpp'])
                        
                        mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib/x64'], ...
                            ['-I ' folder], ['-I"' opencl_path '\include"'], ['-I' af_path '\include'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    elseif ~isempty(getenv('CUDA_PATH'))
                        try
                            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib\x64'], ...
                                ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/improved_Siddon_openCL.cpp'])
                            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
                        end
                        
                        mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib\x64'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/OpenCL_matrixfree.cpp'], ...
                            [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/precomp.cpp'], [folder '/find_lors.cpp'])
                        
                        mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib\x64'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    elseif ~isempty(getenv('AMDAPPSDKROOT'))
                        try
                            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                                ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/improved_Siddon_openCL.cpp'])
                            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
                        end
                        
                        mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/OpenCL_matrixfree.cpp'],...
                            [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/precomp.cpp'], [folder '/find_lors.cpp'])
                        
                        mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    elseif ~isempty(getenv('OCL_ROOT'))
                        try
                            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                                ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/improved_Siddon_openCL.cpp'])
                            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
                        end
                        
                        mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/OpenCL_matrixfree.cpp'],...
                            [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/precomp.cpp'], [folder '/find_lors.cpp'])
                        
                        mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-L' opencl_path '\lib'], ...
                            ['-I ' folder], ['-I' opencl_path '\include'], ['-I' af_path '\include'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    else
                        try
                            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-I ' folder], ['-I' af_path '\include'], ...
                                [folder '/improved_Siddon_openCL.cpp'])
                            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
                        end
                        
                        mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-I ' folder], ['-I' af_path '\include'], ...
                            [folder '/OpenCL_matrixfree.cpp'],[folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], ...
                            [folder '/precomp.cpp'], [folder '/find_lors.cpp'])
                        
                        mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '\lib'],['-I ' folder], ['-I' af_path '\include'], ...
                            [folder '/ArrayFire_OpenCL_device_info.cpp'])
                    end
                    movefile('ArrayFire_OpenCL_device_info.mex', [folder '/ArrayFire_OpenCL_device_info.mex'],'f');
                    movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree.mex'],'f');
                    disp('Implementation 2 built')
                catch ME
                    if verbose
                        warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                            'Compiler error:']);
                        disp(ME.message)
                    else
                        warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                            'Compiler error is shown with install_mex(1)']);
                    end
                end
            end
            try
                %%%%%%%%%%%%%%%%%%%%%%%%% Implementation 3 %%%%%%%%%%%%%%%%%%%%%%%%
                mkoctfile('--mex', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                    ['-I"' opencl_include_path '"'], [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
                
                mkoctfile('--mex', '-lOpenCL', ['-L"' opencl_lib_path '"'], ['-I ' folder], ...
                    ['-I"' opencl_include_path '"'], [folder '/OpenCL_matrixfree_multi_gpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], ...
                    [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'],[folder '/precomp.cpp'],[folder '/forward_backward_projections.cpp'], ...
                    [folder '/multigpu_OSEM.cpp'])
                
                movefile('OpenCL_device_info.mex', [folder '/OpenCL_device_info.mex'],'f');
                movefile('OpenCL_matrixfree_multi_gpu.mex', [folder '/OpenCL_matrixfree_multi_gpu.mex'],'f');
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
    else
        %%% Linux/MAC %%%        
        if use_CUDA
            try
                disp('Attemping to build CUDA code.')
                mkoctfile('--mex', '-lafcuda', '-lcuda', '-lnvrtc', ['-L"' af_path '/lib64"'],...
                    ['-L"' cuda_path '/lib64"'], ['-I ' folder], ['-I"' cuda_path '/include"'], ['-I"' af_path '/include"'], [folder '/CUDA_matrixfree.cpp'],...
                    [folder '/functions.cpp'], [folder '/reconstruction_AF_matrixfree_CUDA.cpp'], [folder '/AF_cuda_functions.cpp'], ...
                    [folder '/compute_ML_estimates.cpp'], [folder '/compute_OS_estimates_iter.cpp'], ...
                    [folder '/compute_OS_estimates_subiter_CUDA.cpp'])
                movefile('CUDA_matrixfree.mex', [folder '/CUDA_matrixfree.mex'],'f');
                disp('CUDA support enabled.')
            catch
                warning('CUDA support not enabled')
            end
        end

        compflags = 'COMPFLAGS="$COMPFLAGS -std=c++11"';
        cxxflags = 'CXXFLAGS="$CXXFLAGS -DOPENCL"';
        try
            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/AMDAPPSDK-3.0/lib/x86_64' , '-L/opt/amdgpu-pro/lib64',['-I ' folder], ['-I' af_path '/include'], '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', '"-DOPENCL"', [folder '/improved_Siddon_openCL.cpp'], [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/precomp.cpp'],...
                [folder '/AF_opencl_functions.cpp'])
            movefile('improved_Siddon_openCL.mex', [folder '/improved_Siddon_openCL.mex'],'f');
        end
        try
            mkoctfile('--mex', '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64' ,['-I ' folder], ['-I' af_path '/include'], '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', [folder '/ArrayFire_OpenCL_device_info.cpp'])
            
            mkoctfile('--mex', compflags, cxxflags, '-lafopencl', '-lOpenCL', ['-L' af_path '/lib64'], '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], ...
                '-L/opt/amdgpu-pro/lib64', '-L/opt/AMDAPPSDK-3.0/lib/x86_64', ['-I ' folder], ['-I' af_path '/include'], '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], ...
                '-I/opt/AMDAPPSDK-3.0/include', '"-DOPENCL"', [folder '/OpenCL_matrixfree.cpp'], [folder '/functions.cpp'], [folder '/precomp.cpp'], [folder '/opencl_error.cpp'], ...
                [folder '/find_lors.cpp'], [folder '/AF_opencl_functions.cpp'], [folder '/compute_ML_estimates.cpp'], ...
                [folder '/compute_OS_estimates_iter.cpp'], [folder '/compute_OS_estimates_subiter.cpp'])
            
            movefile('ArrayFire_OpenCL_device_info.mex', [folder '/ArrayFire_OpenCL_device_info.mex'],'f');
            movefile('OpenCL_matrixfree.mex', [folder '/OpenCL_matrixfree.mex'],'f');
            disp('Implementation 2 built')
        catch ME
            if verbose
                warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                    'Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build OpenCL files for implementation 2! If you do not need implementation 2 (matrix free OpenCL with ArrayFire) ignore this warning. '...
                    'Compiler error is shown with install_mex(1)']);
            end
        end
        try
            mkoctfile('--mex', '-lOpenCL', '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
            
            mkoctfile('--mex', '-lOpenCL', '-L/usr/local/cuda/lib64', ['-L"' opencl_lib_path '"'], '-L/opt/AMDAPPSDK-3.0/lib/x86_64', '-L/opt/amdgpu-pro/lib64', ...
                '-I/usr/local/cuda/include', ['-I"' opencl_include_path '"'], '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree_multi_gpu.cpp'], ...
                [folder '/functions_multigpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], [folder '/precomp.cpp'], [folder '/opencl_error.cpp'],[folder '/forward_backward_projections.cpp'], ...
                    [folder '/multigpu_OSEM.cpp'])
            movefile('OpenCL_device_info.mex', [folder '/OpenCL_device_info.mex'],'f');
            movefile('OpenCL_matrixfree_multi_gpu.mex', [folder '/OpenCL_matrixfree_multi_gpu.mex'],'f');
            disp('Implementation 3 built')
        catch ME
            if verbose
                warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                    'If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64''). Compiler error:']);
                disp(ME.message)
            else
                warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
                    'Compiler error is shown with install_mex(1). If OpenCL SDK has been installed in a non-standard path, it can be added manually by using '...
                    'install_mex(0, ''/PATH/TO/OPENCL/INCLUDE'', ''/PATH/TO/OPENCL/LIB/X64'')']);
            end
        end
    end
end