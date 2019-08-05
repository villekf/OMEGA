function install_mex(varargin)
%install_mex Build the necessary MEX-files
% This file can be used to build all the necessary MEX-files. ROOT and/or
% OpenCL MEX-files are optional and require additional dependencies. ROOT
% support requires that ROOT is installed and OpenCL MEX-files require
% OpenCL headers and (for some features) ArrayFire library. For ROOT
% support the user should either ensure that the main ROOT folder is either
% on path or specify it in the second input. The installation of these will
% be done automatically if the relevant headers and libraries are found. If
% not, a warning is displayed. 
%
% Example:
%   install_mex(1)
% INPUTS:
%   Input 1: If the value is set to 1, then the (possible) compiler errors
%   are displared. Otherwise only warnings are displayed if the
%   compilations were unsuccessful. This input is optional.
%   Input 2: The path to the ROOT installation folder. Default folder on
%   Windows platforms is C:\Program Files\root, on others /opt/root
%
% OUTPUTS:
%   Various MEX-files needed by OMEGA


numvarargs = length(varargin);
if numvarargs > 2
    error('Requires at most two input parameters')
end
if ispc
    optargs = {false, 'C:/Program Files/root'};
else
    optargs = {false, '/opt/root'};
end

optargs(1:numvarargs) = varargin;

% use_root = optargs{1};
% use_opencl = optargs{2};
verbose = optargs{1};
root_path = optargs{2};

folder = fileparts(which('install_mex.m'));
folder = strrep(folder, '\','/');

try
    mex('-largeArrayDims', '-outdir', folder, ['-I"' folder '"'], [folder '/improved_Siddon_algorithm_discard.cpp'], [folder '/projector_functions.cpp'])
catch
    try
        mex('-largeArrayDims', '-outdir', folder, ['-I"' folder '"'], 'COMPFLAGS="$COMPFLAGS -std=c++11 -w"', [folder '/improved_Siddon_algorithm_discard.cpp'], ...
            [folder '/projector_functions.cpp'])
    catch ME
        if verbose
            disp(ME.message);
        else
            warning('Build failed, make sure you have a C++11 compatible compiler. Use install_mex(1) to see compiler error')
        end
    end
end
try 
   mex('-largeArrayDims', '-outdir', folder, ['-I"' folder '"'], 'COMPFLAGS="$COMPFLAGS -std=c++11 /openmp /Qopenmp -qopenmp"', 'CXXFLAGS="$CXXFLAGS -fopenmp"', 'LDFLAGS="$LDFLAGS -fopenmp"', ...
       [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], ...
       [folder '/sequential_improved_siddon_openmp.cpp'], [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], ...
       [folder '/orth_siddon_no_precompute.cpp'], [folder '/original_siddon_function.cpp']);
    disp('Implementation 4 built with OpenMP (parallel) support')
catch ME
    if verbose
        mex('-largeArrayDims', '-outdir', folder, ['-I"' folder '"'], 'COMPFLAGS="$COMPFLAGS -std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
            [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
            [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], [folder '/orth_siddon_no_precompute.cpp'], ...
            [folder '/original_siddon_function.cpp'])
        warning('Implementation 4 built without OpenMP (parallel) support. Compiler error: ')
        disp(ME.message);
    else
        mex('-largeArrayDims', '-outdir', folder, ['-I"' folder '"'], 'COMPFLAGS="$COMPFLAGS -std=c++11"', [folder '/projector_mex.cpp'], [folder '/projector_functions.cpp'], ...
            [folder '/improved_siddon_precomputed.cpp'], [folder '/orth_siddon_precomputed.cpp'], [folder '/sequential_improved_siddon_openmp.cpp'], ...
            [folder '/sequential_improved_siddon_no_precompute_openmp.cpp'], [folder '/improved_siddon_no_precompute.cpp'], [folder '/orth_siddon_no_precompute.cpp'], ...
            [folder '/original_siddon_function.cpp'])
        warning('Implementation 4 built without OpenMP (parallel) support, Use install_mex(1) to see compiler error')
    end
end
mex('-largeArrayDims', '-outdir', folder, [folder '/gate_lmf_matlab.cpp'])

mex('-largeArrayDims', '-outdir', folder, [folder '/inveon_list2matlab.cpp'])

if ispc
    try
%         mex('-v', '-largeArrayDims', '-outdir', folder, 'COMPFLAGS="$COMPFLAGS -MD -EHsc -GR"', '-lCore', '-lTree', ['-L"' root_path '/lib"'],['-I"' root_path '/include"'], [folder '/GATE_root_matlab.cpp'])
        mex('-v', '-largeArrayDims', '-outdir', folder, 'COMPFLAGS="$COMPFLAGS -MD -EHsc -GR"', '-lCore', '-lTree', ['-L"' root_path '/lib"'],['-I"' root_path '/include"'], [folder '/root_test.cpp'])
%         mex('-v', '-largeArrayDims', '-outdir', folder, ['-L"' root_path '/lib"'], ['-I"' root_path '/include"'], '-lCore', '-lImt', '-lRIO', '-lNet', '-lHist', '-lGraf', '-lGraf3d', ...
%             '-lGpad', '-lROOTDataFrame', '-lROOTVecOps', '-lTree', '-lTreePlayer', '-lRint', '-lPostscript', '-lMatrix', '-lPhysics', '-lMathCore', '-lThread', '-lGui',...
%             '-lEve','-lEG','-lGeom','-lGed','-lRGL','-lImt','COMPFLAGS="$COMPFLAGS -MD -EHsc -GR"', [folder '/GATE_root_matlab.cpp'])

    catch ME
        if verbose
            warning(['Unable to build ROOT support! If you are using VS 2017 this is probably due to a bug in VS, otherwise make sure that ROOT was compiled with '...
            'the same compiler you are using to compile this mex-file and provide the path to ROOT install folder with install_mex(1, ''C:\path\to\ROOT''). Compiler error: '])
            disp(ME.message);
        else
            warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
            'the same compiler you are using to compile this mex-file and provide the path to ROOT install folder with install_mex(0, ''C:\path\to\ROOT''). '...
            'Compiler error is shown with install_mex(1)']);
        end
    end
else
    try
        mex('-largeArrayDims', '-outdir', folder, 'CXXFLAGS="$CXXFLAGS $(root-config --cflags)"', '-lCore', '-lTree', '-ldl', 'LDFLAGS="$LDFLAGS $(root-config --libs)"', ...
            [folder '/GATE_root_matlab.cpp'])
    catch ME
        if verbose
            warning(['Unable to build ROOT support! Make sure that ROOT was compiled with the same compiler you are using to compile this mex-file and provide the path '...
                'to ROOT install folder with install_mex(1, ''/path/to/ROOT''). Compiler error: '])
            disp(ME.message);
        else
            warning(['Unable to build ROOT support! If you do not need ROOT support ignore this warning. Otherwise make sure that ROOT was compiled with '...
            'the same compiler you are using to compile this mex-file and provide the path to ROOT install folder with install_mex(0, ''/path/to/ROOT''). '...
            'Compiler error is shown with install_mex(1)']);
        end
    end
end

% If you get any compilation errors then change the paths if
% necessary. The examples below are for CUDA, AMD APP SDK (untested)
% OCL SDK, and Intel SDK.

if (ispc)
    %%% Windows %%%
    try
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/improved_Siddon_openCL.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/OpenCL_matrixfree.cpp'],...
                [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/reconstruction_custom_matrixfree.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
        elseif ~isempty(getenv('CUDA_PATH'))
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/improved_Siddon_openCL.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/OpenCL_matrixfree.cpp'], ...
                [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/reconstruction_custom_matrixfree.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ...
                ['-I"' folder '"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
        elseif ~isempty(getenv('AMDAPPSDKROOT'))
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('AMDAPPSDKROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('AMDAPPSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/improved_Siddon_openCL.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('AMDAPPSDKROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('AMDAPPSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/OpenCL_matrixfree.cpp'],...
                [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/reconstruction_custom_matrixfree.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('AMDAPPSDKROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('AMDAPPSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
        elseif ~isempty(getenv('OCL_ROOT'))
                        mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('OCL_ROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('OCL_ROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/improved_Siddon_openCL.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('OCL_ROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('OCL_ROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/OpenCL_matrixfree.cpp'],...
                [folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], [folder '/reconstruction_custom_matrixfree.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('OCL_ROOT') '\lib"'], ...
                ['-I"' folder '"'], ['-I"' getenv('OCL_ROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], [folder '/ArrayFire_OpenCL_device_info.cpp'])
        else
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' folder '"'], ['-I"' getenv('AF_PATH') '\include"'], ...
                [folder '/improved_Siddon_openCL.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' folder '"'], ['-I"' getenv('AF_PATH') '\include"'], ...
                [folder '/OpenCL_matrixfree.cpp'],[folder '/functions.cpp'],[folder '/opencl_error.cpp'], [folder '/reconstruction_AF_matrixfree.cpp'], ...
                [folder '/reconstruction_custom_matrixfree.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' folder '"'], ['-I"' getenv('AF_PATH') '\include"'], ...
                [folder '/ArrayFire_OpenCL_device_info.cpp'])
        end
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
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' folder '"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ...
                [folder '/OpenCL_matrixfree_multi_gpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' folder '"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ...
                [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
        elseif ~isempty(getenv('CUDA_PATH'))
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' folder '"'], ['-I"' getenv('CUDA_PATH') '\include"'], ...
                [folder '/OpenCL_matrixfree_multi_gpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' folder '"'], ['-I"' getenv('CUDA_PATH') '\include"'], ...
                [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
        elseif ~isempty(getenv('AMDAPPSDKROOT'))
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('AMDAPPSDKROOT') '\lib"'], ['-I"' folder '"'], ['-I"' getenv('AMDAPPSDKROOT') '\include"'], ...
                [folder '/OpenCL_matrixfree_multi_gpu.cpp'],[folder '/multi_gpu_reconstruction.cpp'], [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('AMDAPPSDKROOT') '\lib"'], ['-I"' folder '"'], ['-I"' getenv('AMDAPPSDKROOT') '\include"'], ...
                [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
        elseif ~isempty(getenv('OCL_ROOT'))
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('OCL_ROOT') '\lib"'], ['-I"' folder '"'], ['-I"' getenv('OCL_ROOT') '\include"'], ...
                [folder '/OpenCL_matrixfree_multi_gpu.cpp'],[folder '/multi_gpu_reconstruction.cpp'], [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', ['-L"' getenv('OCL_ROOT') '\lib"'], ['-I"' folder '"'], ['-I"' getenv('OCL_ROOT') '\include"'], ...
                [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
        else
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', [folder '/OpenCL_matrixfree_multi_gpu.cpp'],[folder '/multi_gpu_reconstruction.cpp'],...
                [folder '/functions_multigpu.cpp'],[folder '/opencl_error.cpp'])
            
            mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
        end
    catch ME
        if verbose
            warning('Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. Compiler error:');
            disp(ME.message)
        else
            warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
            'Compiler error is shown with install_mex(1)']);
        end
    end
else
    %%% Linux/MAC %%%
    try
        mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', '-L/opt/arrayfire/lib', '-L/usr/local/cuda/lib64', '-L/opt/intel/opencl/lib64', ...
            '-L/opt/AMDAPPSDK-3.0/lib' , '-L/opt/amdgpu-pro/lib64',['-I"' folder '"'], '-I/opt/arrayfire/include', '-I/usr/local/cuda/include', '-I/opt/intel/opencl/include', ...
            '-I/opt/AMDAPPSDK-3.0/include', [folder '/improved_Siddon_openCL.cpp'])
        
        mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', '-L/opt/arrayfire/lib', '-L/usr/local/cuda/lib64', '-L/opt/intel/opencl/lib64', '-L/opt/amdgpu-pro/lib64',...
            '-L/opt/AMDAPPSDK-3.0/lib', ['-I"' folder '"'], '-I/opt/arrayfire/include', '-I/usr/local/cuda/include', '-I/opt/intel/opencl/include', '-I/opt/AMDAPPSDK-3.0/include', ...
            [folder '/OpenCL_matrixfree.cpp'], [folder '/functions.cpp'], [folder '/opencl_error.cpp'])
        
        mex('-largeArrayDims', '-outdir', folder, '-lafopencl', '-lOpenCL', '-L/opt/arrayfire/lib', '-L/usr/local/cuda/lib64', '-L/opt/intel/opencl/lib64', '-L/opt/amdgpu-pro/lib64', ...
            '-L/opt/AMDAPPSDK-3.0/lib' ,['-I"' folder '"'], '-I/opt/arrayfire/include', '-I/usr/local/cuda/include', '-I/opt/intel/opencl/include', ...
            '-I/opt/AMDAPPSDK-3.0/include', [folder '/ArrayFire_OpenCL_device_info.cpp'])
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
        mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', '-L/usr/local/cuda/lib64', '-L/opt/intel/opencl/lib64', '-L/opt/AMDAPPSDK-3.0/lib', '-L/opt/amdgpu-pro/lib64', ...
            '-I/usr/local/cuda/include', '-I/opt/intel/opencl/include', '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_matrixfree_multi_gpu.cpp'], ...
            [folder '/functions_multigpu.cpp'], [folder '/multi_gpu_reconstruction.cpp'], [folder '/opencl_error.cpp'])
        
        mex('-largeArrayDims', '-outdir', folder, '-lOpenCL', '-L/usr/local/cuda/lib64', '-L/opt/intel/opencl/lib64', '-L/opt/AMDAPPSDK-3.0/lib', '-L/opt/amdgpu-pro/lib64', ...
            '-I/usr/local/cuda/include', '-I/opt/intel/opencl/include', '-I/opt/AMDAPPSDK-3.0/include', [folder '/OpenCL_device_info.cpp'],[folder '/opencl_error.cpp'])
    catch ME
        if verbose
            warning('Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. Compiler error:');
            disp(ME.message)
        else
            warning(['Unable to build OpenCL files for implementation 3! If you do not need implementation 3 (matrix free OpenCL) ignore this warning. '...
            'Compiler error is shown with install_mex(1)']);
        end
    end
end