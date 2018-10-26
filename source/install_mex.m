function install_mex(varargin)
%% MEX build file
% This file can be used to build all the necessary mex-files

numvarargs = length(varargin);
if numvarargs > 2
    error('Requires at most two input parameters')
end
optargs = {false false};

optargs(1:numvarargs) = varargin;

use_root = optargs{1};
use_opencl = optargs{2};

mex -largeArrayDims improved_Siddon_algorithm_array.cpp
mex -largeArrayDims improved_Siddon_algorithm_discard.cpp
mex -largeArrayDims improved_Siddon_algorithm_raw.cpp
mex -largeArrayDims improved_Siddon_algorithm.cpp
mex -largeArrayDims improved_Siddon_algorithm_array_raw.cpp
mex -largeArrayDims Siddon_algorithm_discard_L.cpp
mex -largeArrayDims gate_lmf_matlab.cpp
mex -largeArrayDims sequential_reconstruction_raw.cpp
mex -largeArrayDims sequential_reconstruction.cpp

if use_root
    if isunix
        mex -largeArrayDims CXXFLAGS='$CXXFLAGS $(root-config --cflags)' -lCore -lTree -ldl LDFLAGS='$LDFLAGS $(root-config --libs)' GATE_root_matlab.cpp
    elseif ispc
        disp('Root files are not supported on Windows currently, due to a bug in VS')
    end
end

if use_opencl
    
    % If you get any compilation errors then change the paths if
    % necessary. The examples below are for CUDA and Intel SDK.
    
    if (ispc)
        %%% Windows %%%
        if isempty(getenv('CUDA_PATH')) && ~isempty(getenv('INTELOCLSDKROOT'))
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_raw.cpp')
            
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('INTELOCLSDKROOT') '\lib\x64"'], ['-I"' getenv('INTELOCLSDKROOT') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU_raw.cpp')
        elseif ~isempty(getenv('CUDA_PATH'))
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_raw.cpp')
            
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-L"' getenv('CUDA_PATH') '\lib\x64"'], ['-I"' getenv('CUDA_PATH') '\include"'], ['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU_raw.cpp')
        else
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_raw.cpp')
            
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU.cpp')
            mex('-largeArrayDims', '-lafopencl', '-lOpenCL', ['-L"' getenv('AF_PATH') '\lib"'],['-I"' getenv('AF_PATH') '\include"'], 'improved_Siddon_openCL_matrixfree_GPU_raw.cpp')
        end
        
    elseif (isunix)
        %%% Linux %%%
        mex -largeArrayDims -lafopencl -lOpenCL -L/opt/arrayfire/lib -L/usr/local/cuda/lib64 -L/opt/intel/opencl/lib64 -I/opt/arrayfire/include -I/usr/local/cuda/include improved_Siddon_openCL.cpp
        mex -largeArrayDims -lafopencl -lOpenCL -L/opt/arrayfire/lib -L/usr/local/cuda/lib64 -L/opt/intel/opencl/lib64 -I/opt/arrayfire/include -I/usr/local/cuda/include improved_Siddon_openCL_raw.cpp
        
        
        mex -largeArrayDims -lafopencl -lOpenCL -L/opt/arrayfire/lib -L/usr/local/cuda/lib64 -L/opt/intel/opencl/lib64 -I/opt/arrayfire/include -I/usr/local/cuda/include improved_Siddon_openCL_matrixfree_GPU.cpp
        mex -largeArrayDims -lafopencl -lOpenCL -L/opt/arrayfire/lib -L/usr/local/cuda/lib64 -L/opt/intel/opencl/lib64 -I/opt/arrayfire/include -I/usr/local/cuda/include improved_Siddon_openCL_matrixfree_GPU_raw.cpp
        
    end
end