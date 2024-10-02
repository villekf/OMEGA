# -*- coding: utf-8 -*-
import subprocess
import os
import argparse

class CommandLine:
    def __init__(self):
        parser = argparse.ArgumentParser(description = "Additional library paths")
        parser.add_argument("-A", "--AF", help = "Example: -A /path/to/ArrayFire", required = False, default = "")
        parser.add_argument("-O", "--OpenCL", help = "Example: -O /path/to/OpenCL", required = False, default = "")
        parser.add_argument("-R", "--Root", help = "Example: -R /path/to/ROOT", required = False, default = "")
        arg = parser.parse_args()
        openclpath = arg.OpenCL
        rpath = arg.Root
        afpath = arg.AF

        if os.name == 'nt':
            homedir = os.path.expanduser('~')
            if len(openclpath) > 0:
                opencllib = os.path.join(openclpath, 'lib', 'x64')
                openclpath = os.path.join(openclpath, 'include')
            elif 'CUDA_PATH' in os.environ:
                openclpath = os.path.join(os.environ['CUDA_PATH'], 'include')
                opencllib = os.path.join(os.environ['CUDA_PATH'] , 'lib', 'x64')
            elif 'OCL_ROOT' in os.environ: 
                openclpath = os.path.join(os.environ['OCL_ROOT'], 'include')
                opencllib = os.path.join(os.environ['OCL_ROOT'], 'lib', 'x86_64')
            elif 'INTELOCLSDKROOT' in os.environ:
                openclpath = os.path.join(os.environ['OCL_ROOT'], 'include')
                opencllib = os.path.join(os.environ['OCL_ROOT'], 'lib', 'x64')
            elif os.path.exists('C:\\Program Files (x86)\\OCL_SDK_Light'):
                openclpath = 'C:\\Program Files (x86)\\OCL_SDK_Light\\include'
                opencllib = 'C:\\Program Files (x86)\\OCL_SDK_Light\\lib\\x86_64'
            elif os.path.exists(homedir + '\\Documents\\OCL_SDK_Light'):
                openclpath = os.path.join(homedir, '\\Documents\\OCL_SDK_Light\\include')
                opencllib = os.path.join(homedir ,'\\Documents\\OCL_SDK_Light\\lib\\x86_64')
            elif os.path.exists(homedir + '\\Downloads\\OCL_SDK_Light'):
                openclpath = os.path.join(homedir, '\\Downloads\\OCL_SDK_Light\\include')
                opencllib = os.path.join(homedir, '\\Downloads\\OCL_SDK_Light\\lib\\x86_64')
            else:
                print('OpenCL not found! Install OpenCL library and headers, for example https://github.com/GPUOpen-LibrariesAndSDKs/OCL-SDK/releases and make sure that OCL_ROOT is set.')
            compiler = "cl.exe"
            if len(afpath) > 0:
                aflib = os.path.join(afpath, 'lib')
                afpath = os.path.join(afpath, 'include')
            elif 'AF_PATH' in os.environ:
                afpath = os.path.join(os.environ['AF_PATH'], "include")
                aflib = os.path.join(os.environ['AF_PATH'], "lib")
            elif os.path.exists(homedir + '\\Documents\\ArrayFire\\v3'):
                afpath = homedir + '\\Documents\\ArrayFire\\v3\\include'
                aflib = homedir + '\\Documents\\ArrayFire\\v3\\lib\\x86_64'
            elif os.path.exists(homedir + '\\Downloads\\ArrayFire\\v3'):
                afpath = homedir + '\\Downloads\\ArrayFire\\v3\\include'
                aflib = homedir + '\\Downloads\\ArrayFire\\v3\\lib\\x86_64'
            elif os.path.exists(homedir + '\\Documents\\ArrayFire'):
                afpath = homedir + '\\Documents\\ArrayFire\\include'
                aflib = homedir + '\\Documents\\ArrayFire\\lib\\x86_64'
            elif os.path.exists(homedir + '\\Downloads\\ArrayFire'):
                afpath = homedir + '\\Downloads\\ArrayFire\\include'
                aflib = homedir + '\\Downloads\\ArrayFire\\lib\\x86_64'
            if len(afpath) == 0:
                raise ValueError('ArrayFire not found! Please install ArrayFire and make sure that AF_PATH is set')
            sdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'cpp'))
            if len(openclpath) > 0:
                try:
                    options = '/O2 /DOPENCL /DAF /DAF_RELEASE /DAF_OPENCL /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + afpath + '"' + ' /I"' + openclpath + '"'
                    libs = '"OpenCL.lib" "afopencl.lib"'
                    link = '/link ' + libs + ' /LIBPATH:"' + aflib + '" /LIBPATH:"' + opencllib + '" /MACHINE:X64 /OUT:"OpenCL_matrixfree_lib.dll"'
                    files = '"' + sdir + '\\omega_maincpp.cpp"'
                    compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                    # compile_command = [compiler, options, files, link]
                    result = subprocess.run(compile_command, check=True)
                    if result.stderr is None:
                        options = '/O2 /DOPENCL /DAF /DMTYPE /DAF_RELEASE /DAF_OPENCL /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + afpath + '"' + ' /I"' + openclpath + '"'
                        link = '/link ' + libs + ' /LIBPATH:"' + aflib + '" /LIBPATH:"' + opencllib + '" /MACHINE:X64 /OUT:"OpenCL_matrixfree_uint16_lib.dll"'
                        compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                        result = subprocess.run(compile_command, check=True)
                        print('OpenCL support compiled successfully!')
                except Exception:
                    print("OpenCL build failed")
            cudapath = ''
            if 'CUDA_PATH' in os.environ:
                cudapath = os.environ['CUDA_PATH'] + '\\include'
                cudalib = os.environ['CUDA_PATH'] + '\\lib\\x64'
            if len(cudapath) > 0:
                try:
                    options = '/O2 /DCUDA /DAF /DAF_RELEASE /DAF_CUDA /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + afpath + '"' + ' /I"' + cudapath + '"'
                    libs = '"cuda.lib" "nvrtc.lib" "afcuda.lib"'
                    link = '/link ' + libs + ' /LIBPATH:"' + aflib + '" /LIBPATH:"' + cudalib + '" /MACHINE:X64 /OUT:"CUDA_matrixfree_lib.dll"'
                    files = '"' + sdir + '\\omega_maincpp.cpp"'
                    compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                    # compile_command = [compiler, options, files, link]
                    result = subprocess.run(compile_command, check=True)
                    if result.stderr is None:
                        options = '/O2 /DCUDA /DAF /DMTYPE /DAF_RELEASE /DAF_CUDA /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + afpath + '"' + ' /I"' + cudapath + '"'
                        link = '/link ' + libs + ' /LIBPATH:"' + aflib + '" /LIBPATH:"' + cudalib + '" /MACHINE:X64 /OUT:"CUDA_matrixfree_uint16_lib.dll"'
                        compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                        result = subprocess.run(compile_command, check=True)
                        print('CUDA support compiled successfully!')
                except Exception:
                    print("CUDA build failed")
            else:
                print('CUDA not found. No CUDA code compiled!')
            
            if len(rpath) > 0:
                rlib = os.path.join(rpath, 'lib')
                rpath = os.path.join(rpath, 'include')
            elif os.path.exists('C:/root'):
                rlib = os.path.join('C:/root', 'lib')
                rpath = os.path.join('C:/root', 'include')
            elif os.path.exists('C:/Program Files/root'):
                rlib = os.path.join('C:/Program Files/root', 'lib')
                rpath = os.path.join('C:/Program Files/root', 'include')
            else:
                for yy in range(6,13):
                    if os.path.exists(rpath):
                        break
                    for xx in range(0,60,2):
                        if os.path.exists(rpath):
                            break
                        for zz in range(0,60):
                            if xx < 10:
                                xc = '0' + str(xx)
                            else:
                                xc = str(xx)
                            if zz < 10:
                                zc = '0' + str(zz)
                            else:
                                zc = str(zz)
                            rpath = 'C:/root_v' + str(yy) + '.' + xc + '.' + zc
                            if os.path.exists(rpath):
                                break
                            rpath = 'C:/Program Files/root_v' + str(yy) + '.' + xc + '.' + zc
                            if os.path.exists(rpath):
                                break
                rlib = os.path.join(rpath, 'lib')
                rpath = os.path.join(rpath, 'include')
                    
                
                
            try:
                options = '/O2 /openmp /DCPU /DAF /DAF_RELEASE /DAF_CPU /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + afpath + '"'
                libs = '"afcpu.lib"'
                link = '/link ' + libs + ' /LIBPATH:"' + aflib + '" /MACHINE:X64 /OUT:"CPU_matrixfree_lib.dll"'
                files = '"' + sdir + '\\omega_maincpp.cpp"'
                compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                result = subprocess.run(compile_command, check=True)
                if result.stderr is None:
                    print('CPU support compiled successfully!')
            except Exception:
                print("CPU build failed")
                
            try:
                options = '/O2 /D__x86_64 /LD /EHsc /I"' + sdir + '"'
                link = '/link ' + '/MACHINE:X64 /OUT:"createSinogram.dll"'
                files = '"' + sdir + '\\createSinogram.cpp"'
                compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                subprocess.run(compile_command, check=True)
            except Exception:
                print("Build failed")
            
            try:
                options = '/O2 /D__x86_64 /LD /EHsc /I"' + sdir + '"'
                link = '/link ' + '/MACHINE:X64 /OUT:"inveon.dll"'
                files = '"' + sdir + '\\inveonMain.cpp"'
                compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                subprocess.run(compile_command, check=True)
                print('Inveon list-mode support compiled successfully!')
            except Exception:
                print("Build failed")
                
            try:
                options = '/std:c++17 /O2 /D__x86_64 /LD /EHsc /I"' + sdir + '" /I"' + rpath + '"'
                libs = '"libCore.lib" "libRIO.lib" "libTree.lib"'
                link = '/link ' + libs + ' /LIBPATH:"' + rlib + '" /MACHINE:X64 /OUT:"libRoot.dll"'
                files = '"' + sdir + '\\loadRootData.cpp"'
                compile_command = compiler + ' ' + options + ' ' + files + ' ' + link
                result = subprocess.run(compile_command, check=True)
                print('ROOT support compiled successfully!')
            except Exception:
                print("ROOT support build failed")
        
        else:
            compiler = 'g++'
            if len(openclpath) > 0:
                opencllib = os.path.join(openclpath, '/lib')
                openclpath = os.path.join(openclpath, '/include')
            elif os.path.exists('/opt/intel/opencl/'):
                openclpath = '/opt/intel/opencl/include'
                opencllib = '/opt/intel/opencl/lib/x64'
            elif os.path.exists('/usr/local/cuda/targets/x86_64-linux'):
                openclpath = '/usr/local/cuda/targets/x86_64-linux/include'
                opencllib = '/usr/local/cuda/lib64'
            elif os.path.exists('/usr/lib/x86_64-linux-gnu/'):
                openclpath = '/usr/include'
                opencllib = '/usr/lib/x86_64-linux-gnu/'
            # if len(openclpath) == 0:
            #     raise ValueError('OpenCL not found! Please install OpenCL headers and library')
            if len(afpath) > 0:
                aflib = os.path.join(afpath, 'lib64')
                afpath = os.path.join(afpath, 'include')
            elif 'AF_PATH' in os.environ:
                afpath = os.environ['AF_PATH'] + "/include"
                aflib = os.environ['AF_PATH'] + "/lib64"
            elif os.path.exists('/opt/ArrayFire/'):
                afpath = '/opt/ArrayFire/include'
                aflib = '/opt/ArrayFire/lib64'
            # if len(afpath) == 0:
            #     raise ValueError('ArrayFire not found! Please install ArrayFire and input install directory with -A /path/to/arrayfire')
                
            sdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'cpp'))
            # options = ['-shared', '-fPIC', '-DOPENCL', '-DAF', '-DAF_OPENCL', '-I' + sdir, '-I' + afpath, '-I' + openclpath]
            # options = '-shared -fPIC -DOPENCL -DAF -DAF_RELEASE -DAF_OPENCL -I' + sdir + ' -I' + afpath + ' -I' + openclpath
            lib1 = '-L' + opencllib
            lib2 = '-L' + aflib
            lib3 = '-lOpenCL'
            lib4 = '-lafopencl'
            link1 = '-Wl,--export-dynamic'
            link = '-oOpenCL_matrixfree_lib.so'
            files = '' + sdir + '/omega_maincpp.cpp'
            # compile_command = options + ' ' + link + ' ' + files + ' ' + libs
            # compile_command = [compiler, options, files, link]
            # print(compile_command)
            try:
                result = subprocess.run([compiler, '-shared', '-fPIC', '-DOPENCL', '-DAF', '-DAF_OPENCL', '-I' + sdir, '-I' + afpath, '-I' + openclpath, 
                                         link1, link, files, lib1, lib2, lib3, lib4], check=True)
                if result.stderr is None:
                    link = '-oOpenCL_matrixfree_uint16_lib.so'
                    result = subprocess.run([compiler, '-shared', '-fPIC', '-DOPENCL', '-DAF', '-DAF_OPENCL', '-DMTYPE', '-I' + sdir, '-I' + afpath, '-I' + openclpath, 
                                             link1, link, files, lib1, lib2, lib3, lib4], check=True)
                    print('OpenCL support compiled successfully!')
            except Exception:
                print("OpenCL build failed")
            cudapath = ''
            if os.path.exists('/usr/local/cuda/targets/x86_64-linux'):
                cudapath = '/usr/local/cuda/targets/x86_64-linux/include'
                cudalib = '/usr/local/cuda/lib64'
            elif os.path.exists('/usr/lib/x86_64-linux-gnu/'):
                cudapath = '/usr/include'
                cudalib = '/usr/lib/x86_64-linux'
            # if len(cudapath) == 0:
            #     raise ValueError('CUDA not found. CUDA version not compiled.')
            lib1 = '-L' + cudalib
            lib3 = '-lcuda'
            lib4 = '-lafcuda'
            lib5 = '-lnvrtc'
            link = '-oCUDA_matrixfree_lib.so'
            files = '' + sdir + '/omega_maincpp.cpp'
            try:
                result = subprocess.run([compiler, '-shared', '-fPIC', '-DCUDA', '-DAF', '-DAF_CUDA', '-I' + sdir, '-I' + afpath, '-I' + cudapath, 
                                         link1, link, files, lib1, lib2, lib3, lib4, lib5], check=True)
                if result.stderr is None:
                    link = '-oCUDA_matrixfree_uint16_lib.so'
                    result = subprocess.run([compiler, '-shared', '-fPIC', '-DCUDA', '-DAF', '-DAF_CUDA', '-DMTYPE', '-I' + sdir, '-I' + afpath, '-I' + cudapath, 
                                             link1, link, files, lib1, lib2, lib3, lib4, lib5], check=True)
                    print('CUDA support compiled successfully!')
            except Exception:
                print("CUDA build failed")
                
            lib4 = '-lafcpu'
            # lib3 = '-liomp5'
            link = '-oCPU_matrixfree_lib.so'
            files = '' + sdir + '/omega_maincpp.cpp'
            try:
                result = subprocess.run([compiler, '-shared','-fopenmp', '-fPIC', '-DCPU', '-DAF', '-DAF_CPU', '-I' + sdir, '-I' + afpath,
                                         link1, link, files, lib2, lib4], check=True)
                if result.stderr is None:
                    print('CPU support compiled successfully!')
            except Exception:
                print("CPU build failed")
                
            link = '-ocreateSinogram.so'
            files = '' + sdir + '/createSinogram.cpp'
            try:
                result = subprocess.run([compiler, '-shared', '-fPIC', '-I' + sdir, link, files], check=True)
            except Exception:
                print("Build failed")
                
            link = '-oinveon.so'
            files = '' + sdir + '/inveonMain.cpp'
            try:
                result = subprocess.run([compiler, '-shared', '-fPIC', '-I' + sdir, link, files], check=True)
                print('Inveon list-mode support compiled successfully!')
            except Exception:
                print("Build failed")
            
            if len(rpath) > 0:
                rlib = os.path.join(rpath, 'lib')
                rpath = os.path.join(rpath, 'include')
            else:
                rpath = result = os.popen('root-config --incdir').read()
                rpath = rpath[:-1]
                rlib = result = os.popen('root-config --libdir').read()
                rlib = rlib[:-1]
                
            lib1 = '-L' + rlib
            lib2 = '-lCore'
            lib3 = '-lRIO'
            lib4 = '-lTree'
            lib5 = '-ldl'
            link = '-olibRoot.so'
            files = '' + sdir + '/loadRootData.cpp'
            try:
                result = subprocess.run([compiler, '-shared', '-fPIC', '-pthread', '-std=c++17', '-m64', '-I' + rpath, link, files, lib1, lib2, lib3, lib4, lib5], check=True)
                print('ROOT support compiled successfully!')
            except Exception:
                print("ROOT support build failed")
            
if __name__ == '__main__':
    CommandLine()