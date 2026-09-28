# -*- coding: utf-8 -*-
"""
Windows DLL search path helper.

Since Python 3.8, ctypes.CDLL on Windows resolves the dependencies of a
library with LOAD_LIBRARY_SEARCH_DEFAULT_DIRS | LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR,
which does not include PATH. The ArrayFire libraries (afopencl, afcpu, afcuda)
and NVRTC/HIP/ROCm (amdhip64, hiprtc, hipblas, hipfft, hipsolver, hipsparse)
that the OMEGA libraries are linked against are normally only found through
PATH, and thus have to be registered explicitly. The same applies to ROOT
(libCore), which libRoot.dll imports; ROOT has no dedicated path environment
variable in OMEGA, so it is located by searching PATH instead.
"""

import os

# Directories that have already been registered
_addedDirs = set()

def _findDirOnPath(fileName):
    """
    Returns the first directory on PATH that contains fileName, or an empty string.

    Used for libraries such as ROOT that have no dedicated path environment variable
    and are found through PATH only.
    """
    for pathDir in os.environ.get('PATH', '').split(os.pathsep):
        if len(pathDir) == 0:
            continue
        try:
            if os.path.isfile(os.path.join(pathDir, fileName)):
                return pathDir
        except OSError:
            continue
    return ''

def addDLLDirectories():
    """
    Adds the ArrayFire, HIP/ROCm, CUDA and ROOT runtime directories to the DLL
    search path.

    Does nothing on non-Windows platforms. The search path entries are
    intentionally never removed, i.e. they remain valid for the lifetime of
    the process. Calling this function several times is safe.
    """
    if os.name != 'nt' or not hasattr(os, 'add_dll_directory'):
        return
    dllDirs = []
    if 'AF_PATH' in os.environ:
        dllDirs.append(os.path.join(os.environ['AF_PATH'], 'lib'))
    hipRoot = ''
    if 'HIP_PATH' in os.environ:
        hipRoot = os.environ['HIP_PATH']
    elif 'ROCM_PATH' in os.environ:
        hipRoot = os.environ['ROCM_PATH']
    if len(hipRoot) > 0:
        hipRoot = os.path.normpath(hipRoot)
        dllDirs.append(os.path.join(hipRoot, 'bin'))
        dllDirs.append(os.path.join(hipRoot, 'lib'))
        # The ROCm Python wheel (TheRock) splits its runtime across sibling
        # directories named _rocm_sdk_*: _rocm_sdk_core holds amdhip64/hiprtc,
        # while the hipblas/hipfft/hipsolver/hipsparse DLLs that afcuda.dll
        # depends on live in _rocm_sdk_libraries\bin. Register all of them.
        parentDir = os.path.dirname(hipRoot)
        try:
            for entry in os.listdir(parentDir):
                if entry.startswith('_rocm_sdk_'):
                    dllDirs.append(os.path.join(parentDir, entry, 'bin'))
        except OSError:
            pass
    if 'CUDA_PATH' in os.environ:
        dllDirs.append(os.path.join(os.environ['CUDA_PATH'], 'bin'))
    # ROOT has no dedicated path environment variable in OMEGA. ROOTSYS is set only when
    # thisroot.bat has been sourced, so fall back to searching PATH for the directory that
    # actually contains libCore.dll, which is what libRoot.dll imports.
    rootDir = ''
    if 'ROOTSYS' in os.environ:
        rootDir = os.path.join(os.environ['ROOTSYS'], 'bin')
    if not os.path.isdir(rootDir):
        rootDir = _findDirOnPath('libCore.dll')
    if not os.path.isdir(rootDir):
        for rootTest in ['C:\\root\\bin', 'C:\\Program Files\\root\\bin']:
            if os.path.isdir(rootTest):
                rootDir = rootTest
                break
    if len(rootDir) > 0:
        dllDirs.append(rootDir)
    for dllDir in dllDirs:
        if dllDir in _addedDirs or not os.path.isdir(dllDir):
            continue
        try:
            os.add_dll_directory(dllDir)
            _addedDirs.add(dllDir)
        except OSError:
            pass
