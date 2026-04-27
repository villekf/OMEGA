# -*- coding: utf-8 -*-
def checkCUDA(deviceNum = 0):
    try:
        import arrayfire as af
        try:
            if af.get_active_backend() != 'opencl':
                af.set_backend('opencl')
            info = af.device.info_str()
            loc = info.find('[' + str(deviceNum) + ']')
            if loc == -1:
                loc = info.find('-' + str(deviceNum) + '-')
            info = info[loc:]
            loc2 = info.find("\n")
            cuda = info[:loc2].find('NVIDIA')
        except:
            af.set_backend(af.BackendType.opencl)
            info = af.info_string()
            loc = info.find('[' + str(deviceNum) + ']')
            if loc == -1:
                loc = info.find('-' + str(deviceNum) + '-')
            info = info[loc:]
            loc2 = info.find("\n")
            cuda = info[:loc2].find('NVIDIA')
        if cuda > 0:
            return True
        else:
            return False
    except ModuleNotFoundError:
        print('ArrayFire not found. Unable to determine whether the GPU is CUDA capable or not!')
        return False