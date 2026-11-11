

def loadProSpectaData(options):
    from pydicom import dcmread
    import numpy as np
    dcm = dcmread(options.fpath)

    # FOV alignment: SPECT image grid (world limits) used when resampling a CT attenuation image onto the SPECT image
    dz = float(dcm.PixelSpacing[0])
    Nz = float(dcm.Rows)
    sZ = float(dcm.DetectorInformationSequence[0].ImagePositionPatient[2])
    cZ = sZ - 0.5 * dz * (Nz - 1)
    tableHeight = float(dcm.RotationInformationSequence[0].TableHeight)
    FOVa_x = float(np.asarray(options.FOVa_x).ravel()[0])
    FOVa_y = float(np.asarray(options.FOVa_y).ravel()[0])
    axial_fov = float(np.asarray(options.axial_fov).ravel()[0])
    options.refSPECT = {'ImageSize': (int(np.asarray(options.Nx).ravel()[0]), int(np.asarray(options.Ny).ravel()[0]), int(np.asarray(options.Nz).ravel()[0])),
                        'XWorldLimits': (-0.5 * FOVa_x, 0.5 * FOVa_x),
                        'YWorldLimits': (-0.5 * FOVa_y - tableHeight, 0.5 * FOVa_y - tableHeight),
                        'ZWorldLimits': (-0.5 * axial_fov + cZ, 0.5 * axial_fov + cZ)}

    # Load projection images
    options.SinM = dcm.pixel_array
    options.SinM = options.SinM.transpose((2, 1, 0))
    
    # Number of rows in a projection image
    options.nRowsD = options.SinM.shape[0]

    # Number of columns in a projection image
    options.nColsD = options.SinM.shape[1]

    # Number of projections
    options.nProjections = options.SinM.shape[2]

    # Number of detector heads
    options.nHeads = 2
    
    # Starting angles for the two detector panels
    startAngle1 = (float)(dcm.DetectorInformationSequence[0].StartAngle)
    startAngle2 = (float)(dcm.DetectorInformationSequence[1].StartAngle)

    # Increment between adjacent projection angles
    angleIncrement = (float)(dcm.RotationInformationSequence[0].AngularStep)

    # Rotation angles
    options.angles = np.concatenate((np.arange(startAngle2, startAngle2 + angleIncrement * (options.nProjections / options.nHeads), angleIncrement),np.arange(startAngle1, startAngle1 + angleIncrement * (options.nProjections / options.nHeads - 1), angleIncrement)))

    # Radial distance of the detector panel from the COR
    options.radiusPerProj = np.concatenate((dcm.DetectorInformationSequence[0].RadialPosition, dcm.DetectorInformationSequence[1].RadialPosition))
    
    return 0