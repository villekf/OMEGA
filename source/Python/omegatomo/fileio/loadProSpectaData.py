

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
    # The frames are stored as (frame, row, column). The image is arranged to (column, reversed row, frame) (the same as
    # rot90(squeeze(dicomread(fpath)), 3) in MATLAB): the DICOM rows run from the head to the feet, the reversed rows run along the
    # axial coordinate of the SPECT image and the CT attenuation image
    options.SinM = np.asfortranarray(dcm.pixel_array.transpose((2, 1, 0))[:, ::-1, :])
    
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

    # Additional offangle for CT (registration of the SPECT and CT images, determined with a Jaszczak phantom)
    additionalOffAngleCT = -5.

    # Rotation angles. The frames of the file are ordered by the detector head (all the frames of head 1 followed by all the
    # frames of head 2) and the angles are in the opposite direction of the DICOM start angles (same as in MATLAB)
    nProjectionsPerHead = options.nProjections // options.nHeads
    options.angles = -np.concatenate((startAngle1 + angleIncrement * np.arange(nProjectionsPerHead),
                                      startAngle2 + angleIncrement * np.arange(nProjectionsPerHead))) + additionalOffAngleCT

    # Radial distance of the detector panel from the COR (the same frame order as the angles)
    options.radiusPerProj = np.concatenate((np.asarray(dcm.DetectorInformationSequence[0].RadialPosition, dtype=np.float64),
                                            np.asarray(dcm.DetectorInformationSequence[1].RadialPosition, dtype=np.float64)))
    
    return 0