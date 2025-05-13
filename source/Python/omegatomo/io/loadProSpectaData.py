from pydicom import dcmread
import numpy as np

def loadProSpectaData(options):
    dcm = dcmread(options.fpath)

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