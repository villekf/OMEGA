# -*- coding: utf-8 -*-
"""
This example contains an example of PyTorch and OMEGA interoperability 
using projection SPECT data. The dataset is Siemens Pro.specta
projection data available at DOI 10.5281/zenodo.17315440.

This version implements a deep-image-prior reconstruction with a 3-D U-net.
The OMEGA forward and backward projectors are used to optimize the network
parameters directly from the measured SPECT projections.

This example uses PyTorch thus requires either a CUDA or Metal compatible device.
"""

# %% Imports and run configuration
from pathlib import Path

import numpy as np
from omegatomo.projector import proj
import torch
import torch.nn as nn
import torch.nn.functional as F
from pymatreader import read_mat

# %% OMEGA data, scanner, and reconstruction configuration
options = proj.projectorClass()

# Path to .mat file
options.fpath = './jaszczak_spectct_projection_data.mat' 


# Set PyTorch backend
options.SPECT = True # Required for SPECT data
options.useTorch = True # Use PyTorch tensors for storing data
if sys.platform == 'darwin':
    options.useMetal = True
    device = torch.device("mps")
else:
    options.useCUDA = True
    options.useCuPy = True # Use CuPy, PyCUDA support is deprecated
    device = torch.device("cuda")


###########################################################################
###########################################################################
###########################################################################
############################### LOAD DATA #################################
###########################################################################
###########################################################################
###########################################################################

data = read_mat(options.fpath)
options.SinM = np.array(data['projection_data'])
options.angles = np.array(data['angular_position']).squeeze()
options.radiusPerProj = np.array(data['radial_position']).squeeze()
options.nRowsD = options.SinM.shape[0]
options.nColsD = options.SinM.shape[1]
options.nProjections = options.SinM.shape[2]
energy_window = data.get("energy_window", None)
pixel_spacing = np.array(data["pixel_spacing"]).squeeze()
detector_thickness = float(np.squeeze(data["detector_thickness"]))

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### Crystal thickness (mm)
options.cr_p = detector_thickness

### Crystal width (mm)
options.dPitchX = float(pixel_spacing[0])
options.dPitchY = float(pixel_spacing[1])

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Prospecta'
 
###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################
 
### Reconstructed image pixel count. U-net uses 64x64x64 as input.
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 64; # X-direction
options.Ny = 64; # Y-direction
options.Nz = 64; # Z-direction

### FOV size [mm]
# NOTE: Non-cubical voxels may not work
options.FOVa_x = options.dPitchX*64; # [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*64; # [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*64; # [mm], z-axis of FOV (axial)

### Flip the image?
options.flipImageX = False
options.flipImageY = False
options.flipImageZ = False


### How much is the image rotated in degrees?
# NOTE: The rotation is done in the detector space (before reconstruction).
# Positive values perform the rotation in counterclockwise direction
options.offangle = 0

# Axial offset; this centers the phantom
options.oOffsetZ = -16 * options.dPitchX

###########################################################################
###########################################################################
###########################################################################
############################## CORRECTIONS ################################
###########################################################################
###########################################################################
###########################################################################

######################### Attenuation correction ##########################
# Currently scaling and resampling is not supported for the attenuation map.
options.attenuation_correction = False


########################### Resolution recovery ##########################
### Collimator-detector response function (CDRF)
# For projector types 2 and 6 you can either input either:
# 1. the collimator parameters (default) for an analytic solution for round (and hexagonal) holes (this may be unoptimal),
# 2. the standard deviations for both transaxial and axial directions or
# 3. the (Gaussian) PSF filter
# 4. the shifts of each ray traced 

# NOTE: For projector type 1 the CDRF is determined by
# options.rayShiftsDetector and options.rayShiftsSource defined in option 
# 4. These can also be calculated automatically when collimator parameters
# (1.) are input.

# NOTE: With projector_type == 2 (orthogonal distance projector), only the
# collimator parameters below are used for CDR calculation i.e. the
# collimator hole is assumed to be a circle. Thus only 1. below is
# supported with projector_type == 2
#
# 1. The collimator parameters (projector types 1, 2 and 6)
# Collimator hole length (mm)
options.colL = float(np.squeeze(data["collimator_thickness"]))
# Collimator hole radius (mm)
options.colR = float(np.squeeze(data["collimator_hole_radius"]))
# Distance from collimator to the detector (mm)
options.colD = 0.0
# Intrinsic resolution (mm)
options.iR = float(np.squeeze(data["detector_intrinsic_resolution"]))
# Focal distance (XY)
options.colFxy = np.inf
# Focal distance (Z)
options.colFz = np.inf

# 2. If you have the standard deviations for transaxial (XY) and axial (Z)
# directions, you can input them here instead of the above values The
# dimensions need to be options.nProjections x options.Nx. Only for
# projector type 6.
# options.sigmaZ = np.ones((options.nProjections, options.Nx), dtype=np.float32)
# options.sigmaXY = np.ones((options.nProjections, options.Nx), dtype=np.float32)

# 3. You can input the filter for the CDRF directly. This should be of the
# size filterSizeXY x filterSizeZ. Only for
# projector type 6.
# options.gFilter = np.ones((1, 1, options.Nx), dtype=np.float32)

# 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
# the shape of the collimator hole. The values of rayShiftsDetector and
# rayShiftsSource represent [shift1XY, shift1Z, shift2XY, ...] in mm. Size
# should be 2*n_rays_axial*n_rays_transaxial x nColsD x nRowsD x nHeads. If not input, values
# are calculated automatically.
options.n_rays_axial = 1
options.n_rays_transaxial = 1
# options.rayShiftsDetector = np.zeros((2*options.n_rays_axial*options.n_rays_transaxial, options.nColsD, options.nRowsD, options.nHeads));
# options.rayShiftsSource = np.zeros((2*options.n_rays_axial*options.n_rays_transaxial, options.nColsD, options.nRowsD, options.nHeads));
 
###########################################################################
###########################################################################
###########################################################################
############################# MISC PROPERTIES #############################
###########################################################################
###########################################################################
###########################################################################

### Name of current datafile/examination
# This is used to name the saved measurement data and also load it in
# future sessions.
options.name = 'spect_DIP_example'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1.  Maximum value of 3 is
# supported.
options.verbose = 1

###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = (Improved) Siddon ray-based projector
# 2 = Orthogonal distance ray tracing
# 6 = Rotation-based projector
# See the documentation on some details on the projectors:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
# NOTE: with rotation-based projector, the sinogram must be resized and
# resampled to match FOV XZ-plane size and resolution.
options.projector_type = 1

### Use images instead of buffers? For rotation-based projector this
# implies hardware texture interpolation, which typically has 8 bit 
# precision. With buffers, software interpolation with 32 bit floats is
# used.
options.useImages = False

###########################################################################
###########################################################################
####################### DEEP IMAGE PRIOR SETTINGS #########################
###########################################################################
###########################################################################

# %% Projector initialization
# DIP uses the complete measurement set in every network update.
options.Niter = 25
options.subsets = 1
options.subsetType = 8

# Intermediate saving
CHECKPOINT_EVERY = 0  # Set to 0 to disable intermediate checkpoints.
OUTPUT_DIR = "spect_DIP_output"

# %% Reconstruction helpers and U-net definition
class ConvBlock3D(nn.Module):
    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        negative_slope: float = 0.2,
    ):
        super().__init__()

        self.block = nn.Sequential(
            nn.Conv3d(
                in_channels,
                out_channels,
                kernel_size=3,
                padding=1,
                bias=False,
            ),
            nn.BatchNorm3d(
                out_channels,
                affine=True,
                track_running_stats=False,
            ),
            nn.LeakyReLU(negative_slope, inplace=False),
            nn.Conv3d(
                out_channels,
                out_channels,
                kernel_size=3,
                padding=1,
                bias=False,
            ),
            nn.BatchNorm3d(
                out_channels,
                affine=True,
                track_running_stats=False,
            ),
            nn.LeakyReLU(negative_slope, inplace=False),
        )

    def forward(self, x):
        return self.block(x)

# Network architecture inspired by 10.1088/1361-6560/ace49c
class UNet(nn.Module):
    def __init__(self):
        super().__init__()

        c1, c2, c3, c4 = 32, 64, 128, 256

        self.enc1 = ConvBlock3D(1, c1)
        self.down1 = nn.Conv3d(c1, c2, kernel_size=3, stride=2, padding=1)

        self.enc2 = ConvBlock3D(c2, c2)
        self.down2 = nn.Conv3d(c2, c3, kernel_size=3, stride=2, padding=1)

        self.enc3 = ConvBlock3D(c3, c3)
        self.down3 = nn.Conv3d(c3, c4, kernel_size=3, stride=2, padding=1)

        self.bottleneck = ConvBlock3D(c4, c4)

        self.up3 = nn.ConvTranspose3d(c4, c3, kernel_size=2, stride=2)
        self.dec3 = ConvBlock3D(c3 + c3, c3)

        self.up2 = nn.ConvTranspose3d(c3, c2, kernel_size=2, stride=2)
        self.dec2 = ConvBlock3D(c2 + c2, c2)

        self.up1 = nn.ConvTranspose3d(c2, c1, kernel_size=2, stride=2)
        self.dec1 = ConvBlock3D(c1 + c1, c1)

        self.output = nn.Conv3d(c1, 1, kernel_size=1)

    def forward(self, z):
        e1 = self.enc1(z)
        e2 = self.enc2(self.down1(e1))
        e3 = self.enc3(self.down2(e2))
        b = self.bottleneck(self.down3(e3))

        x = self.up3(b)
        x = self.dec3(torch.cat((x, e3), dim=1))

        x = self.up2(x)
        x = self.dec2(torch.cat((x, e2), dim=1))

        x = self.up1(x)
        x = self.dec1(torch.cat((x, e1), dim=1))

        return F.relu(self.output(x))

# Minimized objective functions to select from
def objGaussianLS(A, fp, y):
    residual = fp - y
    objective = torch.sum(residual.square())
    gradient_x_vector = 2.0 * (A.T() * residual)
    return objective, gradient_x_vector

def objPoissonLogLikelihood(A, fp, y):
    objective = torch.sum(fp - y * torch.log(fp))
    gradient_x_vector = A.T() * (1.0 - y / fp)
    return objective, gradient_x_vector

def objPoissonKL(A, fp, y):
    objective_terms = fp - y
    positive = y > 0
    objective_terms[positive] += (y[positive] * (torch.log(y[positive]) - torch.log(fp[positive])))
    objective = torch.sum(objective_terms)
    gradient_x_vector = A.T() * (1.0 - y / fp)
    return objective, gradient_x_vector


# %% Reconstruction
options.addProjector() # Initialize projector
options.initProj()

measured = np.asarray(options.SinM, dtype=np.float32).ravel(order="F")
y = torch.as_tensor(measured, dtype=torch.float32, device=device)

# NCDHW = (1,1,z,y,x).
fixed_noise = torch.randn((1, 1, int(options.Nz[0]), int(options.Ny[0]), int(options.Nx[0])), dtype=torch.float32, device=device)

model = UNet().to(device)

optimizer = torch.optim.LBFGS(
    model.parameters(),
    lr=1.0,
    max_iter=80,
    max_eval=100,
    history_size=20,
    tolerance_grad=1e-7,
    tolerance_change=1e-9,
    line_search_fn="strong_wolfe",
)

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
loss_history = []

model.train()
for iteration in range(1, options.Niter + 1):
    def closure():
        # LBFGS can call this function more than once, so every evaluation must start from clean gradients and recompute the current model.
        optimizer.zero_grad(set_to_none=True)

        x_ncdhw = model(fixed_noise)
        x_vector = x_ncdhw.reshape(-1)

        with torch.no_grad():
            options.subset = 0
            dip_fp = options * x_vector.detach() + options.epps # Forward projection of U-net output. Add scatter here if used.
            
            # LS, Gaussian noise:
            #objective, gradient_x_vector = objGaussianLS(options, dip_fp, y)
            
            # Negative poisson log-likelihood
            #objective, gradient_x_vector = objPoissonLogLikelihood(options, dip_fp, y)
            
            # KL divergence
            objective, gradient_x_vector = objPoissonKL(options, dip_fp, y)
            
            gradient_x = gradient_x_vector.reshape_as(x_ncdhw).to(
                device=x_ncdhw.device,
                dtype=x_ncdhw.dtype,
            )

        x_ncdhw.backward(gradient=gradient_x)
        return objective

    objective = optimizer.step(closure)
    objective_value = float(objective.detach().cpu())
    loss_history.append(objective_value)
    print(
        f"DIP iteration {iteration}/{options.Niter}: "
        f"objective={objective_value:.7e}"
    )

    if CHECKPOINT_EVERY and iteration % CHECKPOINT_EVERY == 0:
        torch.save({
            "iteration": iteration,
            "shape_xyz": (int(options.Nx[0]), int(options.Ny[0]), int(options.Nz[0])),
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "fixed_noise": fixed_noise.detach().cpu(),
            "loss_history": list(loss_history),
        }, output_dir / "spect_dip_checkpoint.pt",)

model.eval()
with torch.no_grad():
    f_DIP = model(fixed_noise)

# Convert (z,y,x) C-order storage to (x,y,z) NumPy volume.
f_DIP = f_DIP[0, 0].detach().cpu().numpy()
f_DIP = np.transpose(f_DIP, (2, 1, 0))

np.save(output_dir / "spect_dip_reconstruction.npy", f_DIP)
np.save(output_dir / "spect_dip_loss_history.npy", np.asarray(loss_history))
torch.save({
    "iteration": iteration,
    "shape_xyz": (int(options.Nx[0]), int(options.Ny[0]), int(options.Nz[0])),
    "model_state_dict": model.state_dict(),
    "optimizer_state_dict": optimizer.state_dict(),
    "fixed_noise": fixed_noise.detach().cpu(),
    "loss_history": list(loss_history),
}, output_dir / "spect_dip_checkpoint.pt",)

# %% Plot
from omegatomo.util.volume3Dviewer import volume3Dviewer
volume3Dviewer(f_DIP)
