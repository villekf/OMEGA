#pragma once
#if defined(_MSC_VER)
#ifndef DLL_FUNCTION 
#define DLL_FUNCTION __declspec(dllimport)
#endif
#else
#define DLL_FUNCTION  __attribute__((visibility("default")))
#endif
#include "reconstructionAF.h"

struct inputStruct {
    // Is raw data format used (optional)
    uint8_t use_raw_data = 0;
    // Is list-mode data used (optional)
    uint8_t listmode = 0;
    // Verbosity level (optional)
    int8_t verbose = 0;
    // Number of transaxial rays (optional)
    uint16_t n_rays_transaxial = 1;
    // Number of axial rays (optional)
    uint16_t n_rays_axial = 1;
    // Projector type (optional)
    uint32_t projector_type = 1;
    // If 1, use attenuation correction (PET and SPECT) (optional)
    uint32_t attenuation_correction = 0;
    // If 1, use additional correction (PET and SPECT) (optional)
    uint32_t additionalCorrection = 0;
    // If 1, use normalization correction (PET) (optional)
    uint32_t normalization_correction = 0;
    // If 1, use scatter and/or randoms correction (PET, SPECT and CT) (optional)
    uint32_t randoms_correction = 0;
    // Number of columns in the projection image/sinogram
    uint32_t nColsD;
    // Number of rows in the projection image/sinogram
    uint32_t nRowsD;
    // Number of angular samples in the sinogram
    uint32_t Nang;
    // Number of radial distances in the sinogram
    uint32_t Ndist;
    // Number of subsets (optional)
    uint32_t subsets = 1;
    // Number of detectors per ring (PET) (optional)
    uint32_t det_per_ring = 1;
    // Number of PET rings (PET) (optional)
    uint32_t rings = 1;
    // Original size of the image (matrix) if multi-resolution reconstruction is used (x-axis)
    uint32_t NxOrig;
    // Original size of the image (matrix) if multi-resolution reconstruction is used (y-axis)
    uint32_t NyOrig;
    // Original size of the image (matrix) if multi-resolution reconstruction is used (z-axis)
    uint32_t NzOrig;
    // For multi-resolution reconstruction, the size of the prior computation (equals the original image size) (x-axis)
    uint32_t NxPrior;
    // For multi-resolution reconstruction, the size of the prior computation (equals the original image size) (y-axis)
    uint32_t NyPrior;
    // For multi-resolution reconstruction, the size of the prior computation (equals the original image size) (z-axis)
    uint32_t NzPrior;
    // Number of iterations (optional)
    uint32_t Niter = 1;
    // Number of time steps (optional), for dynamic imaging
    uint32_t Nt = 1;
    // Subset selection type (optional)
    uint32_t subsetType = 1;
    // Number of multi-volumes (optional)
    uint32_t nMultiVolumes = 0;
    // Number of crystal layers (optional)
    uint32_t nLayers = 1;
    // If 1, the PDHG primal and dual step sizes are determined adaptively (optional)
    uint32_t PDAdaptiveType = 0;
    // Number of power method iterations for PDHG primal and dual steps sizes (optional)
    uint32_t powerIterations = 20;
    // Number of deblurring iterations if deblurring is selected (optional)
    uint32_t deblur_iterations = 10;
    // The first iteration to take into account the gradient-based preconditioner (optional)
    uint32_t gradInitIter = 0;
    // Same as above, but the last iteration (the gradient is assumed constant after this) (optional)
    uint32_t gradLastIter = Niter - 1;
    // The number of filtering steps (optional)
    uint32_t filteringIterations = 20;
    // Mean type for weighted mean prior (optional)
    uint32_t mean_type = 1;
    // Neighborhood size for MRP, NLM, GGMRF, Huber, quadratic, weighted mean and hyperbolic prior (optional)
    // x-axis
    uint32_t Ndx = 1;
    // y-axis
    uint32_t Ndy = 1;
    // z-axis
    uint32_t Ndz = 1;
    // Patch size for NLM (optional)
    // x-axis
    uint32_t Nlx = 1;
    // y-axis
    uint32_t Nly = 1;
    // z-axis
    uint32_t Nlz = 1;
    // Size of the Gaussian PSF (optional)
    // x-axis
    uint32_t g_dim_x = 0;
    // y-axis
    uint32_t g_dim_y = 0;
    // z-axis
    uint32_t g_dim_z = 0;
    // Number of iterations with anisotropic diffusion smoothing (optional)
    uint32_t NiterAD = 1;
    // The index for the center voxel for some priors, such as quadratic (optional)
    uint32_t inffi = 0;
    uint32_t Nf;
    // Selected device (optional)
    uint32_t deviceNum = 0;
    // selected platform (unused!)
    uint32_t platform = 0;
    // Numerical derivative type, i.e. 0 is forward difference (unused!)
    uint32_t derivativeType = 0;
    // Type of TV prior (optional)
    uint32_t TVtype = 1;
    uint32_t FluxType;
    uint32_t DiffusionType;
    uint32_t POCS_NgradIter;
    // Number of projections
    int64_t nProjections = 1;
    // Number of TOF bins
    int64_t TOF_bins = 1;
    // Small value for TV type 1 (optional)
    float tau = 0.f;
    // Tube radius when using volume of intersection projector
    float tube_radius;
    // Small constant
    float epps = 1e-5f;
    // TOF Gaussian standard deviation (optional)
    float sigma_x;
    // Tube width for the orthogonal distance-based ray tracer
    float tube_width_z;
    // Unused
    float tube_width_xy;
    // volume of intersection projector
    float bmin;
    float bmax;
    float Vmax;
    // Global correction factor (scalar)
    float global_factor = 1.f;
    // Interpolation length with interpolation-based ray-tracer
    float dL;
    // Flat value (CT only)
    float flat;
    // Upper bound for MBSREM and MRAMLA
    float U = 1000.f;
    // Acceleration factor for ACOSEM
    float h_ACOSEM;
    // Detector pitch (x/y-axis)
    float dPitchX;
    // Detector pitch (z-axis)
    float dPitchY;
    // Detector pitch (x/y-axis)
    float cr_p;
    // Detector pitch (z-axis)
    float cr_pz;
    // Gaussian standard devation for NLM
    float NLMsigma;
    // Sum of the weights for quadratic, Huber and hyperbolic priors
    float w_sum;
    // Anisotropic diffusion smoothing parameters
    float KAD;
    float TimeStepAD;
    // RDP scaling value
    float RDP_gamma;
    // Huber bound
    float huber_delta;
    // Bounds for the gradient-based preconditioner
    float gradV1;
    float gradV2;
    // Regularization parameters for TGV
    float alpha0TGV;
    float alpha1TGV;
    // Adjustable parameter for GGMRF
    float GGMRF_p;
    float GGMRF_q;
    float GGMRF_c;
    // Regularization parameter
    float beta = 0.f;
    // Unused
    float T;
    // Dimension sizes for BDD
    float dSizeXBP;
    float dSizeZBP;
    // TV smoothing value
    float TVsmoothing;
    // Anatomical weight
    float C;
    // Lange function weight
    float SATVPhi;
    // Smoothing value for APLS
    float eta;
    // Smoothing value for APLS
    float APLSsmoothing;
    // Hyperbolic prior weight
    float hyperbolicDelta;
    // Source to detector distance for FDK weights
    float sourceToCRot = 0.f;
    float POCS_alpha = 0.2f;
    float POCS_rMax = 0.95f;
    float POCS_alphaRed = 0.95f;
    float POCSepps = 1e-4f;
    // Is PSF smoothing used
    bool use_psf = false;
    // Is TOF used
    bool TOF = false;
    // Is the detector panel tilted
    bool pitch = false;
    // Is SPECT data used
    bool SPECT = false;
    // Is full sinogram-based PET data used
    bool PET = false;
    // Is CT-data used
    bool CT = false;
    // Use high-dimensional reconstruction
    bool largeDim = false;
    // Load TOF data
    bool loadTOF = true;
    // If true, stores the primal-dual gap values
    bool storeResidual = false;
    // Remove the DC-component from the BDD forward projection
    bool meanFP = false;
    // Remove the DC-component from the BDD backward projection
    bool meanBP = false;
    // Use forward projection mask
    bool useMaskFP = false;
    // Use backward projection mask
    bool useMaskBP = false;
    // Use transaxial orthogonal distance-based ray tracer
    bool orthTransaxial = false;
    // Use axial orthogonal distance-based ray tracer
    bool orthAxial = false;
    // Enforce positivity with non-Poisson algorithms
    bool enforcePositivity = false;
    // Use multi-resolution reconstruction
    bool useMultiResolutionVolumes = false;
    // Save all iterations to host
    bool save_iter = false;
    // Use deblurring
    bool deblurring = false;
    // Use FMAD
    bool useMAD = false;
    // Use OpenCL images/CUDA textures
    bool useImages = false;
    // Use extended FOV
    bool useEFOV = false;
    // Use CT-based attenuation images
    bool CTAttenuation = true;
    // Use offset correction
    bool offsetCorrection = false;
    // Scale the relaxation parameters
    bool relaxationScaling = false;
    // Compute the relaxation parameters
    bool computeRelaxationParameters = false;
    // Store all forward projections
    bool storeFP = false;
    // Use 2D TGV
    bool use2DTGV = false;
    // Use MRP without normalization
    bool med_no_norm = false;
    // Use NLM with MRP type smoothing
    bool NLM_MRP = false;
    // Use NLTV
    bool NLTV = false;
    // Use NLRDP
    bool NLRD = false;
    // Use NL-Lange
    bool NLLange = false;
    // Use NLGGMRF
    bool NLGGMRF = false;
    // Use reference image for NLM
    bool NLM_use_anatomical = false;
    // Use anatomical weighting for TV
    bool TV_use_anatomical = false;
    // Include neighboring corners with RDP
    bool RDPIncludeCorners = false;
    // Use L2 ball with PDHG
    bool useL2Ball = true;
    // Save the sensitivity image
    bool saveSens = false;
    // Use 64-bit atomic functions
    bool use_64bit_atomics = false;
    // Use 32-bit atomic functions
    bool use_32bit_atomics = false;
    // Compute the full sensitity image when using listmode data
    bool compute_sensitivity_image = false;
    // Whether FDK weights are used or not
    bool useFDKWeights = true;
    // If true, index-based listmode reconstruction is used
    bool useIndexBasedReconstruction = false;
    // Compute the selected algorithm/prior
    bool OSEM = false;
    bool LSQR = false;
    bool CGLS = false;
    bool SART = false;
    bool FISTA = false;
    bool FISTAL1 = false;
    bool MRAMLA = false;
    bool RAMLA = false;
    bool ROSEM = false;
    bool RBI = false;
    bool DRAMA = false;
    bool COSEM = false;
    bool ECOSEM = false;
    bool ACOSEM = false;
    bool OSL_OSEM = false;
    bool MBSREM = false;
    bool BSREM = false;
    bool ROSEM_MAP = false;
    bool OSL_RBI = false;
    bool OSL_COSEM = false;
    bool PKMA = false;
    bool SPS = false;
    bool PDHG = false;
    bool PDHGKL = false;
    bool PDHGL1 = false;
    bool PDDY = false;
    bool CV = false;
    bool POCS = false;
    bool FDK = false;
    bool MRP = false;
    bool quad = false;
    bool Huber = false;
    bool L = false;
    bool FMH = false;
    bool weighted_mean = false;
    bool TV = false;
    bool hyperbolic = false;
    bool AD = false;
    bool APLS = false;
    bool TGV = false;
    bool NLM = false;
    bool RDP = false;
    bool GGMRF = false;
    bool ProxTV = false;
    bool ProxRDP = false;
    bool ProxNLM = false;
    bool MAP = false;
    bool custom = false;
    uint64_t mDim;
    uint64_t nIterSaved;
    uint64_t sizeScat;
    uint64_t eFOV;
    uint64_t sizeX;
    uint64_t sizeZ;
    uint64_t sizeAtten;
    uint64_t sizeNorm;
    uint64_t sizePSF;
    uint64_t sizeXYind;
    uint64_t sizeZind;
    uint64_t xCenterSize;
    uint64_t yCenterSize;
    uint64_t zCenterSize;
    uint64_t sizeV;
    uint64_t measElem;
    // Detector coordinates (CT) or transaxial detector coordinates (PET) or all detector coordinates (custom scanner/listmode data)
    float* x;
    // 
    float* z;
    // Unused
    float* uV;
    float* dx;
    float* dy;
    float* dz;
    float* bx;
    float* by;
    float* bz;
    float* atten;
    float* norm;
    int64_t* pituus;
    uint32_t* xy_index;
    uint16_t* z_index;
    float* x_center;
    float* y_center;
    float* z_center;
    float* V;
    float* gaussPSF;
    uint32_t* saveNiter;
    uint32_t* Nx;
    uint32_t* Ny;
    uint32_t* Nz;
//#ifdef MTYPE
//    uint16_t* randoms;
//#else
    float* randoms;
//#endif
    float* corrVector;
    float* x0;
    float* offsetVal;
    float* dScaleX4;
    float* dScaleY4;
    float* dScaleZ4;
    float* dSizeX;
    float* dSizeY;
    float* dScaleX;
    float* dScaleY;
    float* dScaleZ;
    float* kerroin4;
    float* lam_drama;
    uint8_t* maskFP;
    uint8_t* maskBP;
    uint8_t* eFOVIndices;
    uint8_t* maskPrior;
    float* angles;
    uint32_t* blurPlanes;
    float* gFilter;
    uint64_t* gFSize;
    bool* precondTypeImage;
    bool* precondTypeMeas;
    float* referenceImage;
    float* filterIm;
    float* filter;
    float* filter2;
    float* Ffilter;
    float* s;
    float* weights_quad;
    float* weights_huber;
    float* weighted_weights;
    float* APLS_ref_image;
    float* lambdaN;
    float* lambdaFiltered;
    float* alpha_PKMA;
    float* alphaPrecond;
    float* NLM_ref;
    float* tauCP;
    float* tauCPFilt;
    float* sigmaCP;
    float* sigma2CP;
    float* thetaCP;
    float* TOFCenter;
    float* TV_ref;
    uint16_t* trIndices;
    uint16_t* axIndices;
};

void copyStruct(inputStruct& options, structForScalars& inputScalars, Weighting& w_vec, RecMethods& MethodList) {
    MethodList.OSEM = options.OSEM;
    MethodList.RAMLA = options.RAMLA;
    MethodList.MRAMLA = options.MRAMLA;
    MethodList.ROSEM = options.ROSEM;
    MethodList.RBI = options.RBI;
    MethodList.DRAMA = options.DRAMA;
    MethodList.COSEM = options.COSEM;
    MethodList.ECOSEM = options.ECOSEM;
    MethodList.ACOSEM = options.ACOSEM;
    MethodList.LSQR = options.LSQR;
    MethodList.CGLS = options.CGLS;
    MethodList.SART = options.SART;
    MethodList.FISTA = options.FISTA;
    MethodList.FISTAL1 = options.FISTAL1;
    if (MethodList.LSQR || MethodList.CGLS)
        MethodList.initAlg = true;

    // Priors
    MethodList.MRP = options.MRP;
    MethodList.Quad = options.quad;
    MethodList.Huber = options.Huber;
    MethodList.L = options.L;
    MethodList.FMH = options.FMH;
    MethodList.WeightedMean = options.weighted_mean;
    MethodList.TV = options.TV;
    MethodList.hyperbolic = options.hyperbolic;
    MethodList.AD = options.AD;
    MethodList.APLS = options.APLS;
    MethodList.TGV = options.TGV;
    MethodList.NLM = options.NLM;
    MethodList.RDP = options.RDP;
    MethodList.GGMRF = options.GGMRF;
    MethodList.ProxTV = options.ProxTV;
    MethodList.ProxTGV = options.TGV;
    MethodList.ProxRDP = options.ProxRDP;
    MethodList.ProxNLM = options.ProxNLM;

    // MAP/prior-based algorithms
    MethodList.OSLOSEM = options.OSL_OSEM;
    MethodList.BSREM = options.BSREM;
    MethodList.MBSREM = options.MBSREM;
    MethodList.ROSEMMAP = options.ROSEM_MAP;
    MethodList.RBIOSL = options.OSL_RBI;
    MethodList.OSLCOSEM = options.OSL_COSEM;
    MethodList.PKMA = options.PKMA;
    MethodList.SPS = options.SPS;
    MethodList.PDHG = options.PDHG;
    MethodList.PDHGKL = options.PDHGKL;
    MethodList.PDHGL1 = options.PDHGL1;
    MethodList.CV = options.CV;
    MethodList.PDDY = options.PDDY;
    MethodList.POCS = options.POCS;

    // Whether MAP/prior-based algorithms are used
    MethodList.MAP = options.MAP;

    // Custom prior
    MethodList.CUSTOM = options.custom;

    MethodList.FDK = options.FDK;

    // Primal-dual algorithms
    if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
        MethodList.CPType = true;

    inputScalars.PET = options.PET;
    inputScalars.CT = options.CT;
    inputScalars.SPECT = options.SPECT;
    inputScalars.pitch = options.pitch;
    inputScalars.enforcePositivity = options.enforcePositivity;
    inputScalars.multiResolution = options.useMultiResolutionVolumes;
    inputScalars.nMultiVolumes = options.nMultiVolumes;
    inputScalars.indexBased = options.useIndexBasedReconstruction;

    // The number of x-voxels in the estimated image
    uint32_t* Nx = options.Nx;
    inputScalars.Nx = std::vector<uint32_t>(Nx, Nx + inputScalars.nMultiVolumes + 1);

    // The number of y-voxels in the estimated image
    uint32_t* Ny = options.Ny;
    inputScalars.Ny = std::vector<uint32_t>(Ny, Ny + inputScalars.nMultiVolumes + 1);

    // The number of z-voxels in the estimated image
    uint32_t* Nz = options.Nz;
    inputScalars.Nz = std::vector<uint32_t>(Nz, Nz + inputScalars.nMultiVolumes + 1);

    // The size of x-voxels in the estimated image
    float* dx = options.dx;
    inputScalars.dx = std::vector<float>(dx, dx + inputScalars.nMultiVolumes + 1);

    // The size of y-voxels in the estimated image
    float* dy = options.dy;
    inputScalars.dy = std::vector<float>(dy, dy + inputScalars.nMultiVolumes + 1);

    // The size of z-voxels in the estimated image
    float* dz = options.dz;
    inputScalars.dz = std::vector<float>(dz, dz + inputScalars.nMultiVolumes + 1);

    // The distance from the origin to the corner of the image (x-direction)
    float* bx = options.bx;
    inputScalars.bx = std::vector<float>(bx, bx + inputScalars.nMultiVolumes + 1);

    // The distance from the origin to the corner of the image (y-direction)
    float* by = options.by;
    inputScalars.by = std::vector<float>(by, by + inputScalars.nMultiVolumes + 1);

    // The distance from the origin to the corner of the image (z-direction)
    float* bz = options.bz;
    inputScalars.bz = std::vector<float>(bz, bz + inputScalars.nMultiVolumes + 1);

    // The size of the first dimension in the input sinogram/projection
    inputScalars.nRowsD = options.nRowsD;

    inputScalars.verbose = options.verbose;

    // Detector pair numbers, for raw data
    //const uint16_t* L = options.prhs;
    //const size_t numRows = mxGetM(prhs[ind]);
    //inputScalars.sizeL = mxGetNumberOfElements(prhs[ind]);

    // Is TOF data used?
    inputScalars.TOF = options.TOF;

    // Variance of the Gaussian TOF
    inputScalars.sigma_x = options.sigma_x;

    // Centers of the TOF-bins
    inputScalars.TOFCenter = options.TOFCenter;

    // Index offset for TOF subsets
    inputScalars.nBins = options.TOF_bins;

    inputScalars.raw = options.use_raw_data;

    inputScalars.use_psf = options.use_psf;

    // Is the attenuation correction included
    inputScalars.attenuation_correction = options.attenuation_correction;

    // Is the normalization correction included
    inputScalars.normalization_correction = options.normalization_correction;

    inputScalars.Niter = options.Niter;

    inputScalars.subsets = options.subsets;

    inputScalars.epps = options.epps;


    inputScalars.tube_width = options.tube_width_z;

    // Center coordinates of voxels in the X-dimension
    inputScalars.x_center = options.x_center;
    inputScalars.size_center_x = options.xCenterSize;

    // Center coordinates of voxels in the Y-dimension
    inputScalars.y_center = options.y_center;
    inputScalars.size_center_y = options.yCenterSize;

    // Center coordinates of voxels in the Z-dimension
    inputScalars.z_center = options.z_center;
    inputScalars.size_center_z = options.zCenterSize;

    // Randoms corrections
    inputScalars.randoms_correction = options.randoms_correction;

    // The type of projector used (Siddon or orthogonal)
    inputScalars.projector_type = options.projector_type;

    // Number of rays in Siddon
    inputScalars.n_rays = options.n_rays_transaxial;

    // Number of rays in Siddon (axial)
    inputScalars.n_rays3D = options.n_rays_axial;

    // Number of time steps
    inputScalars.Nt = options.Nt;

    // Use 64-bit integer atomic functions if possible
    inputScalars.atomic_64bit = options.use_64bit_atomics;

    inputScalars.bmin = options.bmin;

    inputScalars.bmax = options.bmax;

    inputScalars.Vmax = options.Vmax;

    inputScalars.V = options.V;
    inputScalars.size_V = options.sizeV;

    inputScalars.gaussian = options.gaussPSF;

    if (inputScalars.verbose >= 3) {
        mexPrint("Copied inputs");
    }

    inputScalars.saveIter = options.save_iter;
    inputScalars.saveNIter = options.saveNiter;
    inputScalars.saveIterationsMiddle = options.nIterSaved;
    size_t Ni = 0ULL;
    if (inputScalars.saveIter)
        Ni = static_cast<size_t>(inputScalars.Niter);
    else if (inputScalars.saveIterationsMiddle > 0)
        Ni = inputScalars.saveIterationsMiddle;
    const size_t outSize = static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * static_cast<size_t>(inputScalars.Nz[0]);
    const size_t outSize2 = Ni + 1ULL;

    // Output dimensions
    const size_t dim[5] = { static_cast<size_t>(inputScalars.Nx[0]), static_cast<size_t>(inputScalars.Ny[0]), static_cast<size_t>(inputScalars.Nz[0]), static_cast<size_t>(outSize2), static_cast<size_t>(inputScalars.Nt) };

    if (inputScalars.projector_type == 1 || inputScalars.projector_type == 11 || inputScalars.projector_type == 14 || inputScalars.projector_type == 15 || inputScalars.projector_type == 12 || inputScalars.projector_type == 13)
        inputScalars.FPType = 1;
    else if (inputScalars.projector_type == 2 || inputScalars.projector_type == 21 || inputScalars.projector_type == 24 || inputScalars.projector_type == 22 || inputScalars.projector_type == 23)
        inputScalars.FPType = 2;
    else if (inputScalars.projector_type == 3 || inputScalars.projector_type == 31 || inputScalars.projector_type == 34 || inputScalars.projector_type == 32 || inputScalars.projector_type == 33)
        inputScalars.FPType = 3;
    else if (inputScalars.projector_type == 41 || inputScalars.projector_type == 4 || inputScalars.projector_type == 42 || inputScalars.projector_type == 43 || inputScalars.projector_type == 45)
        inputScalars.FPType = 4;
    else if (inputScalars.projector_type == 51 || inputScalars.projector_type == 5 || inputScalars.projector_type == 54)
        inputScalars.FPType = 5;
    else if (inputScalars.projector_type == 6)
        inputScalars.FPType = 6;

    if (inputScalars.projector_type == 11 || inputScalars.projector_type == 41 || inputScalars.projector_type == 51 || inputScalars.projector_type == 21 || inputScalars.projector_type == 31 || inputScalars.projector_type == 1)
        inputScalars.BPType = 1;
    else if (inputScalars.projector_type == 12 || inputScalars.projector_type == 42 || inputScalars.projector_type == 22 || inputScalars.projector_type == 32 || inputScalars.projector_type == 2)
        inputScalars.BPType = 2;
    else if (inputScalars.projector_type == 13 || inputScalars.projector_type == 43 || inputScalars.projector_type == 23 || inputScalars.projector_type == 33 || inputScalars.projector_type == 3)
        inputScalars.BPType = 3;
    else if (inputScalars.projector_type == 14 || inputScalars.projector_type == 4 || inputScalars.projector_type == 24 || inputScalars.projector_type == 34 || inputScalars.projector_type == 54)
        inputScalars.BPType = 4;
    else if (inputScalars.projector_type == 15 || inputScalars.projector_type == 5 || inputScalars.projector_type == 45)
        inputScalars.BPType = 5;
    else if (inputScalars.projector_type == 6)
        inputScalars.BPType = 6;

    inputScalars.det_per_ring = options.det_per_ring;
    inputScalars.rings = options.rings;
    inputScalars.global_factor = options.global_factor;
    inputScalars.flat = options.flat;
    inputScalars.useMAD = options.useMAD;
    inputScalars.useImages = options.useImages;
    inputScalars.listmode = options.listmode;
    inputScalars.relaxScaling = options.relaxationScaling;
    inputScalars.computeRelaxation = options.computeRelaxationParameters;
    inputScalars.computeSensImag = options.compute_sensitivity_image;
    inputScalars.CT = options.CT;
    inputScalars.atomic_32bit = options.use_32bit_atomics;
    inputScalars.scatter = options.additionalCorrection;
    inputScalars.CTAttenuation = options.CTAttenuation;
    inputScalars.largeDim = options.largeDim;
    inputScalars.loadTOF = options.loadTOF;
    inputScalars.storeResidual = options.storeResidual;
    if (inputScalars.scatter == 1U) {
        inputScalars.size_scat = options.sizeScat;
    }
    if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
        inputScalars.meanFP = options.meanFP;
        inputScalars.meanBP = options.meanBP;
    }
    if (inputScalars.BPType == 5 || inputScalars.BPType == 4) {
        inputScalars.maskBP = options.useMaskBP;
    }
    inputScalars.maskFP = options.useMaskFP;
    inputScalars.offset = options.offsetCorrection;
    if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        inputScalars.orthXY = options.orthTransaxial;
        inputScalars.orthZ = options.orthAxial;
        if (inputScalars.BPType == 3 || inputScalars.FPType == 3)
            inputScalars.cylRadiusProj3 = options.tube_radius;
    }
    if (inputScalars.offset)
        inputScalars.T = options.offsetVal;
    //inputScalars.T = options.OffsetLimit;
    inputScalars.nProjections = options.nProjections;
    inputScalars.subsetType = options.subsetType;
    if (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 4 || inputScalars.BPType == 5) {
        inputScalars.dL = options.dL;
        inputScalars.d_Scale4.resize(inputScalars.nMultiVolumes + 1);
        inputScalars.dSize.resize(inputScalars.nMultiVolumes + 1);
        inputScalars.d_Scale.resize(inputScalars.nMultiVolumes + 1);
        float* dScaleX4 = options.dScaleX4;
        float* dScaleY4 = options.dScaleY4;
        float* dScaleZ4 = options.dScaleZ4;
        float* dSizeX = nullptr, * dSizeY = nullptr, * dScaleX = nullptr, * dScaleY = nullptr, * dScaleZ = nullptr;
        if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
            dSizeX = options.dSizeX;
            dSizeY = options.dSizeY;
            //float* dSizeZ = options.dSizeZ;
            dScaleX = options.dScaleX;
            dScaleY = options.dScaleY;
            dScaleZ = options.dScaleZ;
        }
        for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
            inputScalars.d_Scale4[ii].x = dScaleX4[ii];
            inputScalars.d_Scale4[ii].y = dScaleY4[ii];
            inputScalars.d_Scale4[ii].z = dScaleZ4[ii];
            //inputScalars.d_Scale.s[3] = 0.f;
            if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
                inputScalars.dSize[ii].x = dSizeX[ii];
                inputScalars.dSize[ii].y = dSizeY[ii];
                //inputScalars.dSize[ii].z = dSizeZ[ii];
                inputScalars.d_Scale[ii].x = dScaleX[ii];
                inputScalars.d_Scale[ii].y = dScaleY[ii];
                inputScalars.d_Scale[ii].z = dScaleZ[ii];
                if (ii == 0) {
                    inputScalars.dSizeBP.x = options.dSizeXBP;
                    inputScalars.dSizeBP.y = options.dSizeZBP;
                }
            }
        }
    }
    if (!inputScalars.CT && !inputScalars.SPECT)
        inputScalars.nLayers = options.nLayers;
    inputScalars.Nf = options.Nf;
    inputScalars.useExtendedFOV = options.useEFOV;
    if (inputScalars.useExtendedFOV)
        inputScalars.eFOV = options.eFOV > 1;
    inputScalars.NxOrig = options.NxOrig;
    inputScalars.NyOrig = options.NyOrig;
    inputScalars.NzOrig = options.NzOrig;
    inputScalars.NxPrior = options.NxPrior;
    inputScalars.NyPrior = options.NyPrior;
    inputScalars.NzPrior = options.NzPrior;
    inputScalars.TGV2D = options.use2DTGV;
    inputScalars.adaptiveType = options.PDAdaptiveType;
    inputScalars.storeFP = options.storeFP;
    //const uint32_t* devPointer = options.use_device;
    //size_t devLength = mxGetNumberOfElements(mxGetField(options, 0, "use_device"));
    //inputScalars.usedDevices = std::vector<uint32_t>(devPointer, devPointer + devLength);
    if (inputScalars.CT) {
        inputScalars.nColsD = options.nColsD;
        inputScalars.nRowsD = options.nRowsD;
    }
    else {
        inputScalars.nColsD = options.Nang;
        inputScalars.nRowsD = options.Ndist;
    }
    if (inputScalars.use_psf) {
        inputScalars.g_dim_x = options.g_dim_x;
        inputScalars.g_dim_y = options.g_dim_y;
        inputScalars.g_dim_z = options.g_dim_z;
        inputScalars.deconvolution = options.deblurring;
    }
    if (inputScalars.use_psf && inputScalars.deconvolution) {
        inputScalars.deblur_iterations = options.deblur_iterations;
    }

    inputScalars.Nxy = inputScalars.Nx[0] * inputScalars.Ny[0];
    inputScalars.im_dim[0] = static_cast<int64_t>(inputScalars.Nxy) * static_cast<int64_t>(inputScalars.Nz[0]);
    if (inputScalars.multiResolution) {
        for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++)
            inputScalars.im_dim[ii] = static_cast<int64_t>(inputScalars.Nx[ii]) * static_cast<int64_t>(inputScalars.Ny[ii]) * static_cast<int64_t>(inputScalars.Nz[ii]);
    }
	
    Ni = 1;
    // Number of voxels
    if (inputScalars.saveIter)
        Ni = inputScalars.Niter + 1U;
    // Load the necessary variables if the corresponding reconstruction method is used
    int yy = 0;

    if (MethodList.DRAMA) {
        // Relaxation parameter
        w_vec.lambda = options.lam_drama;
    }

    // Load regularization parameter
    w_vec.beta = options.beta;
    w_vec.betaReg = w_vec.beta;

    // Masks
    if (inputScalars.maskFP)
        w_vec.maskFP = options.maskFP;
    if (inputScalars.maskBP && (inputScalars.BPType == 4 || inputScalars.BPType == 5)) {
        w_vec.maskBP = options.maskBP;
        //const size_t nMask = mxGetNumberOfElements(getField(options, 0, "maskBP"));
        //if (DEBUG) {
        //    mexPrintBase("nMask = %d\n", nMask);
        //    mexEval();
        //}
    }
    //if (inputScalars.offset)
    //	w_vec.maskOffset = options.maskOffset;
    if (inputScalars.eFOV && inputScalars.useExtendedFOV && !inputScalars.multiResolution)
        w_vec.eFOVIndices = options.eFOVIndices;
    if (inputScalars.useExtendedFOV && !inputScalars.multiResolution)
        w_vec.maskPrior = options.maskPrior;
    else if (inputScalars.maskBP)
        w_vec.maskPrior = options.maskBP;
    //if (MethodList.NLM || MethodList.MRP || MethodList.CPTV || MethodList.CPTVL1 || MethodList.CPTVKL || MethodList.RDP || (MethodList.TV && !w_vec.data.TV_use_anatomical)
    //	|| MethodList.CPTGV || MethodList.CPTGVKL || MethodList.CPTGVL1)
    // CT-related variables such as number of projection images
    if (inputScalars.CT) {
        w_vec.nProjections = options.nProjections;
        w_vec.dPitchX = options.dPitchX;
        w_vec.dPitchY = options.dPitchY;
    }
    else {
        w_vec.nProjections = options.nProjections;
        // Detector pitch
        w_vec.dPitchX = options.cr_p;
        w_vec.dPitchY = options.cr_pz;
    }
    if (inputScalars.FPType == 4 || inputScalars.BPType == 4)
        w_vec.kerroin4 = options.kerroin4;

#ifdef AF
    if (inputScalars.projector_type == 6U) {
        const uint64_t* ind = options.gFSize;
        if (DEBUG) {
            mexPrintBase("indX = %d\n", ind[0]);
            mexPrintBase("indY = %d\n", ind[1]);
            mexEval();
        }
        w_vec.angles = options.angles;
        w_vec.gFilter = af::array(ind[0], ind[1], ind[2], ind[3], options.gFilter);
        w_vec.distInt = options.blurPlanes;
        if (DEBUG) {
            mexPrintBase("w_vec.gFilter.dims(0) = %d\n", w_vec.gFilter.dims(0));
            mexPrintBase("w_vec.gFilter.dims(1) = %d\n", w_vec.gFilter.dims(1));
            mexPrintBase("w_vec.gFilter.dims(2) = %d\n", w_vec.gFilter.dims(2));
            mexPrintBase("w_vec.gFilter.dims(3) = %d\n", w_vec.gFilter.dims(3));
            mexPrintBase("w_vec.distInt[0] = %d\n", w_vec.distInt[0]);
            mexEval();
        }
        if (DEBUG) {
            mexPrint("SPECT vars loaded");
        }
    }
#endif

    // True value means the preconditioner is included
    // precondTypeIm[0] = Diagonal normalization preconditioner (division with the sensitivity image 1 / (A^T1), A is the system matrix)
    // precondTypeIm[1] = EM preconditioner (f / (A^T1), where f is the current estimate)
    // precondTypeIm[2] = IEM preconditioner (max(n, fhat, f)/ (A^T1), where fhat is an estimate of the final image and n is a small positive number)
    // precondTypeIm[3] = Momentum-like preconditioner (basically a step size inclusion)
    // precondTypeIm[4] = Gradient-based preconditioner (Uses the normalized divergence (sum of the gradient) of the current estimate)
    // precondTypeIm[5] = Filtering-based preconditioner
    // precondTYpeIm[6] = Curvature-based preconditioner
    // The first three are mutually exclusive, but otherwise they can be mixed and matched
    const bool* precondIm = options.precondTypeImage;
    w_vec.precondTypeIm.assign(precondIm, precondIm + w_vec.precondTypeIm.size());
    const bool* precondMe = options.precondTypeMeas;
    w_vec.precondTypeMeas.assign(precondMe, precondMe + w_vec.precondTypeMeas.size());
    w_vec.filteringOrig = w_vec.precondTypeMeas[1];
    // precondTypeMeas[0] = Diagonal normalization preconditioner (1 / (A1))
    // precondTypeMeas[1] = Filtering-based preconditioner
#ifdef AF
    if (w_vec.precondTypeIm[2]) {
        if (DEBUG) {
            mexPrintBase("w_vec.precondTypeIm[2] = %u\n", w_vec.precondTypeIm[2]);
            mexEval();
        }
        w_vec.preRef.resize(inputScalars.nMultiVolumes + 1);
        const float* ref = options.referenceImage;
        w_vec.preRef[0] = af::array(inputScalars.im_dim[0], ref);
        if (inputScalars.multiResolution) {
            for (int kk = 1; kk <= inputScalars.nMultiVolumes; kk++) {
                size_t dimApu = inputScalars.im_dim[kk - 1];
                w_vec.preRef[kk] = af::array(inputScalars.im_dim[kk], &ref[dimApu]);
            }
        }
        if (DEBUG) {
            mexPrint("Precond im [2] loaded");
        }
    }

    if (w_vec.precondTypeIm[4]) {
        if (DEBUG) {
            mexPrintBase("w_vec.precondTypeIm[4] = %u\n", w_vec.precondTypeIm[4]);
            mexEval();
        }
        w_vec.gradV1 = options.gradV1;
        w_vec.gradV2 = options.gradV2;
        w_vec.gradInitIter = options.gradInitIter;
        w_vec.gradFinalIter = options.gradLastIter;
        w_vec.gradF.resize(inputScalars.nMultiVolumes + 1);
        if (DEBUG) {
            mexPrint("Precond im [4] loaded");
        }
    }

    if (w_vec.precondTypeIm[5]) {
        if (DEBUG) {
            mexPrintBase("w_vec.precondTypeIm[5] = %u\n", w_vec.precondTypeIm[5]);
            mexEval();
        }
        w_vec.filterIm = af::array(inputScalars.Nf, inputScalars.Nf, options.filterIm);
        w_vec.filterIter = options.filteringIterations;
        if (DEBUG) {
            mexPrint("Precond im [5] loaded");
        }
    }

    if (w_vec.precondTypeMeas[1]) {
        if (DEBUG) {
            mexPrintBase("w_vec.precondTypeMeas[1] = %u\n", w_vec.precondTypeMeas[1]);
            mexEval();
        }
        if (FINVERSE)
            w_vec.filter = af::array(inputScalars.Nf, options.filter);
        else
            if (inputScalars.subsetType == 5)
                w_vec.filter = af::array(inputScalars.nColsD, options.filter2);
            else
                w_vec.filter = af::array(inputScalars.nRowsD, options.filter2);
        w_vec.Ffilter = af::array(inputScalars.Nf, options.Ffilter);
        w_vec.filterIter = options.filteringIterations;
        if (DEBUG) {
            mexPrint("Precond meas [1] loaded");
        }
    }

    // The complete sensitivity image is computed
    if (MethodList.RBI || MethodList.RBIOSL || MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || (w_vec.precondTypeIm[0]
        || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2]))
        w_vec.computeD = true;

    if (w_vec.precondTypeMeas[0] || MethodList.SART || MethodList.POCS)
        w_vec.computeM = true;
#endif

    // Load TV related input data
    if (MethodList.TV && MethodList.MAP) {
        // Is anatomical reference image used
        w_vec.data.TV_use_anatomical = options.TV_use_anatomical;
        // Tau-value
        w_vec.data.tau = options.tau;
        // "Smoothing" parameter, prevents zero values in the square root
        w_vec.data.TVsmoothing = options.TVsmoothing;
        // The type of TV prior used
        w_vec.data.TVtype = options.TVtype;
        // If anatomical prior is used, load the necessary coefficients
#ifdef AF
        if (w_vec.data.TV_use_anatomical) {
            if (w_vec.data.TVtype == 1) {
                w_vec.data.refIm = af::array(inputScalars.im_dim[0], options.s, afHost);
            }
            else {
                w_vec.data.refIm = af::array(inputScalars.im_dim[0], options.TV_ref, afHost);
            }
            w_vec.data.T = options.T;
            w_vec.data.C = options.C;
        }
#endif
        if (w_vec.data.TVtype == 4 || w_vec.data.TVtype == 6)
            w_vec.data.SATVPhi = options.SATVPhi;
        if (DEBUG) {
            mexPrintBase("w_vec.data.TVtype = %u\n", w_vec.data.TVtype);
            mexPrintBase("w_vec.data.SATVPhi = %f\n", w_vec.data.SATVPhi);
            mexEval();
        }
    }
    // General variables for neighborhood-based methods
    if ((MethodList.L || MethodList.FMH || MethodList.WeightedMean || MethodList.Quad || MethodList.Huber || MethodList.MRP || MethodList.NLM || MethodList.ProxNLM || MethodList.hyperbolic || MethodList.RDP)
        && MethodList.MAP) {
        // Neighborhood size
        w_vec.Ndx = options.Ndx;
        w_vec.Ndy = options.Ndy;
        w_vec.Ndz = options.Ndz;
        // Is normalization used in MRP, FMH, L, weighted mean or AD
        w_vec.med_no_norm = options.med_no_norm;
        w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
        if (DEBUG) {
            mexPrint("Neighborhood loaded");
        }
    }
//#ifdef AF
//    if ((MethodList.L || MethodList.FMH) && MethodList.MAP) {
//        // Index values for the neighborhood
//        w_vec.tr_offsets = af::array(inputScalars.im_dim[0], w_vec.dimmu, getUint32s(options, "tr_offsets", 0), afHost);
//    }
//#endif
    if (MethodList.FMH || MethodList.Quad || MethodList.Huber)
        w_vec.inffi = options.inffi;
    // Weights for the quadratic prior
#ifdef AF
    if (MethodList.Quad && MethodList.MAP) {
        w_vec.weights_quad = af::array(w_vec.dimmu - 1, options.weights_quad, afHost);
        int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
        af::array weight = w_vec.weights_quad * -1.f;
        af::array w_quad = af::constant(0.f, pituus);
        w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
        w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
        w_quad(pituus / 2) = af::abs(af::sum(weight));
        w_vec.weights_quad = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
        if (DEBUG) {
            mexPrint("Quad loaded");
        }
    }
    if (MethodList.hyperbolic && MethodList.MAP) {
        w_vec.weights = options.weights_quad;;
        w_vec.data.SATVPhi = options.hyperbolicDelta;
    }
#endif
    if ((MethodList.RDP && MethodList.MAP) || MethodList.ProxRDP) {
        w_vec.RDP_gamma = options.RDP_gamma;
        w_vec.RDPLargeNeighbor = options.RDPIncludeCorners;
        if (DEBUG) {
            mexPrint("RDP loaded");
        }
    }
    if (MethodList.GGMRF) {
        w_vec.GGMRF_p = options.GGMRF_p;
        w_vec.GGMRF_q = options.GGMRF_q;
        w_vec.GGMRF_c = options.GGMRF_c;
        w_vec.GGMRF_pqc = (w_vec.GGMRF_p - w_vec.GGMRF_q) / std::pow(w_vec.GGMRF_c, w_vec.GGMRF_p - w_vec.GGMRF_q);
        w_vec.weights = options.weights_quad;
        if (DEBUG) {
            mexPrint("GGMRF loaded");
        }
    }
#ifdef AF
    if (MethodList.Huber && MethodList.MAP) {
        w_vec.weights_huber = af::array(w_vec.dimmu - 1, options.weights_huber, afHost);
        int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
        af::array weight = w_vec.weights_huber * -1.f;
        af::array w_quad = af::constant(0.f, pituus);
        w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
        w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
        w_quad(pituus / 2) = af::abs(af::sum(weight));
        w_vec.weights_huber = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
        w_vec.huber_delta = options.huber_delta;
        if (DEBUG) {
            mexPrint("Huber loaded");
        }
    }
    //if (MethodList.L && MethodList.MAP)
    //    w_vec.a_L = af::array(w_vec.dimmu, options.a_L, afHost);
    //if (MethodList.FMH && MethodList.MAP) {
    //    if (inputScalars.Nz[0] == 1 || w_vec.Ndz == 0)
    //        w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, options.fmh_weights, afHost);
    //    else
    //        w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, options.fmh_weights, afHost);
    //    w_vec.alku_fmh = options.inffi;
    //}
    if (MethodList.WeightedMean && MethodList.MAP) {
        w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, options.weighted_weights, afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
        // Type of mean used (arithmetic, harmonic or geometric)
        w_vec.mean_type = options.mean_type;
        // Sum of the weights
        w_vec.w_sum = options.w_sum;
        if (DEBUG) {
            mexPrint("WMean loaded");
        }
    }
    if (MethodList.AD && MethodList.MAP) {
        // Time-step value
        w_vec.TimeStepAD = options.TimeStepAD;
        // Conductance (edge value)
        w_vec.KAD = options.KAD;
        // Number of AD iterations
        w_vec.NiterAD = options.NiterAD;
        // Flux type
        uint32_t Flux = options.FluxType;
        // Diffusion type
        uint32_t Diffusion = options.DiffusionType;
        if (Flux == 2U)
            w_vec.FluxType = AF_FLUX_QUADRATIC;
        else
            w_vec.FluxType = AF_FLUX_EXPONENTIAL;
        if (Diffusion == 2U)
            w_vec.DiffusionType = AF_DIFFUSION_MCDE;
        else
            w_vec.DiffusionType = AF_DIFFUSION_GRAD;
        if (!MethodList.L && !MethodList.FMH && !MethodList.WeightedMean && !MethodList.Quad && !MethodList.MRP)
            w_vec.med_no_norm = options.med_no_norm;
    }
    // Asymmetric parallel level sets
    if (MethodList.APLS && MethodList.MAP) {
        // Eta value
        w_vec.data.eta = options.eta;
        // Tau-value
        if (!MethodList.TV)
            w_vec.data.tau = options.tau;
        // Smoothing value
        w_vec.data.TVsmoothing = options.APLSsmoothing;
        // Anatomical reference image
        w_vec.data.refIm = af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], options.APLS_ref_image, afHost);
        w_vec.data.TV_use_anatomical = true;
        w_vec.data.TVtype = 5;
    }
#endif
    if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS) {
        // Relaxation parameter
        w_vec.lambda = options.lambdaN;
        // Upper bound
        w_vec.U = options.U;
    }
    // Relaxation parameters
    if (MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.SART || MethodList.POCS)
        w_vec.lambda = options.lambdaN;
    if (MethodList.PKMA) {
        w_vec.alphaM = options.alpha_PKMA;
        w_vec.lambda = options.lambdaN;
    }
    if ((w_vec.precondTypeIm[5] || w_vec.precondTypeMeas[1]) && (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA)) {
        w_vec.lambdaFiltered = w_vec.lambda;
        w_vec.lambda = options.lambdaFiltered;
    }
    if (DEBUG && (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA)) {
        mexPrintBase("w_vec.lambda[0] = %f\n", w_vec.lambda[0]);
        mexEval();
    }
    if (w_vec.precondTypeIm[3])
        w_vec.alphaPrecond = options.alphaPrecond;
    // Power factor for ACOSEM
    if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1) {
        w_vec.h_ACOSEM = options.h_ACOSEM;
        w_vec.h_ACOSEM_2 = 1.f / w_vec.h_ACOSEM;
    }
    //if (MethodList.TGV && MethodList.MAP) {
    //	w_vec.data.TGVAlpha = options.alpha1TGV;
    //	w_vec.data.TGVBeta = options.beta;
    //	w_vec.data.NiterTGV = options.NiterTGV;
    //}
    if ((MethodList.NLM && MethodList.MAP) || MethodList.ProxNLM) {
        w_vec.NLM_anatomical = options.NLM_use_anatomical;
        w_vec.NLTV = options.NLTV;
        w_vec.NLRD = options.NLRD;
        w_vec.NLM_MRP = options.NLM_MRP;
        w_vec.NLLange = options.NLLange;
        w_vec.NLGGMRF = options.NLGGMRF;
        if (w_vec.NLRD)
            w_vec.RDP_gamma = options.RDP_gamma;
        else if (w_vec.NLLange)
            w_vec.RDP_gamma = options.SATVPhi;
        else if (w_vec.NLGGMRF) {
            w_vec.GGMRF_p = options.GGMRF_p;
            w_vec.GGMRF_q = options.GGMRF_q;
            w_vec.GGMRF_c = options.GGMRF_c;
            w_vec.RDP_gamma = (w_vec.GGMRF_p - w_vec.GGMRF_q) / std::pow(w_vec.GGMRF_c, w_vec.GGMRF_p - w_vec.GGMRF_q);
        }
        if (w_vec.NLM_anatomical)
            w_vec.NLM_ref = options.NLM_ref;;
        w_vec.h2 = options.NLMsigma;
        w_vec.h2 = w_vec.h2 * w_vec.h2;
        w_vec.Nlx = options.Nlx;
        w_vec.Nly = options.Nly;
        w_vec.Nlz = options.Nlz;
#ifdef AF
        w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), options.gaussPSF, afHost);
#endif
        if (DEBUG) {
            mexPrint("NLM loaded");
        }
    }
    if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1) {
        w_vec.tauCP = options.tauCPFilt;
        w_vec.sigmaCP = options.sigmaCP;
        w_vec.powerIterations = options.powerIterations;
        if (DEBUG) {
            mexPrint("PIter loaded");
        }
    }
    if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
        if (w_vec.precondTypeMeas[1])
            w_vec.tauCP2 = options.tauCP;
        else
            w_vec.tauCP = options.tauCP;
        w_vec.sigma2CP = options.sigma2CP;
        w_vec.betaReg = options.beta;
        w_vec.thetaCP = options.thetaCP;;
        w_vec.alpha0CPTGV = options.alpha0TGV;
        w_vec.alpha1CPTGV = options.alpha1TGV;
        w_vec.UseL2Ball = options.useL2Ball;
        if (inputScalars.adaptiveType == 1) {
            for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
                //w_vec.alphaCP.emplace_back(.3f);
                w_vec.alphaCP.emplace_back(1.f);
        }
        if (DEBUG) {
            mexPrint("CPType loaded");
        }
    }
    inputScalars.useFDKWeights = options.useFDKWeights;
    if (MethodList.FDK && inputScalars.CT && inputScalars.useFDKWeights) {
        w_vec.angles = options.angles;
        inputScalars.DSC = options.sourceToCRot;
    }
    w_vec.derivType = options.derivativeType;
    if (DEBUG) {
        mexPrintBase("options.useFDKWeights = %d\n", options.useFDKWeights);
        mexPrintBase("options.compute_sensitivity_image = %d\n", options.compute_sensitivity_image);
        mexPrintBase("options.OSEM = %d\n", options.OSEM);
        mexEval();
    }
    if (MethodList.POCS) {
        w_vec.ng = options.POCS_NgradIter;
        w_vec.alphaPOCS = options.POCS_alpha;
        w_vec.rMaxPOCS = options.POCS_rMax;
        w_vec.POCSalphaRed = options.POCS_alphaRed;
        w_vec.POCSepps = options.POCSepps;
    }
}

// Transfers the device data to host
// Transfer the ArrayFire arrays from the device to the host pointers
void device_to_host(const RecMethods& MethodList, AF_im_vectors& vec, int64_t& oo, float* output, float* FPoutput, const scalarStruct& inputScalars,
    std::vector<std::vector<std::vector<float>>>& FPEstimates) {
    if (inputScalars.storeFP) {
        size_t dim = 0ULL;
        for (uint32_t ii = 0; ii < inputScalars.subsets * inputScalars.Niter; ii++) {
            const uint32_t jj = ii % inputScalars.subsets;
            const uint32_t kk = ii / inputScalars.subsets;
            std::copy(FPEstimates[kk][jj].begin(), FPEstimates[kk][jj].end(), FPoutput + dim);
            dim += FPEstimates[kk][jj].size();
        }
		if (DEBUG) {
			mexPrintBase("dim = %d\n", dim);
                mexEval();
		}
    }
    // Transfer data back to host
    if (CELL) {
        for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
            if (DEBUG) {
                mexPrintBase("inputScalars.Nx[ii] = %d\n", inputScalars.Nx[ii]);
                mexPrintBase("inputScalars.Ny[ii] = %d\n", inputScalars.Ny[ii]);
                mexPrintBase("inputScalars.Nz[ii] = %d\n", inputScalars.Nz[ii]);
                mexEval();
            }
            if (inputScalars.saveIter || inputScalars.saveIterationsMiddle > 0) {
            }
            else {
                if (MethodList.FDK)
                    vec.rhs_os[ii].host(&output[oo]);
                else
                    vec.im_os[ii].host(&output[oo]);
                if (inputScalars.verbose >= 3)
                    mexPrint("Data transfered to host");
                oo += inputScalars.im_dim[ii];
            }
        }
    }
    else {
        if (inputScalars.saveIter || inputScalars.saveIterationsMiddle > 0) {
        }
        else {
            if (MethodList.FDK && inputScalars.largeDim) {
                //vec.rhs_os[0].host(&output[oo]);
            }
            else if (MethodList.FDK && !inputScalars.largeDim) {
                vec.rhs_os[0].host(&output[oo]);
            }
            else
                vec.im_os[0].host(&output[oo]);
            if (inputScalars.verbose >= 3)
                mexPrint("Data transfered to host");
            oo += inputScalars.im_dim[0];
        }
    }
    af::sync();
}

extern "C" DLL_FUNCTION
#ifdef MTYPE
int omegaMain(inputStruct options, const char* header_directory, const uint16_t * Sino, float* outputPtr, float* FPptr = nullptr, float* residual = nullptr);
#else
int omegaMain(inputStruct options, const char* header_directory, const float* Sino, float* outputPtr, float* FPptr = nullptr, float* residual = nullptr);
#endif