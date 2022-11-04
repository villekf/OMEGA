#pragma once
#ifdef CUDA

#else
#define CTYPE3 cl_float3
#define CTYPE2 cl_float2
#include "opencl_error.hpp"
#endif

typedef struct structForScalars {
	uint32_t projector_type = 1, attenuation_correction = 0, randoms_correction = 0, scatter = 0, normalization_correction = 0, Nx, Ny, Nz,
		size_y, size_x, size_z, subsets = 1, TotSinos = 1, det_per_ring, pRows, dec = 1, Niter = 1, Nt = 1, cSizeX = 0, cSizeY = 0, subsetType = 0,
		nLayers = 0, nRekos = 0, osa_iter0 = 0, nRekos2 = 0, NSinos = 0, subsetsUsed = 1, TOFsubsets = 1, im_dim, Nxy;
	float crystal_size_z = 0.f, cr_pz = 0.f, dx, dy, dz, dPitch = 0.f, epps = 1e-8f, bx = 0.f, by = 0.f, bz = 0.f, bzb = 0.f, maxxx = 0.f,
		maxyy = 0.f, zmax = 0.f, NSlices, sigma_x = 0.f, tube_width = 0.f, bmin = 0.f, bmax = 0.f, Vmax = 0.f, dc_z = 0.f, global_factor = 1.f,
		dL = 0.f, tStep = 0.f, tStart = 0.f, cThickness = 0.f, detY = 0.f, meanV = 0.f, crystal_size_xy = 0.f;
	bool precompute = false, use_psf = false, TOF = false, SPECT = false, PITCH = false, PET = false, meanFP = false, meanBP = false,
		maskFP = false, maskBP = false, orthXY = false, orthZ = false, CT = false, atomic_64bit = false, atomic_32bit = false, loadTOF = true, 
		saveIter = false;
	CTYPE3 d_Scale = { 0.f, 0.f, 0.f, 0.f };
	CTYPE3 dSize = { 0.f, 0.f, 0.f, 0.f };
	CTYPE2 dSizeBP = { 0.f, 0.f };
	uint8_t raw = 0, fp = 0, listmode = 0, verbose = 0;
	uint16_t n_rays = 1, n_rays3D = 1;
	int64_t nBins = 1, nProjections = 0, numelY = 0, numelZ = 0, TOFSize = 0;
	size_t size_of_x, size_atten = 1, size_norm = 1, size_center_x, size_center_y, size_center_z, size_V = 1, size_scat = 1, koko = 0;
	//float* meanVals = nullptr;
} scalarStruct;

void precomp_siddon(const cl_uint& num_devices_context, const cl::Context& context, const std::vector<cl::CommandQueue>& commandQueues, uint16_t* lor1, const float* z_det,
	const float* x, const float* y, scalarStruct inputScalars, const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, const uint32_t* pseudos,
	const uint16_t* L, const uint32_t im_dim, const cl::Kernel& kernel, const size_t numel_x, const size_t local_size[]);