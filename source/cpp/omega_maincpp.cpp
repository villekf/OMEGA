#if defined(_MSC_VER)
#define DLL_FUNCTION __declspec(dllexport)
#endif
#include "libHeader.h"

#if defined(MTYPE) && !defined(MTYPE2)
int omegaMain(inputStruct options, const char* header_directory, const uint16_t* Sino, float* outputPtr, float* FPptr, float* residual) {
#elif !defined(MTYPE) && defined(MTYPE2)
int omegaMain(inputStruct options, const char* header_directory, const uint8_t * Sino, float* outputPtr, float* FPptr, float* residual) {
#else
int omegaMain(inputStruct options, const char* header_directory, const float* Sino, float* outputPtr, float* FPptr, float* residual) {
#endif

	scalarStruct inputScalars;

	af::setDevice(options.deviceNum);

	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;

	// Struct containing the necessary variables for the priors
	Weighting w_vec;


	if (DEBUG) {
		mexPrint("Copying data\n");
	}

	copyStruct(options, inputScalars, w_vec, MethodList);

	// Coordinates of the detectors in z-direction (PET/SPECT) or the directional vectors for the detector panel pixels (CT)
	const float* z_det = options.z;
	inputScalars.size_z = options.sizeZ;

	// Coordinates of the detectors in x/y-directions
	float* x = options.x;
	inputScalars.size_of_x = options.sizeX;

	// attenuation values
	const float* atten = options.atten;
	inputScalars.size_atten = options.sizeAtten;

	// Normalization coefficients
	const float* norm = options.norm;
	inputScalars.size_norm = options.sizeNorm;

	// Number of measurements/LORs
	const int64_t* pituus = options.pituus;

	// XY-indices of the detector coordinates of each LOR
	const uint32_t* xy_index = options.xy_index;
	inputScalars.sizeXY = options.sizeXYind;

	// Z-indices of the detector coordinates of each LOR
	const uint16_t* z_index = options.z_index;
	inputScalars.sizeZ = options.sizeZind;

	const float* x0 = options.x0;

	// The device used
	const uint32_t device = options.deviceNum;

	const size_t size_gauss = options.sizePSF;

	const float* extraCorr = options.corrVector;

//#ifdef MTYPE
//	const uint16_t* randoms = options.randoms;
//#else
	const float* randoms = options.randoms;
//#endif

	size_t mDim = options.measElem / static_cast<size_t>(inputScalars.Nt);

	if (inputScalars.listmode) {
		if (inputScalars.indexBased) {
			w_vec.trIndex = options.trIndices;
			w_vec.axIndex = options.axIndices;
		}
		else {
			w_vec.listCoord = options.x;
			x = options.uV;
		}
		if (inputScalars.TOF)
			w_vec.TOFIndices = options.TOFIndices;
		inputScalars.kokoNonTOF = mDim;
	}
	else
		inputScalars.kokoNonTOF = mDim / inputScalars.nBins;
	inputScalars.kokoTOF = mDim;
	
	
	if (DEBUG) {
		mexPrintBase("kokoTOF = %u\n", inputScalars.kokoTOF);
		mexPrintBase("kokoNonTOF = %u\n", inputScalars.kokoNonTOF);
		mexPrintBase("size_z = %u\n", inputScalars.size_z);
		mexPrintBase("inputScalars.largeDim = %u\n", inputScalars.largeDim);
		mexPrintBase("inputScalars.maskBP = %u\n", inputScalars.maskBP);
		mexPrintBase("inputScalars.maskFP = %u\n", inputScalars.maskFP);
		mexPrintBase("inputScalars.offset = %u\n", inputScalars.offset);
		mexPrintBase("inputScalars.projector_type = %u\n", inputScalars.projector_type);
		mexPrintBase("inputScalars.FPType = %u\n", inputScalars.FPType);
		mexPrintBase("inputScalars.BPType = %u\n", inputScalars.BPType);
		mexPrintBase("inputScalars.useExtendedFOV = %u\n", inputScalars.useExtendedFOV);
		mexPrintBase("inputScalars.eFOV = %u\n", inputScalars.eFOV);
		mexPrintBase("inputScalars.TGV2D = %u\n", inputScalars.TGV2D);
		mexPrintBase("inputScalars.NxPrior = %u\n", inputScalars.NxPrior);
		mexPrintBase("inputScalars.NyPrior = %u\n", inputScalars.NyPrior);
		mexPrintBase("inputScalars.NzPrior = %u\n", inputScalars.NzPrior);
		mexPrintBase("inputScalars.im_dim = %u\n", inputScalars.im_dim[0]);
		mexPrintBase("inputScalars.Nx = %u\n", inputScalars.Nx[0]);
		mexPrintBase("inputScalars.Ny = %u\n", inputScalars.Ny[0]);
		mexPrintBase("inputScalars.Nz = %u\n", inputScalars.Nz[0]);
		mexPrintBase("inputScalars.Nf = %u\n", inputScalars.Nf);
		mexPrintBase("inputScalars.Nt = %u\n", inputScalars.Nt);
		mexPrintBase("options.measElem = %u\n", options.measElem);
		mexPrintBase("mDim = %u\n", mDim);
		mexPrintBase("inputScalars.nColsD = %u\n", inputScalars.nColsD);
		mexPrintBase("inputScalars.nRowsD = %u\n", inputScalars.nRowsD);
		mexPrintBase("inputScalars.bmin = %f\n", inputScalars.bmin);
		mexPrintBase("inputScalars.bmax = %f\n", inputScalars.bmax);
		mexPrintBase("inputScalars.Vmax = %f\n", inputScalars.Vmax);
		mexPrintBase("inputScalars.size_V = %u\n", inputScalars.size_V);
		mexPrintBase("MethodList.FDK = %u\n", MethodList.FDK);
		mexPrintBase("w_vec.dPitchX = %f\n", w_vec.dPitchX);
		mexEval();
	}

	if (inputScalars.verbose >= 3) {
		mexPrint("Loaded struct values. Starting reconstruction itself...");
	}
	int status = reconstructionAF(z_det, x, Sino, randoms, inputScalars, device, pituus, w_vec, MethodList, header_directory, x0,
		outputPtr, FPptr, atten, norm, extraCorr, size_gauss, xy_index, z_index, residual);

	if (status != 0)
		mexPrint("Reconstruction failed!");

	fflush(stdout);

	return 0;
}