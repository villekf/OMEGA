#include "functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeForwardStep(const RecMethods& MethodList, af::array& y, af::array& input, const int64_t length, const scalarStruct& inputScalars,
	Weighting& w_vec, const af::array& randomsData, AF_im_vectors& vec) {

	if (inputScalars.randoms_correction)
		input += randomsData;
	//af::array y = constant(0.f, length, 1);
	if (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.ECOSEM ||
		MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.DRAMA) {
		if (inputScalars.CT) {
			input = exp(-input) / y;
		}
		else
			input = y / (input + inputScalars.epps);
	}
	else if (MethodList.RAMLA || MethodList.BSREM || MethodList.RBI || MethodList.RBIOSL) {
		if (inputScalars.CT) {
			input = exp(-input) / y - 1.f;
		}
		else
			input = y / (input + inputScalars.epps) - 1.f;
	}
	else if (MethodList.PKMA) {
		if (inputScalars.CT) {
			input = 1.f - exp(-input) / y;
		}
		else
			input = 1.f - y / (input + inputScalars.epps);
	}
	else if (MethodList.MBSREM || MethodList.MRAMLA) {

	}
	else if (MethodList.LSQR) {
		input -= w_vec.alphaLSQR * y;
		w_vec.betaLSQR = norm(input);
		input = input / w_vec.betaLSQR;
		af::array Dy = input;
		y = Dy.copy();
	}
	else if (MethodList.CGLS) {
		const float normi = af::norm(input);
		w_vec.alphaCGLS = w_vec.gammaCGLS / (normi * normi);
		input = vec.rCGLS - w_vec.alphaCGLS * input;
		af::array Dy = input;
		vec.rCGLS = Dy.copy();
	}
	else if (MethodList.CPLS || MethodList.CPTV) {
		af::array res = input - y;
		const af::array Dy = (vec.pCP + w_vec.sigmaCP * res) / (1.f + w_vec.sigmaCP);
		input = Dy.copy();
		vec.pCP = Dy.copy();
	}
	else if (MethodList.CPLSKL || MethodList.CPTVKL) {
		//af::array res = input - y;
		const af::array Dy = .5f * (1.f + vec.pCP + w_vec.sigmaCP * input - af::sqrt(af::pow(vec.pCP + w_vec.sigmaCP * input - 1.f, 2.) + 4.f * w_vec.sigmaCP * y));
		input = Dy.copy();
		vec.pCP = Dy.copy();
	}

	af::sync();
	//d_Sino = cl::Buffer(*y.device<cl_mem>(), true);
}

void computeIntegralImage(const scalarStruct& inputScalars, const Weighting& w_vec, const int64_t length, af::array& outputFP, af::array& meanBP) {
	if (inputScalars.projector_type == 5) {
		outputFP = af::moddims(outputFP, w_vec.size_x, w_vec.size_y, length);
		if (inputScalars.meanBP) {
			meanBP = af::mean(af::mean(outputFP, 0), 1);
			outputFP -= af::tile(meanBP, w_vec.size_x, w_vec.size_y, 1);
		}
		outputFP = af::sat(outputFP);
		outputFP = af::join(0, af::constant(0.f, 1, outputFP.dims(1), outputFP.dims(2)), outputFP);
		//if (DEBUG) {
		//	mexPrintf("outputFP.dims(0) = %u\n", outputFP.dims(0));
		//	mexPrintf("outputFP.dims(1) = %u\n", outputFP.dims(1));
		//	mexEvalString("pause(.0001);");
		//}
		outputFP = af::flat(af::join(1, af::constant(0.f, outputFP.dims(0), 1, outputFP.dims(2)), outputFP));
	}
}