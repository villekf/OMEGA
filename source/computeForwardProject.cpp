#include "structs.h"
// Use ArrayFire namespace for convenience
using namespace af;

void computeForwardStep(const RecMethods& MethodList, af::array& y, af::array& input, const int64_t length, const scalarStruct& inputScalars,
	Weighting& w_vec, const af::array& randomsData) {

	if (inputScalars.randoms_correction)
		input += randomsData;
	//af::array y = constant(0.f, length, 1);
	if (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.RBI || MethodList.ECOSEM ||
		MethodList.ROSEM || MethodList.RAMLA || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.RBIOSL ||
		MethodList.PKMA) {
		if (inputScalars.CT) {
			input = exp(-input) / y;
		}
		else
			input = y / (input + inputScalars.epps);
	}
	else if (MethodList.RAMLA || MethodList.BSREM) {
		if (inputScalars.CT) {
			input = exp(-input) / y - 1.f;
		}
		else
			input = y / (input + inputScalars.epps) - 1.f;
	}
	else if (MethodList.LSQR) {
		input -= w_vec.alphaLSQR * y;
		w_vec.betaLSQR = norm(input);
		input = input / w_vec.betaLSQR;
		af::array Dy = input;
		y = Dy.copy();
	}
	else if (MethodList.CGLS) {

	}
	else if (MethodList.CPLS) {

	}
	else if (MethodList.CPTV) {

	}
	else if (MethodList.MBSREM) {

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