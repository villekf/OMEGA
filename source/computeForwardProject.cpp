#include "AF_opencl_functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void forwardProjection(const RecMethods& MethodList, array& y, array& input, const size_t length, const scalarStruct& inputScalars, 
	const Weighting& w_vec) {


	//af::array y = constant(0.f, length, 1);
	if (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.RBI || MethodList.ECOSEM ||
		MethodList.ROSEM || MethodList.RAMLA || MethodList.BSREM || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.RBIOSL ||
		MethodList.PKMA) {
		if (inputScalars.CT) {
			input = exp(-input) / y;
		}
		else
			input = y / input;
	}
	else if (MethodList.LSQR) {

	}
	else if (MethodList.CP) {

	}
	else if (MethodList.MBSREM) {

	}

	af::sync();
	//d_Sino = cl::Buffer(*y.device<cl_mem>(), true);
}