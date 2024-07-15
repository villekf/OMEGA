#include "arrayfire.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	//af::array A = af::constant(1.f, 3, 3, 2);
	//A = af::sat(A);
	//const char* b = af::toString("A", A, 4, false);
	const char * c = af::infoString();
	//mexPrintf("%s\n", b);

	plhs[0] = mxCreateString(c);
}