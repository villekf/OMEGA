#include "arrayfire.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	const char * c = af::infoString();

	mexPrintf("\n%s\n", c);

	plhs[0] = mxCreateString(c);
}