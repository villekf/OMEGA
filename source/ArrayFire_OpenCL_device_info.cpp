#include "arrayfire.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	const char * c = af::infoString();

	plhs[0] = mxCreateString(c);
}