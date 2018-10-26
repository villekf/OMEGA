#include <fstream>
#include "mex.h"
#include <cmath>


FILE *streami;


void histogram(uint16_t * LL1, uint16_t * LL2, uint64_t * tpoints, char *argv, const double vali, const double alku, const double loppu, const size_t outsize2, 
	const uint32_t detectors, const int header_bytes, const int R_bits, const int M_bits, const int S_bits, const int C_bits, const int L_bits, 
	int data_bytes, const int R_length, const int M_length, const int C_length, const double coincidence_window, const bool source, const uint32_t linear_multip, 
	const uint32_t cryst_per_block, const uint32_t blocks_per_ring, const uint32_t det_per_ring, int16_t* S, const size_t pituus, uint64_t* output)
{

	static int64_t qb;
	unsigned char ew1[8];
	static int tag;

	int i1 = -1;
	int64_t i = 1;
	int ll = -1;

	uint64_t ms = 0;
	uint64_t ms2 = 0;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	fopen_s(&streami, argv, "rb");
#else
	streami = fopen(argv, "rb");
#endif
	if (streami == NULL) {
		mexErrMsgIdAndTxt("MATLAB:gate_lmf_matlab:invalidFile",
			"Error opening file");
	}

	bool warning = false;
	uint32_t R;    //detector 1
	uint32_t M;    //detector 2
	uint32_t C;
	uint64_t time = 0;
	int oo = 0;
	data_bytes = data_bytes - 8;

	uint64_t window = static_cast<uint64_t>(coincidence_window);
	fseek(streami, header_bytes, SEEK_SET);
	while (i = fread(&ew1, 1, 8, streami) != 0) {

		i1++;
		ll++;


		uint64_t num = (uint64_t)ew1[7] | (uint64_t)ew1[6] << 8 | (uint64_t)ew1[5] << 16 | (uint64_t)ew1[4] << 24 | (uint64_t)ew1[3] << 32 | (uint64_t)ew1[2] << 40 | (uint64_t)ew1[1] << 48 | (uint64_t)ew1[0] << 56;
		//uint64_t num = (uint64_t)ew1[0] | (uint64_t)ew1[1] << 8 | (uint64_t)ew1[2] << 16 | (uint64_t)ew1[3] << 24 | (uint64_t)ew1[4] << 32 | (uint64_t)ew1[5] << 40 | (uint64_t)ew1[6] << 48 | (uint64_t)ew1[7] << 56;

		tag = ((num << (63)) & 1);



		if (tag != 0 && warning == false) {
			mexPrintf("Tag %d\n", tag);
			mexWarnMsgTxt("Tag bit not zero, make sure header and data bytes are correct");
			warning = true;
		}

		int16_t X, Y, Z;
		
		if (tag == 0) {
			if (ll == 0) {
				ms = (num & 0xFFFFFFFFFFFFFFFF);
				i = fread(&ew1, 1, data_bytes, streami);
				uint64_t num = (uint64_t)ew1[7] | (uint64_t)ew1[6] << 8 | (uint64_t)ew1[5] << 16 | (uint64_t)ew1[4] << 24 | (uint64_t)ew1[3] << 32 | (uint64_t)ew1[2] << 40 | (uint64_t)ew1[1] << 48 | (uint64_t)ew1[0] << 56;
				R = ((num >> 64 - R_bits) & R_length);
				M = (num >> 64 - R_bits - M_bits) & M_length;
				C = (num >> 64 - R_bits - M_bits - S_bits - C_bits) & C_length;
				if (source) {
					X = (num >> 32) & 0xFFFF;
					Y = (num >> 16) & 0xFFFF;
					Z = num & 0xFFFF;
				}
			}
			else {
				ms2 = (num & 0xFFFFFFFFFFFFFFFF);
				i = fread(&ew1, 1, data_bytes, streami);
				uint64_t num = (uint64_t)ew1[7] | (uint64_t)ew1[6] << 8 | (uint64_t)ew1[5] << 16 | (uint64_t)ew1[4] << 24 | (uint64_t)ew1[3] << 32 | (uint64_t)ew1[2] << 40 | (uint64_t)ew1[1] << 48 | (uint64_t)ew1[0] << 56;
				if (ms2 - ms <= window) {
					uint32_t ring_number1 = (M % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(C) / static_cast<double>(cryst_per_block)));
					uint32_t ring_pos1 = (R % blocks_per_ring)*cryst_per_block + (C % cryst_per_block);
					R = ((num >> 64 - R_bits) & R_length);
					M = (num >> 64 - R_bits - M_bits) & M_length;
					C = (num >> 64 - R_bits - M_bits - S_bits - C_bits) & C_length;
					uint32_t ring_number2 = (M % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(C) / static_cast<double>(cryst_per_block)));
					uint32_t ring_pos2 = (R % blocks_per_ring)*cryst_per_block + (C % cryst_per_block);
					uint32_t L1 = ring_number1*det_per_ring + ring_pos1;
					uint32_t L2 = ring_number2*det_per_ring + ring_pos2;
					if (outsize2 == 1) {
						LL1[L1*detectors + L2] = LL1[L1*detectors + L2] + static_cast<uint16_t>(1);
					}
					else {
						LL1[i1] = static_cast<uint16_t>(L1 + 1);
						LL2[i1] = static_cast<uint16_t>(L2 + 1);
					}
					ll = -1;
					if (source) {
						S[i1] = (num >> 32) & 0xFFFF;
						S[i1 + pituus] = (num >> 16) & 0xFFFF;
						S[i1 + pituus * 2] = num & 0xFFFF;
						S[i1 + pituus * 3] = X;
						S[i1 + pituus * 4] = Y;
						S[i1 + pituus * 5] = Z;
					}
				}
				else {
					R = ((num >> 64 - R_bits) & R_length);
					M = (num >> 64 - R_bits - M_bits) & M_length;
					C = (num >> 64 - R_bits - M_bits - S_bits - C_bits) & C_length;
					if (source) {
						X = (num >> 32) & 0xFFFF;
						Y = (num >> 16) & 0xFFFF;
						Z = num & 0xFFFF;
					}
					ms = ms2;
				}
				time = ms2;
			}
		}
	}
	fclose(streami);
	return;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 22) {
		mexErrMsgIdAndTxt("MATLAB:gate_lmf_matlab:invalidNumInputs",
			"22 input arguments required.");
	}
	else if (nlhs > 5) {
		mexErrMsgIdAndTxt("MATLAB:gate_lmf_matlab:maxlhs",
			"Too many output arguments.");
	}

	/* Create a matrix for the return argument */
	double vali = (double)mxGetScalar(prhs[1]);
	double alku = (double)mxGetScalar(prhs[2]);
	double loppu = (double)mxGetScalar(prhs[3]);
	size_t pituus = (size_t)mxGetScalar(prhs[4]);
	uint32_t detectors = (uint32_t)mxGetScalar(prhs[5]);
	uint32_t blocks_per_ring = (uint32_t)mxGetScalar(prhs[6]);
	uint32_t cryst_per_block = (uint32_t)mxGetScalar(prhs[7]);
	uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[8]);
	uint32_t linear_multp = (uint32_t)mxGetScalar(prhs[9]);
	int header_bytes = (int)mxGetScalar(prhs[10]);
	int data_bytes = (int)mxGetScalar(prhs[11]);
	int R_bits = (int)mxGetScalar(prhs[12]);
	int M_bits = (int)mxGetScalar(prhs[13]);
	int S_bits = (int)mxGetScalar(prhs[14]);
	int C_bits = (int)mxGetScalar(prhs[15]);
	int L_bits = (int)mxGetScalar(prhs[16]);
	int R_length = (int)mxGetScalar(prhs[17]);
	int M_length = (int)mxGetScalar(prhs[18]);
	int C_length = (int)mxGetScalar(prhs[19]);
	bool source = (bool)mxGetScalar(prhs[20]);
	double coincidence_window = (double)mxGetScalar(prhs[21]);
	size_t outsize2 = (loppu - alku) / vali;
	if (outsize2 == 1) {
		plhs[0] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}
	else {
		plhs[0] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
	}
	plhs[2] = mxCreateNumericMatrix(outsize2, 1, mxUINT64_CLASS, mxREAL);

	if (source)
		plhs[3] = mxCreateNumericMatrix(pituus, 6, mxINT16_CLASS, mxREAL);
	else
		plhs[3] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);

	/* Assign pointers to the various parameters */
	uint16_t * LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t * LL2 = (uint16_t*)mxGetData(plhs[1]);
	uint64_t * tpoints = (uint64_t*)mxGetData(plhs[2]);
	int16_t* S = 0;
	uint64_t * output = 0;
	if (source)
		S = (int16_t*)mxGetData(plhs[3]);
	else
		output = (uint64_t*)mxGetData(plhs[3]);

	// Check for char type
	char *argv;

	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	argv = mxArrayToString(prhs[0]);

	/* Do the actual computations in a subroutine */
	histogram(LL1, LL2, tpoints, argv, vali, alku, loppu, outsize2, detectors, header_bytes, R_bits, M_bits, S_bits, C_bits, L_bits, data_bytes, R_length, M_length,
		C_length, coincidence_window, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, pituus, output);

	return;
}
