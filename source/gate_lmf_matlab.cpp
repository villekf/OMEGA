/**************************************************************************
* This file loads the LMF binary data and converts it into the detector
* pair numbers. 
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include <fstream>
#include "mex.h"
#include <cmath>


FILE *streami;


void histogram(uint16_t * LL1, uint16_t * LL2, uint64_t * tpoints, char *argv, const double vali, const double alku, const double loppu, const size_t outsize2, 
	const uint32_t detectors, const int header_bytes, const int R_bits, const int M_bits, const int S_bits, const int C_bits, const int L_bits, 
	int data_bytes, const int R_length, const int M_length, const int C_length, const double coincidence_window, const bool source, const uint32_t linear_multip, 
	const uint32_t cryst_per_block, const uint32_t blocks_per_ring, const uint32_t det_per_ring, int16_t* S, const size_t pituus, uint64_t* output, const double time_step,
	int* int_loc, const double* time_intervals, bool obtain_trues, bool store_randoms, bool store_scatter, uint16_t* Ltrues, uint16_t* Lscatter, uint16_t* Lrandoms,
	uint8_t* trues_loc)
{

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
	uint32_t R; // R-sector
	uint32_t M; // Modules
	uint32_t C; // Crystals
	uint64_t time = 0;
	uint32_t eventID1;
	uint32_t eventID2;
	uint8_t n_ComptonP1;
	uint8_t n_ComptonP2;
	int oo = 0;
	data_bytes = data_bytes - 8;
	double aika = alku;
	bool begin = false;
	if (outsize2 > 1)
		begin = true;
	int pa = 0;
	int jj = 0;
	tpoints[jj] = 0;
	jj++;
	uint64_t num;

	uint64_t window = static_cast<uint64_t>(coincidence_window);
	fseek(streami, header_bytes, SEEK_SET);
	while (i = fread(&ew1, 1, 8, streami) != 0) {

		i1++;
		ll++;

		int bytes_left = data_bytes;


		num = (uint64_t)ew1[7] | (uint64_t)ew1[6] << 8 | (uint64_t)ew1[5] << 16 | (uint64_t)ew1[4] << 24 | (uint64_t)ew1[3] << 32 | (uint64_t)ew1[2] << 40 | (uint64_t)ew1[1] << 48 | (uint64_t)ew1[0] << 56;
		//uint64_t num = (uint64_t)ew1[0] | (uint64_t)ew1[1] << 8 | (uint64_t)ew1[2] << 16 | (uint64_t)ew1[3] << 24 | (uint64_t)ew1[4] << 32 | (uint64_t)ew1[5] << 40 | (uint64_t)ew1[6] << 48 | (uint64_t)ew1[7] << 56;

		tag = ((num << (63)) & 1);



		if (tag != 0 && warning == false) {
			mexPrintf("Tag %d\n", tag);
			mexWarnMsgTxt("Tag bit not zero, make sure header and data bytes are correct");
			warning = true;
		}

		int16_t X, Y, Z;
		
		if (tag == 0) {
			// First single
			if (ll == 0) {
				// Time
				ms = (num & 0xFFFFFFFFFFFFFFFF);
				// Continue if the time is less than the start time
				if ((static_cast<double>(ms) * time_step) < alku) {
					ll = -1;
					i = fread(&ew1, 1, data_bytes, streami);
					continue;
				}
				// Stop if the time is more than the end time
				else if ((static_cast<double>(ms) * time_step) > loppu)
					break;
				//i = fread(&ew1, 1, data_bytes, streami);
				i = fread(&ew1, 1, 2, streami);
				bytes_left -= 2;
				//if (data_bytes == 8)
					//num = (uint64_t)ew1[7] | (uint64_t)ew1[6] << 8 | (uint64_t)ew1[5] << 16 | (uint64_t)ew1[4] << 24 | (uint64_t)ew1[3] << 32 | (uint64_t)ew1[2] << 40 | (uint64_t)ew1[1] << 48 | (uint64_t)ew1[0] << 56;
				num = (uint64_t)ew1[1] | (uint64_t)ew1[0] << 8;
				R = ((num >> (16 - R_bits)) & R_length);
				M = (num >> (16 - R_bits - M_bits)) & M_length;
				C = (num >> (16 - R_bits - M_bits - S_bits - C_bits)) & C_length;
				if (obtain_trues || store_randoms || store_scatter) {
					i = fread(&ew1, 1, 4, streami);
					eventID1 = (uint32_t)ew1[3] | (uint32_t)ew1[2] << 8 | (uint32_t)ew1[1] << 16 | (uint32_t)ew1[0] << 24;
					bytes_left -= 4;
				}
				if (source) {
					i = fread(&ew1, 1, 6, streami);
					num = (uint64_t)ew1[5] | (uint64_t)ew1[4] << 8 | (uint64_t)ew1[3] << 16 | (uint64_t)ew1[2] << 24 | (uint64_t)ew1[1] << 32 | (uint64_t)ew1[0] << 40;
					X = (num >> 32) & 0xFFFF;
					Y = (num >> 16) & 0xFFFF;
					Z = num & 0xFFFF;
					bytes_left -= 6;
				}
				if (store_scatter || obtain_trues) {
					i = fread(&ew1, 1, 1, streami);
					n_ComptonP1 = (uint8_t)ew1[0];
					bytes_left -= 1;
				}
				if (bytes_left > 0)
					i = fread(&ew1, 1, bytes_left, streami);
			}
			// If a single has already been detected, check the second single
			else {
				ms2 = (num & 0xFFFFFFFFFFFFFFFF);
				i = fread(&ew1, 1, 2, streami);
				bytes_left -= 2;
				num = (uint64_t)ew1[1] | (uint64_t)ew1[0] << 8;
				// If within coincidence, form a coincidence event
				if (ms2 - ms <= window) {
					uint32_t ring_number1 = (M % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(C) / static_cast<double>(cryst_per_block)));
					uint32_t ring_pos1 = (R % blocks_per_ring)*cryst_per_block + (C % cryst_per_block);
					R = ((num >> (16 - R_bits)) & R_length);
					M = (num >> (16 - R_bits - M_bits)) & M_length;
					C = (num >> (16 - R_bits - M_bits - S_bits - C_bits)) & C_length;
					uint32_t ring_number2 = (M % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(C) / static_cast<double>(cryst_per_block)));
					uint32_t ring_pos2 = (R % blocks_per_ring)*cryst_per_block + (C % cryst_per_block);
					uint32_t L1 = ring_number1*det_per_ring + ring_pos1;
					uint32_t L2 = ring_number2*det_per_ring + ring_pos2;
					if (obtain_trues || store_randoms || store_scatter) {
						i = fread(&ew1, 1, 4, streami);
						bytes_left -= 4;
						eventID2 = (uint32_t)ew1[3] | (uint32_t)ew1[2] << 8 | (uint32_t)ew1[1] << 16 | (uint32_t)ew1[0] << 24;
						if (outsize2 == 1 && store_randoms) {
							if (eventID1 != eventID2)
								Lrandoms[L1*detectors + L2] = Lrandoms[L1*detectors + L2] + static_cast<uint16_t>(1);
						}
						else if (store_randoms) {
							if (eventID1 != eventID2)
								Lrandoms[i1] = 1;
						}
					}
					if (source) {
						i = fread(&ew1, 1, 6, streami);
						bytes_left -= 6;
						num = (uint64_t)ew1[5] | (uint64_t)ew1[4] << 8 | (uint64_t)ew1[3] << 16 | (uint64_t)ew1[2] << 24 | (uint64_t)ew1[1] << 32 | (uint64_t)ew1[0] << 40;
						S[i1] = X;
						S[i1 + pituus] = Y;
						S[i1 + pituus * 2] = Z;
						S[i1 + pituus * 3] = (num >> 32) & 0xFFFF;
						S[i1 + pituus * 4] = (num >> 16) & 0xFFFF;
						S[i1 + pituus * 5] = num & 0xFFFF;
					}
					if (obtain_trues || store_scatter) {
						i = fread(&ew1, 1, 1, streami);
						n_ComptonP2 = (uint8_t)ew1[0];
						bytes_left -= 1;
						if (outsize2 == 1) {
							if (obtain_trues && eventID1 == eventID2 && n_ComptonP1 == 0 && n_ComptonP2 == 0) {
								Ltrues[L1*detectors + L2] = Ltrues[L1*detectors + L2] + static_cast<uint16_t>(1);
								if (source)
									trues_loc[i1] = 1;
							}
							else if (store_scatter && eventID1 == eventID2 && (n_ComptonP1 > 0 || n_ComptonP2 > 0))
								Lscatter[L1*detectors + L2] = Lscatter[L1*detectors + L2] + static_cast<uint16_t>(1);
						}
						else {
							if (obtain_trues && eventID1 == eventID2 && n_ComptonP1 == 0 && n_ComptonP2 == 0)
								Ltrues[i1] = 1;
							else if (store_scatter && eventID1 == eventID2 && (n_ComptonP1 > 0 || n_ComptonP2 > 0))
								Lscatter[i1] = 1;
						}
					}
					if (outsize2 == 1) {
						LL1[L1*detectors + L2] = LL1[L1*detectors + L2] + static_cast<uint16_t>(1);
						int_loc[0] = pa;
					}
					else {
						if (begin) {
							while ((static_cast<double>(ms2) * time_step) >= time_intervals[pa])
								pa++;
							//pa--;
							begin = false;
							aika = time_intervals[pa];
							int_loc[0] = pa;
						}
						if ((static_cast<double>(ms2) * time_step) >= aika) {
							tpoints[jj++] = i1;
							aika = time_intervals[++pa];
						}
						LL1[i1] = static_cast<uint16_t>(L1 + 1);
						LL2[i1] = static_cast<uint16_t>(L2 + 1);
					}
					ll = -1;
					if (bytes_left > 0)
						i = fread(&ew1, 1, bytes_left, streami);
				}
				// Otherwise discard the earlier single and replace it with this one
				else {
					R = ((num >> (16 - R_bits)) & R_length);
					M = (num >> (16 - R_bits - M_bits)) & M_length;
					C = (num >> (16 - R_bits - M_bits - S_bits - C_bits)) & C_length;
					if (obtain_trues || store_randoms || store_scatter) {
						i = fread(&ew1, 1, 4, streami);
						eventID1 = (uint32_t)ew1[3] | (uint32_t)ew1[2] << 8 | (uint32_t)ew1[1] << 16 | (uint32_t)ew1[0] << 24;
						bytes_left -= 4;
					}
					if (source) {
						i = fread(&ew1, 1, 6, streami);
						bytes_left -= 6;
						num = (uint64_t)ew1[5] | (uint64_t)ew1[4] << 8 | (uint64_t)ew1[3] << 16 | (uint64_t)ew1[2] << 24 | (uint64_t)ew1[1] << 32 | (uint64_t)ew1[0] << 40;
					}
					if (obtain_trues || store_scatter) {
						bytes_left -= 1;
						i = fread(&ew1, 1, 1, streami);
						n_ComptonP1 = (uint8_t)ew1[0];
					}
					if (source) {
						X = (num >> 32) & 0xFFFF;
						Y = (num >> 16) & 0xFFFF;
						Z = num & 0xFFFF;
					}
					if (bytes_left > 0)
						i = fread(&ew1, 1, bytes_left, streami);
					ms = ms2;
				}
				time = ms2;
			}
		}
	}
	int_loc[1] = pa;
	tpoints[jj] = i1;
	if (begin) {
		int_loc[0] = 0;
		int_loc[1] = 0;
	}
	fclose(streami);
	return;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 27) {
		mexErrMsgIdAndTxt("MATLAB:gate_lmf_matlab:invalidNumInputs",
			"27 input arguments required.");
	}
	else if (nlhs > 9) {
		mexErrMsgIdAndTxt("MATLAB:gate_lmf_matlab:maxlhs",
			"Too many output arguments.");
	}

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
	double time_step = (double)mxGetScalar(prhs[22]);
	double* time_intervals = (double*)mxGetData(prhs[23]);
	bool obtain_trues = (bool)mxGetScalar(prhs[24]);
	bool store_randoms = (bool)mxGetScalar(prhs[25]);
	bool store_scatter = (bool)mxGetScalar(prhs[26]);
	size_t outsize2 = (loppu - alku) / vali;
	if (outsize2 == 1) {
		plhs[0] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);

		if (obtain_trues) {
			plhs[5] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
			if (source)
				plhs[8] = mxCreateNumericMatrix(pituus, 1, mxUINT8_CLASS, mxREAL);
		}
		else {
			plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[8] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
		}
		if (store_randoms)
			plhs[6] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		else
			plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_scatter)
			plhs[7] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		else
			plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}
	else {
		plhs[0] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);

		if (obtain_trues) {
			plhs[5] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
			plhs[8] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
		}
		else {
			plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[8] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
		}
		if (store_randoms)
			plhs[6] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_scatter)
			plhs[7] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}
	plhs[2] = mxCreateNumericMatrix(outsize2 + 2, 1, mxUINT64_CLASS, mxREAL);

	if (source)
		plhs[3] = mxCreateNumericMatrix(pituus, 6, mxINT16_CLASS, mxREAL);
	else
		plhs[3] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
	plhs[4] = mxCreateNumericMatrix(2, 1, mxINT32_CLASS, mxREAL);

	/* Assign pointers to the various parameters */
	// Detector numbers
	uint16_t * LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t * LL2 = (uint16_t*)mxGetData(plhs[1]);
	// Time points
	uint64_t * tpoints = (uint64_t*)mxGetData(plhs[2]);
	int16_t* S = 0;
	uint64_t * output = 0;
	uint8_t * trues_loc;
	if (source)
		S = (int16_t*)mxGetData(plhs[3]);
	else
		output = (uint64_t*)mxGetData(plhs[3]);
	int *int_loc = (int*)mxGetData(plhs[4]);
	uint16_t * Ltrues = (uint16_t*)mxGetData(plhs[5]);
	uint16_t * Lscatter = (uint16_t*)mxGetData(plhs[6]);
	uint16_t * Lrandoms = (uint16_t*)mxGetData(plhs[7]);
	if (obtain_trues && source)
		trues_loc = (uint8_t*)mxGetData(plhs[8]);
	else
		trues_loc = 0;

	// Check for char type
	char *argv;

	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	argv = mxArrayToString(prhs[0]);

	/* Do the actual computations in a subroutine */
	histogram(LL1, LL2, tpoints, argv, vali, alku, loppu, outsize2, detectors, header_bytes, R_bits, M_bits, S_bits, C_bits, L_bits, data_bytes, R_length, M_length,
		C_length, coincidence_window, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, pituus, output, time_step, int_loc, time_intervals, obtain_trues,
		store_randoms, store_scatter, Ltrues, Lscatter, Lrandoms, trues_loc);

	return;
}
