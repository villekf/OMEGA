#include "mex.h"
#include <stdint.h>
#include <cmath>
#include <TROOT.h>
#include "TChain.h"


void histogram(uint16_t * LL1, uint16_t * LL2, uint32_t * tpoints, double vali, const double alku, const double loppu, const size_t outsize2, 
	const uint32_t detectors, const bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring, 
	const uint32_t det_per_ring, float* S, float* output, TTree *Coincidences, const int64_t Nentries)
{

	Int_t crystalID1, crystalID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
	Float_t sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2;
	Double_t time1, time2;

	Coincidences->SetBranchAddress("crystalID1", &crystalID1);
	Coincidences->SetBranchAddress("crystalID2", &crystalID2);
	Coincidences->SetBranchAddress("moduleID1", &moduleID1);
	Coincidences->SetBranchAddress("moduleID2", &moduleID2);
	Coincidences->SetBranchAddress("rsectorID1", &rsectorID1);
	Coincidences->SetBranchAddress("rsectorID2", &rsectorID2);
	Coincidences->SetBranchAddress("sourcePosX1", &sourcePosX1);
	Coincidences->SetBranchAddress("sourcePosX2", &sourcePosX2);
	Coincidences->SetBranchAddress("sourcePosY1", &sourcePosY1);
	Coincidences->SetBranchAddress("sourcePosY2", &sourcePosY2);
	Coincidences->SetBranchAddress("sourcePosZ1", &sourcePosZ1);
	Coincidences->SetBranchAddress("sourcePosZ2", &sourcePosZ2);
	Coincidences->SetBranchAddress("time1", &time1);
	Coincidences->SetBranchAddress("time2", &time2);

	Int_t nbytes = 0;

	for (int64_t kk = 0; kk < Nentries; kk++) {

		nbytes += Coincidences->GetEntry(kk);

		uint32_t ring_number1 = (moduleID1 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
		uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring)*cryst_per_block + (crystalID1 % cryst_per_block);
		uint32_t ring_number2 = (moduleID2 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
		uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring)*cryst_per_block + (crystalID2 % cryst_per_block);
		uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
		uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
		if (outsize2 == 1) {
			LL1[L1*detectors + L2] = LL1[L1*detectors + L2] + static_cast<uint16_t>(1);
		}
		else {
			LL1[kk] = static_cast<uint16_t>(L1 + 1);
			LL2[kk] = static_cast<uint16_t>(L2 + 1);
		}
		if (time2 >= alku + vali) {
			tpoints[kk] = kk;
			vali += vali;
		}
		if (source) {
			S[kk] = sourcePosX1;
			S[kk + Nentries] = sourcePosY1;
			S[kk + Nentries * 2] = sourcePosZ1;
			S[kk + Nentries * 3] = sourcePosX2;
			S[kk + Nentries * 4] = sourcePosY2;
			S[kk + Nentries * 5] = sourcePosZ2;
		}
	}
	return;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 10) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:invalidNumInputs",
			"10 input arguments required.");
	}
	else if (nlhs > 5) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:maxlhs",
			"Too many output arguments.");
	}

	/* Create a matrix for the return argument */
	double vali = (double)mxGetScalar(prhs[1]);
	double alku = (double)mxGetScalar(prhs[2]);
	double loppu = (double)mxGetScalar(prhs[3]);
	uint32_t detectors = (uint32_t)mxGetScalar(prhs[4]);
	uint32_t blocks_per_ring = (uint32_t)mxGetScalar(prhs[5]);
	uint32_t cryst_per_block = (uint32_t)mxGetScalar(prhs[6]);
	uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[7]);
	uint32_t linear_multp = (uint32_t)mxGetScalar(prhs[8]);
	bool source = (bool)mxGetScalar(prhs[9]);
	size_t outsize2 = (loppu - alku) / vali;

	// Count inputs and check for char type
	char *argv;
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	argv = mxArrayToString(prhs[0]);

	TChain *Coincidences = new TChain("Coincidences");
	Coincidences->Add(argv);


	int64_t Nentries = Coincidences->GetEntries();

	if (outsize2 == 1) {
		plhs[0] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}
	else {
		plhs[0] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
	}
	plhs[2] = mxCreateNumericMatrix(outsize2, 1, mxUINT32_CLASS, mxREAL);
	if (source)
		plhs[3] = mxCreateNumericMatrix(Nentries, 6, mxSINGLE_CLASS, mxREAL);
	else
		plhs[3] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	/* Assign pointers to the various parameters */
	uint16_t * LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t * LL2 = (uint16_t*)mxGetData(plhs[1]);
	uint32_t * tpoints = (uint32_t*)mxGetData(plhs[2]);
	float* S = 0;
	float * output = 0;
	if (source)
		S = (float*)mxGetData(plhs[3]);
	else
		output = (float*)mxGetData(plhs[3]);


	histogram(LL1, LL2, tpoints, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, 
		output, Coincidences, Nentries);

	delete Coincidences;
	gROOT->Reset();

	return;
}
