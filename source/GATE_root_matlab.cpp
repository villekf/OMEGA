#include "mex.h"
#include <stdint.h>
#include <cmath>
#include <TROOT.h>
#include "TChain.h"


void histogram(uint16_t * LL1, uint16_t * LL2, uint32_t * tpoints, double vali, const double alku, const double loppu, const size_t outsize2, 
	const uint32_t detectors, const bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring, 
	const uint32_t det_per_ring, float* S, float* output, TTree *Coincidences, const int64_t Nentries, const double* time_intervals, int* int_loc, 
	const bool obtain_trues, const bool store_scatter, const bool store_randoms, const bool* scatter_components, uint16_t* Ltrues, uint16_t* Lscatter, 
	uint16_t* Lrandoms, bool* trues_loc, const int64_t Ndelays, bool randoms_correction, TTree *delay, uint16_t * Ldelay1, uint16_t * Ldelay2,
	int* int_loc_delay, uint32_t * tpoints_delay, bool* randoms_loc, bool* scatter_loc)
{

	Int_t crystalID1, crystalID2, moduleID1, moduleID2, rsectorID1, rsectorID2, eventID1, eventID2, comptonPhantom1, comptonPhantom2, 
		comptonCrystal1, comptonCrystal2, RayleighPhantom1, RayleighPhantom2, RayleighCrystal1, RayleighCrystal2;
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
	//Coincidences->SetBranchAddress("time1", &time1);
	Coincidences->SetBranchAddress("time2", &time2);
	if (obtain_trues || store_scatter || store_randoms) {
		Coincidences->SetBranchAddress("eventID1", &eventID1);
		Coincidences->SetBranchAddress("eventID2", &eventID2);
	}
	if (obtain_trues || store_scatter) {
		//if (scatter_components[0]) {
			Coincidences->SetBranchAddress("comptonPhantom1", &comptonPhantom1);
			Coincidences->SetBranchAddress("comptonPhantom2", &comptonPhantom2);
		//}
		//if (scatter_components[1]) {
			Coincidences->SetBranchAddress("comptonCrystal1", &comptonCrystal1);
			Coincidences->SetBranchAddress("comptonCrystal2", &comptonCrystal2);
		//}
		//if (scatter_components[2]) {
			Coincidences->SetBranchAddress("RayleighPhantom1", &RayleighPhantom1);
			Coincidences->SetBranchAddress("RayleighPhantom2", &RayleighPhantom2);
		//}
		//if (scatter_components[3]) {
			Coincidences->SetBranchAddress("RayleighCrystal1", &RayleighCrystal1);
			Coincidences->SetBranchAddress("RayleighCrystal2", &RayleighCrystal2);
		//}
	}

	Int_t nbytes = 0;
	int ll = 0, jj = 0;
	bool begin = false;
	if (outsize2 > 1)
		begin = true;
	int pa = 0;
	double aika = 0.;

	for (int64_t kk = 0; kk < Nentries; kk++) {

		jj++;
		nbytes += Coincidences->GetEntry(kk);

		if (time2 < alku)
			continue;
		else if (time2 > loppu)
			break;
		uint32_t ring_number1, ring_number2;

		if (linear_multip == 1) {
			ring_number1 = moduleID1;
			ring_number2 = moduleID2;
		}
		else {
			ring_number1 = (moduleID1 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
			ring_number2 = (moduleID2 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
		}
		const uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring)*cryst_per_block + (crystalID1 % cryst_per_block);
		const uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring)*cryst_per_block + (crystalID2 % cryst_per_block);
		const uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
		const uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
		if (begin) {
			while (time2 >= time_intervals[pa])
				pa++;
			pa--;
			begin = false;
			aika = time_intervals[pa];
			int_loc[0] = pa;
		}
		if (obtain_trues || store_scatter || store_randoms) {
			bool event_true = true;
			bool event_scattered = false;
			if (eventID1 != eventID2) {
				event_true = false;
			}
			if (event_true && (obtain_trues || store_scatter)) {
				//if (event_true && scatter_components[0]) {
				if (comptonPhantom1 > 0 || comptonPhantom2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[0])
						event_scattered = true;
				}
				//}
				//if (event_true && scatter_components[1]) {
				else if (comptonCrystal1 > 0 || comptonCrystal2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[1])
						event_scattered = true;
				}
				//}
				//if (event_true && scatter_components[2]) {
				else if (RayleighPhantom1 > 0 || RayleighPhantom2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[2])
						event_scattered = true;
				}
				//}
				//if (event_true && scatter_components[3]) {
				else if (RayleighCrystal1 > 0 || RayleighCrystal2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[3])
						event_scattered = true;
				}
				//}
				if (outsize2 == 1ULL) {
					if (event_true && obtain_trues) {
						Ltrues[L1*detectors + L2] = Ltrues[L1*detectors + L2] + static_cast<uint16_t>(1);
						if (source)
							trues_loc[kk] = true;
					}
					else if (event_scattered && store_scatter) {
						Lscatter[L1*detectors + L2] = Lscatter[L1*detectors + L2] + static_cast<uint16_t>(1);
						if (source)
							scatter_loc[kk] = true;
					}
				}
				else {
					if (event_true && obtain_trues)
						Ltrues[kk] = 1u;
					else if (event_scattered && store_scatter)
						Lscatter[kk] = 1u;
				}
			}
			else if (!event_true && store_randoms) {
				if (outsize2 == 1ULL) {
					Lrandoms[L1*detectors + L2] = Lrandoms[L1*detectors + L2] + static_cast<uint16_t>(1);
					if (source)
						randoms_loc[kk] = true;
				}
				else
					Lrandoms[kk] = 1u;
			}
		}
		if (outsize2 == 1) {
			LL1[L1*detectors + L2] = LL1[L1*detectors + L2] + static_cast<uint16_t>(1);
		}
		else {
			LL1[kk] = static_cast<uint16_t>(L1 + 1);
			LL2[kk] = static_cast<uint16_t>(L2 + 1);
		}
		if (time2 >= aika) {
			tpoints[ll++] = kk;
			aika = time_intervals[++pa];
			//vali += vali;
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
	int_loc[1] = pa;
	tpoints[ll] = jj;
	if (begin) {
		int_loc[0] = -1;
		int_loc[1] = -1;
	}


	if (randoms_correction) {
		delay->SetBranchAddress("crystalID1", &crystalID1);
		delay->SetBranchAddress("crystalID2", &crystalID2);
		delay->SetBranchAddress("moduleID1", &moduleID1);
		delay->SetBranchAddress("moduleID2", &moduleID2);
		delay->SetBranchAddress("rsectorID1", &rsectorID1);
		delay->SetBranchAddress("rsectorID2", &rsectorID2);
		delay->SetBranchAddress("time2", &time2);

		nbytes = 0;
		bool begin = false;
		if (outsize2 > 1)
			begin = true;
		int pa = 0;
		double aika = 0.;
		int ll = 0, jj = 0;

		for (int64_t kk = 0; kk < Ndelays; kk++) {
			jj++;
			nbytes += delay->GetEntry(kk);

			if (time2 < alku)
				continue;
			else if (time2 > loppu)
				break;

			uint32_t ring_number1, ring_number2;
			if (linear_multip == 1) {
				ring_number1 = moduleID1;
				ring_number2 = moduleID2;
			}
			else {
				ring_number1 = (moduleID1 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
				ring_number2 = (moduleID2 % linear_multip)*cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
			}
			const uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring)*cryst_per_block + (crystalID1 % cryst_per_block);
			const uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring)*cryst_per_block + (crystalID2 % cryst_per_block);
			const uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
			const uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
			if (begin) {
				while (time2 >= time_intervals[pa])
					pa++;
				pa--;
				begin = false;
				aika = time_intervals[pa];
				int_loc_delay[0] = pa;
			}
			if (outsize2 == 1) {
				Ldelay1[L1*detectors + L2] = Ldelay1[L1*detectors + L2] + static_cast<uint16_t>(1);
			}
			else {
				Ldelay1[kk] = static_cast<uint16_t>(L1 + 1);
				Ldelay2[kk] = static_cast<uint16_t>(L2 + 1);
			}
			if (time2 >= aika) {
				tpoints_delay[ll++] = kk;
				aika = time_intervals[++pa];
			}
		}
		int_loc_delay[1] = pa;
		tpoints_delay[ll] = jj;
		if (begin) {
			int_loc_delay[0] = -1;
			int_loc_delay[1] = -1;
		}
	}
	return;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 16) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:invalidNumInputs",
			"16 input arguments required.");
	}
	else if (nlhs > 16) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:maxlhs",
			"Too many output arguments.");
	}

	/* Create a matrix for the return argument */
	//int * outsize2 = (int *)mxGetData(prhs[1]);
	double vali = (double)mxGetScalar(prhs[1]);
	double alku = (double)mxGetScalar(prhs[2]);
	double loppu = (double)mxGetScalar(prhs[3]);
	uint32_t detectors = (uint32_t)mxGetScalar(prhs[4]);
	uint32_t blocks_per_ring = (uint32_t)mxGetScalar(prhs[5]);
	uint32_t cryst_per_block = (uint32_t)mxGetScalar(prhs[6]);
	uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[7]);
	uint32_t linear_multp = (uint32_t)mxGetScalar(prhs[8]);
	bool source = (bool)mxGetScalar(prhs[9]);
	double* time_intervals = (double*)mxGetData(prhs[10]);
	bool obtain_trues = (bool)mxGetScalar(prhs[11]);
	bool store_scatter = (bool)mxGetScalar(prhs[12]);
	bool store_randoms = (bool)mxGetScalar(prhs[13]);
	bool* scatter_components = (bool*)mxGetData(prhs[14]);
	bool randoms_correction = (bool)mxGetScalar(prhs[15]);
	size_t outsize2 = (loppu - alku) / vali;

	// Count inputs and check for char type
	char *argv;
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	argv = mxArrayToString(prhs[0]);

	TChain *Coincidences = new TChain("Coincidences");
	Coincidences->Add(argv);

	TChain *delay;
	int64_t Ndelays = 0LL;

	if (randoms_correction) {
		delay = new TChain("delay");
		delay->Add(argv);
		Ndelays = delay->GetEntries();
	}


	int64_t Nentries = Coincidences->GetEntries();

	if (outsize2 == 1) {
		plhs[0] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (obtain_trues) {
			plhs[5] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		}
		else
			plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_randoms)
			plhs[6] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		else
			plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_scatter)
			plhs[7] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
		else
			plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (randoms_correction) {
			plhs[9] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
			plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
		else {
			plhs[9] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
	}
	else {
		plhs[0] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);

		if (obtain_trues)
			plhs[5] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_randoms)
			plhs[6] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (store_scatter)
			plhs[7] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		if (randoms_correction) {
			plhs[9] = mxCreateNumericMatrix(Ndelays, 1, mxUINT16_CLASS, mxREAL);
			plhs[10] = mxCreateNumericMatrix(Ndelays, 1, mxUINT16_CLASS, mxREAL);
		}
		else {
			plhs[9] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
	}
	plhs[2] = mxCreateNumericMatrix(outsize2 + 2, 1, mxUINT32_CLASS, mxREAL);
	if (randoms_correction)
		plhs[12] = mxCreateNumericMatrix(outsize2 + 2, 1, mxUINT32_CLASS, mxREAL);
	else
		plhs[12] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);

	if (source)
		plhs[3] = mxCreateNumericMatrix(Nentries, 6, mxSINGLE_CLASS, mxREAL);
	else
		plhs[3] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	plhs[4] = mxCreateNumericMatrix(2, 1, mxINT32_CLASS, mxREAL);

	plhs[11] = mxCreateNumericMatrix(2, 1, mxINT32_CLASS, mxREAL);

	/* Assign pointers to the various parameters */
	uint16_t * LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t * LL2 = (uint16_t*)mxGetData(plhs[1]);
	uint32_t * tpoints = (uint32_t*)mxGetData(plhs[2]);
	float* S = 0;
	float * output = 0;
	bool* trues_loc, *randoms_loc, *scatter_loc;
	if (source)
		S = (float*)mxGetData(plhs[3]);
	else
		output = (float*)mxGetData(plhs[3]);
	int *int_loc = (int*)mxGetData(plhs[4]);
	uint16_t * Ltrues = (uint16_t*)mxGetData(plhs[5]);
	uint16_t * Lscatter = (uint16_t*)mxGetData(plhs[6]);
	uint16_t * Lrandoms = (uint16_t*)mxGetData(plhs[7]);
	if (obtain_trues && source && outsize2 == 1) {
		plhs[8] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
		trues_loc = (bool*)mxGetData(plhs[8]);
	}
	else {
		trues_loc = 0;
		plhs[8] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
	}
	if (store_randoms && source && outsize2 == 1) {
		plhs[13] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
		randoms_loc = (bool*)mxGetData(plhs[13]);
	}
	else {
		randoms_loc = 0;
		plhs[13] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
	}
	if (store_scatter && source && outsize2 == 1) {
		plhs[14] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
		scatter_loc = (bool*)mxGetData(plhs[14]);
	}
	else {
		scatter_loc = 0;
		plhs[14] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
	}
	uint16_t * Ldelay1 = (uint16_t*)mxGetData(plhs[9]);
	uint16_t * Ldelay2 = (uint16_t*)mxGetData(plhs[10]);
	int *int_loc_delay = (int*)mxGetData(plhs[11]);
	uint32_t * tpoints_delay = (uint32_t*)mxGetData(plhs[12]);

	histogram(LL1, LL2, tpoints, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, 
		output, Coincidences, Nentries, time_intervals, int_loc, obtain_trues, store_scatter, store_randoms, scatter_components, Ltrues, Lscatter, 
		Lrandoms, trues_loc, Ndelays, randoms_correction, delay, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_loc, scatter_loc);


	delete Coincidences;
	if (randoms_correction)
		delete delay;
	mexEvalString("pause(.001);");
	gROOT->Reset();

	return;
}
