/**************************************************************************
* ROOT file import into MATLAB. This file contains the C implementation.
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
#include "mex.h"
#include <stdint.h>
#include <cmath>
#include <TROOT.h>
#include "TChain.h"


void histogram(uint16_t * LL1, uint16_t * LL2, uint32_t * tpoints, double vali, const double alku, const double loppu, const size_t outsize2, 
	const uint32_t detectors, bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring, 
	const uint32_t det_per_ring, float* S, float* output, TTree *Coincidences, const int64_t Nentries, const double* time_intervals, int* int_loc, 
	bool obtain_trues, bool store_scatter, bool store_randoms, const bool* scatter_components, uint16_t* Ltrues, uint16_t* Lscatter, 
	uint16_t* Lrandoms, bool* trues_loc, const int64_t Ndelays, bool randoms_correction, TTree *delay, uint16_t * Ldelay1, uint16_t * Ldelay2,
	int* int_loc_delay, uint32_t * tpoints_delay, bool* randoms_loc, bool* scatter_loc, float* x1, float* x2, float* y1, float* y2, float* z1, float* z2, 
	bool store_coordinates, const bool dynamic)
{

	Int_t crystalID1, crystalID2, moduleID1, moduleID2, rsectorID1, rsectorID2, eventID1, eventID2, comptonPhantom1 = 0, comptonPhantom2 = 0,
		comptonCrystal1 = 0, comptonCrystal2 = 0, RayleighPhantom1 = 0, RayleighPhantom2 = 0, RayleighCrystal1 = 0, RayleighCrystal2 = 0;
	Float_t sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
	Double_t time1 = alku, time2 = alku;

	int any = 0;
	int next = 0;
	bool no_time = false;

	if (Coincidences->GetBranchStatus("crystalID1"))
		Coincidences->SetBranchAddress("crystalID1", &crystalID1);
	else {
		mexPrintf("No crystal location information was found from file. Aborting.\n");
		return;
	}
	Coincidences->SetBranchAddress("crystalID2", &crystalID2);
	Coincidences->SetBranchAddress("moduleID1", &moduleID1);
	Coincidences->SetBranchAddress("moduleID2", &moduleID2);
	uint64_t summa = 0ULL;
	for (uint64_t kk = 0ULL; kk < Nentries; kk++) {
		Coincidences->GetEntry(kk);
		summa += moduleID1;
	}
	bool no_modules = false;
	if (summa == 0ULL)
		no_modules = true;
	Coincidences->SetBranchAddress("rsectorID1", &rsectorID1);
	Coincidences->SetBranchAddress("rsectorID2", &rsectorID2);
	if (source) {
		if (Coincidences->GetBranchStatus("sourcePosX1"))
			Coincidences->SetBranchAddress("sourcePosX1", &sourcePosX1);
		else {
			mexPrintf("No X-source coordinates saved for first photon, unable to save source coordinates\n");
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosX2"))
			Coincidences->SetBranchAddress("sourcePosX2", &sourcePosX1);
		else {
			mexPrintf("No X-source coordinates saved for second photon, unable to save source coordinates\n");
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosY1"))
			Coincidences->SetBranchAddress("sourcePosY1", &sourcePosY1);
		else {
			mexPrintf("No Y-source coordinates saved for first photon, unable to save source coordinates\n");
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosY2"))
			Coincidences->SetBranchAddress("sourcePosY2", &sourcePosY2);
		else {
			mexPrintf("No Y-source coordinates saved for second photon, unable to save source coordinates\n");
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosZ1"))
			Coincidences->SetBranchAddress("sourcePosZ1", &sourcePosZ1);
		else {
			mexPrintf("No Z-source coordinates saved for first photon, unable to save source coordinates\n");
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosZ2"))
			Coincidences->SetBranchAddress("sourcePosZ2", &sourcePosZ2);
		else {
			mexPrintf("No Z-source coordinates saved for second photon, unable to save source coordinates\n");
			source = false;
		}
	}
	//Coincidences->SetBranchAddress("time1", &time1);
	if (Coincidences->GetBranchStatus("time2"))
		Coincidences->SetBranchAddress("time2", &time2);
	else {
		if (dynamic) {
			mexPrintf("Dynamic examination selected, but no time information was found from file. Aborting.\n");
			return;
		}
		no_time = true;
	}
	if (store_coordinates) {
		if (Coincidences->GetBranchStatus("globalPosX1"))
			Coincidences->SetBranchAddress("globalPosX1", &globalPosX1);
		else {
			mexPrintf("No X-source coordinates saved for first photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosX2"))
			Coincidences->SetBranchAddress("globalPosX2", &globalPosX2);
		else {
			mexPrintf("No X-source coordinates saved for second photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosY1"))
			Coincidences->SetBranchAddress("globalPosY1", &globalPosY1);
		else {
			mexPrintf("No Y-source coordinates saved for first photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosY2"))
			Coincidences->SetBranchAddress("globalPosY2", &globalPosY2);
		else {
			mexPrintf("No Y-source coordinates saved for second photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosZ1"))
			Coincidences->SetBranchAddress("globalPosZ1", &globalPosZ1);
		else {
			mexPrintf("No Z-source coordinates saved for first photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosZ2"))
			Coincidences->SetBranchAddress("globalPosZ2", &globalPosZ2);
		else {
			mexPrintf("No Z-source coordinates saved for second photon interaction, unable to save interaction coordinates\n");
			store_coordinates = false;
		}
	}
	if (obtain_trues || store_scatter || store_randoms) {
		if (Coincidences->GetBranchStatus("eventID1"))
			Coincidences->SetBranchAddress("eventID1", &eventID1);
		else {
			mexPrintf("No event IDs saved for first photon, unable to save trues/scatter/randoms\n");
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
		if (Coincidences->GetBranchStatus("eventID2"))
			Coincidences->SetBranchAddress("eventID2", &eventID2);
		else {
			mexPrintf("No event IDs saved for second photon, unable to save trues/scatter/randoms\n");
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
	}
	if (obtain_trues || store_scatter) {
		//if (scatter_components[0]) {
		if (Coincidences->GetBranchStatus("comptonPhantom1"))
			Coincidences->SetBranchAddress("comptonPhantom1", &comptonPhantom1);
		else
			any++;
		if (Coincidences->GetBranchStatus("comptonPhantom2"))
			Coincidences->SetBranchAddress("comptonPhantom2", &comptonPhantom2);
		else
			any++;

		if (store_scatter && any == 2 && scatter_components[0] == 1) {
			mexPrintf("Compton phantom selected, but no scatter data was found from ROOT-file\n");
		}
		if (any == 2)
			next++;
		//}
		//if (scatter_components[1]) {
		if (Coincidences->GetBranchStatus("comptonCrystal1"))
			Coincidences->SetBranchAddress("comptonCrystal1", &comptonCrystal1);
		else
			any++;
		if (Coincidences->GetBranchStatus("comptonCrystal2"))
			Coincidences->SetBranchAddress("comptonCrystal2", &comptonCrystal2);
		else
			any++;

		if (store_scatter && ((any == 4 && next == 1) || (any == 2 && next == 0)) && scatter_components[1] == 1) {
			mexPrintf("Compton crystal selected, but no scatter data was found from ROOT-file\n");
		}

		if ((any == 4 && next == 1) || (any == 2 && next == 0))
			next++;
		//}
		//if (scatter_components[2]) {
		if (Coincidences->GetBranchStatus("RayleighPhantom1"))
			Coincidences->SetBranchAddress("RayleighPhantom1", &RayleighPhantom1);
		else
			any++;
		if (Coincidences->GetBranchStatus("RayleighPhantom2"))
			Coincidences->SetBranchAddress("RayleighPhantom2", &RayleighPhantom2);
		else
			any++;

		if (store_scatter && ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1)) && scatter_components[2] == 1) {
			mexPrintf("Rayleigh phantom selected, but no scatter data was found from ROOT-file\n");
		}

		if ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1))
			next++;
		//}
		//if (scatter_components[3]) {
		if (Coincidences->GetBranchStatus("RayleighCrystal1"))
			Coincidences->SetBranchAddress("RayleighCrystal1", &RayleighCrystal1);
		else
			any++;
		if (Coincidences->GetBranchStatus("RayleighCrystal2"))
			Coincidences->SetBranchAddress("RayleighCrystal2", &RayleighCrystal2);
		else
			any++;
		//}

		if (store_scatter && ((any == 8 && next == 3) || (any == 2 && next == 0) || (any == 4 && next == 1) || (any == 6 && next == 2)) && scatter_components[3] == 1) {
			mexPrintf("Rayleigh crystal selected, but no scatter data was found from ROOT-file\n");
		}


		if (store_scatter && any == 8) {
			mexPrintf("Store scatter selected, but no scatter data was found from ROOT-file\n");
		}

	}

	Int_t nbytes = 0;
	int ll = 0, jj = -1;
	bool begin = false;
	if (outsize2 > 1)
		begin = true;
	int pa = 0;
	double aika = 0.;
	int_loc[0] = 1;
	tpoints[ll] = 0;

	for (uint64_t kk = 0; kk < Nentries; kk++) {

		jj++;
		nbytes += Coincidences->GetEntry(kk);

		if (time2 < alku)
			continue;
		else if (time2 > loppu) {
			//int_loc[1] = pa;
			//tpoints[ll] = jj;
			//if (begin) {
			//	int_loc[0] = -1;
			//	int_loc[1] = -1;
			//}
			break;
		}
		uint32_t ring_number1, ring_number2;

		if (linear_multip == 1 && !no_modules) {
			ring_number1 = moduleID1;
			ring_number2 = moduleID2;
		}
		else if (linear_multip == 1 && no_modules) {
			ring_number1 = rsectorID1;
			ring_number2 = rsectorID2;
		}
		else if (no_modules) {
			ring_number1 = static_cast<uint32_t>(std::floor(static_cast<double>(rsectorID1) / static_cast<double>(blocks_per_ring))) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
			ring_number2 = static_cast<uint32_t>(std::floor(static_cast<double>(rsectorID2) / static_cast<double>(blocks_per_ring))) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
		}
		else {
			ring_number1 = (moduleID1 % linear_multip) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
			ring_number2 = (moduleID2 % linear_multip) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
		}
		const uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring) * cryst_per_block + (crystalID1 % cryst_per_block);
		const uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring) * cryst_per_block + (crystalID2 % cryst_per_block);
		uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
		uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
		if (L2 > L1) {
			const uint32_t L3 = L1;
			L1 = L2;
			L2 = L3;
		}
		if (begin) {
			while (time2 >= time_intervals[pa])
				pa++;
			begin = false;
			tpoints[ll++] = 0u;
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
						Ltrues[L1 * detectors + L2] = Ltrues[L1 * detectors + L2] + static_cast<uint16_t>(1);
						if (source)
							trues_loc[kk] = true;
					}
					else if (event_scattered && store_scatter) {
						Lscatter[L1 * detectors + L2] = Lscatter[L1 * detectors + L2] + static_cast<uint16_t>(1);
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
					Lrandoms[L1 * detectors + L2] = Lrandoms[L1 * detectors + L2] + static_cast<uint16_t>(1);
					if (source)
						randoms_loc[kk] = true;
				}
				else
					Lrandoms[kk] = 1u;
			}
		}
		if (outsize2 == 1ULL) {
			LL1[L1 * detectors + L2] = LL1[L1 * detectors + L2] + static_cast<uint16_t>(1);
		}
		else {
			LL1[kk] = static_cast<uint16_t>(L1 + 1);
			LL2[kk] = static_cast<uint16_t>(L2 + 1);
		}
		if (outsize2 > 1ULL && time2 >= aika) {
			tpoints[ll++] = jj;
			aika = time_intervals[++pa];
		}
		if (source) {
			S[kk] = sourcePosX1;
			S[kk + Nentries] = sourcePosY1;
			S[kk + Nentries * 2] = sourcePosZ1;
			S[kk + Nentries * 3] = sourcePosX2;
			S[kk + Nentries * 4] = sourcePosY2;
			S[kk + Nentries * 5] = sourcePosZ2;
		}
		if (store_coordinates) {
			x1[kk] = globalPosX1;
			x2[kk] = globalPosX2;
			y1[kk] = globalPosY1;
			y2[kk] = globalPosY2;
			z1[kk] = globalPosZ1;
			z2[kk] = globalPosZ2;
		}
	}
	if (pa == 0)
		pa++;
	if (ll == 0)
		ll++;
	int_loc[1] = pa;
	tpoints[ll] = jj;
	if (begin) {
		int_loc[0] = 0;
		int_loc[1] = 0;
	}


	if (randoms_correction) {
		delay->SetBranchAddress("crystalID1", &crystalID1);
		delay->SetBranchAddress("crystalID2", &crystalID2);
		delay->SetBranchAddress("moduleID1", &moduleID1);
		delay->SetBranchAddress("moduleID2", &moduleID2);
		delay->SetBranchAddress("rsectorID1", &rsectorID1);
		delay->SetBranchAddress("rsectorID2", &rsectorID2);
		if (delay->GetBranchStatus("time2"))
			delay->SetBranchAddress("time2", &time2);

		nbytes = 0;
		bool begin = false;
		if (outsize2 > 1)
			begin = true;
		int pa = 0;
		double aika = 0.;
		int ll = 0, jj = -1;
		tpoints_delay[ll] = 0;

		for (int64_t kk = 0; kk < Ndelays; kk++) {
			jj++;
			nbytes += delay->GetEntry(kk);

			if (time2 < alku)
				continue;
			else if (time2 > loppu) {
				//int_loc_delay[1] = pa;
				//tpoints_delay[ll] = jj;
				//if (begin) {
				//	int_loc_delay[0] = -1;
				//	int_loc_delay[1] = -1;
				//}
				break;
			}

			uint32_t ring_number1, ring_number2;
			if (linear_multip == 1 && !no_modules) {
				ring_number1 = moduleID1;
				ring_number2 = moduleID2;
			}
			else if (linear_multip == 1 && no_modules) {
				ring_number1 = rsectorID1;
				ring_number2 = rsectorID2;
			}
			else if (no_modules) {
				ring_number1 = static_cast<uint32_t>(std::floor(static_cast<double>(rsectorID1) / static_cast<double>(blocks_per_ring))) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
				ring_number2 = static_cast<uint32_t>(std::floor(static_cast<double>(rsectorID2) / static_cast<double>(blocks_per_ring))) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
			}
			else {
				ring_number1 = (moduleID1 % linear_multip) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID1) / static_cast<double>(cryst_per_block)));
				ring_number2 = (moduleID2 % linear_multip) * cryst_per_block + static_cast<uint32_t>(std::floor(static_cast<double>(crystalID2) / static_cast<double>(cryst_per_block)));
			}
			const uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring)*cryst_per_block + (crystalID1 % cryst_per_block);
			const uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring)*cryst_per_block + (crystalID2 % cryst_per_block);
			uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
			uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
			if (L2 > L1) {
				const uint32_t L3 = L1;
				L1 = L2;
				L2 = L3;
			}
			if (begin) {
				while (time2 >= time_intervals[pa])
					pa++;
				begin = false;
				tpoints_delay[ll++] = 0u;
				aika = time_intervals[pa];
				int_loc_delay[0] = pa;
			}
			if (outsize2 == 1ULL) {
				Ldelay1[L1*detectors + L2] = Ldelay1[L1*detectors + L2] + static_cast<uint16_t>(1);
			}
			else {
				Ldelay1[kk] = static_cast<uint16_t>(L1 + 1);
				Ldelay2[kk] = static_cast<uint16_t>(L2 + 1);
			}
			if (time2 >= aika && outsize2 > 1ULL) {
				tpoints_delay[ll++] = jj;
				aika = time_intervals[++pa];
			}
		}
		if (pa == 0)
			pa++;
		if (ll == 0)
			ll++;
		int_loc_delay[1] = pa;
		tpoints_delay[ll] = jj;
		if (begin) {
			int_loc_delay[0] = 0;
			int_loc_delay[1] = 0;
		}
	}
	return;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 17) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:invalidNumInputs",
			"17 input arguments required.");
	}
	else if (nlhs > 21) {
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
	bool store_coordinates = (bool)mxGetScalar(prhs[16]);
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
	bool dynamic = false;

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
		dynamic = true;
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

	if (store_coordinates) {
		plhs[15] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		plhs[16] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		plhs[17] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		plhs[18] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		plhs[19] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		plhs[20] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
	}
	else {
		plhs[15] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[16] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[17] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[18] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[19] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[20] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	}

	/* Assign pointers to the various parameters */
	uint16_t * LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t * LL2 = (uint16_t*)mxGetData(plhs[1]);
	uint32_t * tpoints = (uint32_t*)mxGetData(plhs[2]);
	float* S = nullptr;
	float * output = nullptr;
	float* x1 = nullptr, * x2 = nullptr, * y1 = nullptr, * y2 = nullptr, * z1 = nullptr, * z2 = nullptr;
	bool* trues_loc, *randoms_loc, *scatter_loc;
	if (source)
		S = (float*)mxGetData(plhs[3]);
	else
		output = (float*)mxGetData(plhs[3]);
	int *int_loc = (int*)mxGetData(plhs[4]);
	uint16_t* Ltrues = (uint16_t*)mxGetData(plhs[5]);
	uint16_t* Lrandoms = (uint16_t*)mxGetData(plhs[6]);
	uint16_t* Lscatter = (uint16_t*)mxGetData(plhs[7]);
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
	if (store_coordinates) {
		x1 = (float*)mxGetData(plhs[15]);
		x2 = (float*)mxGetData(plhs[16]);
		y1 = (float*)mxGetData(plhs[17]);
		y2 = (float*)mxGetData(plhs[18]);
		z1 = (float*)mxGetData(plhs[19]);
		z2 = (float*)mxGetData(plhs[20]);
	}
	uint16_t * Ldelay1 = (uint16_t*)mxGetData(plhs[9]);
	uint16_t * Ldelay2 = (uint16_t*)mxGetData(plhs[10]);
	int *int_loc_delay = (int*)mxGetData(plhs[11]);
	uint32_t * tpoints_delay = (uint32_t*)mxGetData(plhs[12]);

	histogram(LL1, LL2, tpoints, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, 
		output, Coincidences, Nentries, time_intervals, int_loc, obtain_trues, store_scatter, store_randoms, scatter_components, Ltrues, Lscatter, 
		Lrandoms, trues_loc, Ndelays, randoms_correction, delay, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_loc, scatter_loc, 
		x1, x2, y1, y2, z1, z2, store_coordinates, dynamic);


	delete Coincidences;
	if (randoms_correction)
		delete delay;
	mexEvalString("pause(.001);");
	gROOT->Reset();

	return;
}
