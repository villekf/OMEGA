#include <octave/oct.h>
#include <TROOT.h>
#include "TChain.h"
#include "saveSinogram.h"


void histogram(octave_uint16* LL1, octave_uint16* LL2, octave_uint32* tpoints, double vali, const double alku, const double loppu, const size_t outsize2,
	const uint32_t detectors, bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring,
	const uint32_t det_per_ring, float* S, TTree* Coincidences, const int64_t Nentries, const double* time_intervals, octave_int32* int_loc,
	bool obtain_trues, bool store_scatter, bool store_randoms, const bool* scatter_components, octave_uint16* Ltrues, octave_uint16* Lscatter,
	octave_uint16* Lrandoms, bool* trues_loc, const int64_t Ndelays, bool randoms_correction, TTree* delay, octave_uint16* Ldelay1, octave_uint16* Ldelay2,
	octave_int32* int_loc_delay, octave_uint32* tpoints_delay, bool* randoms_loc, bool* scatter_loc, float* x1, float* x2, float* y1, float* y2, 
	float* z1, float* z2, bool store_coordinates, const bool dynamic, const uint32_t cryst_per_block_z, const uint32_t transaxial_multip, const uint32_t rings, 
	const uint64_t sinoSize, const uint32_t Ndist, const uint32_t Nang, const uint32_t ringDifference, const uint32_t span, const octave_uint32* seg,
	const uint64_t NT, const uint64_t TOFSize, const int32_t nDistSide, const bool storeRawData, octave_uint16* Sino, octave_uint16* SinoT, octave_uint16* SinoC, 
	octave_uint16* SinoR, octave_uint16* SinoD, const int32_t detWPseudo, const int32_t nPseudos, const double binSize)
{

	Int_t crystalID1 = 0, crystalID2 = 0, moduleID1 = 0, moduleID2 = 0, submoduleID1 = 0, submoduleID2 = 0, rsectorID1, rsectorID2, eventID1, eventID2, comptonPhantom1 = 0, comptonPhantom2 = 0,
		comptonCrystal1 = 0, comptonCrystal2 = 0, RayleighPhantom1 = 0, RayleighPhantom2 = 0, RayleighCrystal1 = 0, RayleighCrystal2 = 0;
	Float_t sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
	Double_t time1 = alku, time2 = alku;

	int any = 0;
	int next = 0;
	bool no_time = false;

	if (Coincidences->GetBranchStatus("crystalID1"))
		Coincidences->SetBranchAddress("crystalID1", &crystalID1);
	else {
		octave_stdout << "No crystal location information was found from file. Aborting.\n";
		return;
	}
	if (Coincidences->GetBranchStatus("crystalID2"))
		Coincidences->SetBranchAddress("crystalID2", &crystalID2);
	else {
		octave_stdout << "No crystal location information was found from file. Aborting.\n";
		return;
	}
	if (Coincidences->GetBranchStatus("moduleID1"))
		Coincidences->SetBranchAddress("moduleID1", &moduleID1);
	if (Coincidences->GetBranchStatus("moduleID2"))
		Coincidences->SetBranchAddress("moduleID2", &moduleID2);
	if (Coincidences->GetBranchStatus("moduleID1"))
		Coincidences->SetBranchAddress("submoduleID1", &submoduleID1);
	if (Coincidences->GetBranchStatus("submoduleID2"))
		Coincidences->SetBranchAddress("submoduleID2", &submoduleID2);
	uint64_t summa = 0ULL;
	uint64_t summaS = 0ULL;
	for (uint64_t kk = 0ULL; kk < Nentries; kk++) {
		Coincidences->GetEntry(kk);
		if (summa == 0ULL)
			summa += moduleID1;
		if (summaS == 0ULL)
			summaS += submoduleID1;
		if (summa > 0 && summaS > 0)
			break;
	}
	bool no_modules = false;
	bool no_submodules = true;
	if (summa == 0ULL)
		no_modules = true;
	if (summaS > 0ULL)
		no_submodules = false;
	Coincidences->SetBranchAddress("rsectorID1", &rsectorID1);
	Coincidences->SetBranchAddress("rsectorID2", &rsectorID2);
	if (source) {
		if (Coincidences->GetBranchStatus("sourcePosX1"))
			Coincidences->SetBranchAddress("sourcePosX1", &sourcePosX1);
		else {
			octave_stdout << "No X-source coordinates saved for first photon, unable to save source coordinates\n";
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosX2"))
			Coincidences->SetBranchAddress("sourcePosX2", &sourcePosX2);
		else {
			octave_stdout << "No X-source coordinates saved for second photon, unable to save source coordinates\n";
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosY1"))
			Coincidences->SetBranchAddress("sourcePosY1", &sourcePosY1);
		else {
			octave_stdout << "No Y-source coordinates saved for first photon, unable to save source coordinates\n";
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosY2"))
			Coincidences->SetBranchAddress("sourcePosY2", &sourcePosY2);
		else {
			octave_stdout << "No Y-source coordinates saved for second photon, unable to save source coordinates\n";
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosZ1"))
			Coincidences->SetBranchAddress("sourcePosZ1", &sourcePosZ1);
		else {
			octave_stdout << "No Z-source coordinates saved for first photon, unable to save source coordinates\n";
			source = false;
		}
		if (Coincidences->GetBranchStatus("sourcePosZ2"))
			Coincidences->SetBranchAddress("sourcePosZ2", &sourcePosZ2);
		else {
			octave_stdout << "No Z-source coordinates saved for second photon, unable to save source coordinates\n";
			source = false;
		}
	}
	if (Coincidences->GetBranchStatus("time1"))
		Coincidences->SetBranchAddress("time1", &time1);
	if (Coincidences->GetBranchStatus("time2"))
		Coincidences->SetBranchAddress("time2", &time2);
	else {
		if (dynamic) {
			octave_stdout << "Dynamic examination selected, but no time information was found from file. Aborting.\n";
			return;
		}
		no_time = true;
	}
	if (store_coordinates) {
		if (Coincidences->GetBranchStatus("globalPosX1"))
			Coincidences->SetBranchAddress("globalPosX1", &globalPosX1);
		else {
			octave_stdout << "No X-source coordinates saved for first photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosX2"))
			Coincidences->SetBranchAddress("globalPosX2", &globalPosX2);
		else {
			octave_stdout << "No X-source coordinates saved for second photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosY1"))
			Coincidences->SetBranchAddress("globalPosY1", &globalPosY1);
		else {
			octave_stdout << "No Y-source coordinates saved for first photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosY2"))
			Coincidences->SetBranchAddress("globalPosY2", &globalPosY2);
		else {
			octave_stdout << "No Y-source coordinates saved for second photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosZ1"))
			Coincidences->SetBranchAddress("globalPosZ1", &globalPosZ1);
		else {
			octave_stdout << "No Z-source coordinates saved for first photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
		if (Coincidences->GetBranchStatus("globalPosZ2"))
			Coincidences->SetBranchAddress("globalPosZ2", &globalPosZ2);
		else {
			octave_stdout << "No Z-source coordinates saved for second photon interaction, unable to save interaction coordinates\n";
			store_coordinates = false;
		}
	}
	if (obtain_trues || store_scatter || store_randoms) {
		if (Coincidences->GetBranchStatus("eventID1"))
			Coincidences->SetBranchAddress("eventID1", &eventID1);
		else {
			octave_stdout << "No event IDs saved for first photon, unable to save trues/scatter/randoms\n";
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
		if (Coincidences->GetBranchStatus("eventID2"))
			Coincidences->SetBranchAddress("eventID2", &eventID2);
		else {
			octave_stdout << "No event IDs saved for second photon, unable to save trues/scatter/randoms\n";
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
	}
	if (obtain_trues || store_scatter) {
		if (Coincidences->GetBranchStatus("comptonPhantom1"))
			Coincidences->SetBranchAddress("comptonPhantom1", &comptonPhantom1);
		else
			any++;
		if (Coincidences->GetBranchStatus("comptonPhantom2"))
			Coincidences->SetBranchAddress("comptonPhantom2", &comptonPhantom2);
		else
			any++;

		if (store_scatter && any == 2 && scatter_components[0]) {
			octave_stdout << "Compton phantom selected, but no scatter data was found from ROOT-file\n";
		}
		if (Coincidences->GetBranchStatus("comptonCrystal1"))
			Coincidences->SetBranchAddress("comptonCrystal1", &comptonCrystal1);
		else
			any++;
		if (Coincidences->GetBranchStatus("comptonCrystal2"))
			Coincidences->SetBranchAddress("comptonCrystal2", &comptonCrystal2);
		else
			any++;

		if (store_scatter && any == 2 && scatter_components[0]) {
			octave_stdout << "Compton crystal selected, but no scatter data was found from ROOT-file\n";
		}
		if (Coincidences->GetBranchStatus("RayleighPhantom1"))
			Coincidences->SetBranchAddress("RayleighPhantom1", &RayleighPhantom1);
		else
			any++;
		if (Coincidences->GetBranchStatus("RayleighPhantom2"))
			Coincidences->SetBranchAddress("RayleighPhantom2", &RayleighPhantom2);
		else
			any++;

		if (store_scatter && any == 2 && scatter_components[0]) {
			octave_stdout << "Rayleigh phantom selected, but no scatter data was found from ROOT-file\n";
		}
		if (Coincidences->GetBranchStatus("RayleighCrystal1"))
			Coincidences->SetBranchAddress("RayleighCrystal1", &RayleighCrystal1);
		else
			any++;
		if (Coincidences->GetBranchStatus("RayleighCrystal2"))
			Coincidences->SetBranchAddress("RayleighCrystal2", &RayleighCrystal2);
		else
			any++;

		if (store_scatter && any == 2 && scatter_components[0]) {
			octave_stdout << "Rayleigh crystal selected, but no scatter data was found from ROOT-file\n";
		}

		if (store_scatter && any == 8) {
			octave_stdout << "Store scatter selected, but no scatter data was found from ROOT-file\n";
		}
	}
	const bool pseudoD = detWPseudo > det_per_ring;
	const bool pseudoR = nPseudos > 0;
	int32_t gapSize = 0;
	if (pseudoR) {
		gapSize = rings / (nPseudos + 1);
	}

	Int_t nbytes = 0;
	int ll = 0, jj = 0;
	bool begin = false;
	if (outsize2 > 1)
		begin = true;
	int pa = 0;
	double aika = 0.;
	int_loc[0] = 1;
	tpoints[ll] = 0;
	uint32_t L1 = 0U;
	uint32_t L2 = 0U;
	uint64_t nBins = 0ULL;
	if (TOFSize > sinoSize) {
		nBins = TOFSize / sinoSize;
	}

	for (int64_t kk = 0; kk < Nentries; kk++) {

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
		bool event_true = true;
		bool event_scattered = false;
		if (obtain_trues || store_scatter || store_randoms) {
			if (eventID1 != eventID2) {
				event_true = false;
			}
			if (event_true && (obtain_trues || store_scatter)) {
				if (comptonPhantom1 > 0 || comptonPhantom2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[0])
						event_scattered = true;
				}
				else if (comptonCrystal1 > 0 || comptonCrystal2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[1])
						event_scattered = true;
				}
				else if (RayleighPhantom1 > 0 || RayleighPhantom2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[2])
						event_scattered = true;
				}
				else if (RayleighCrystal1 > 0 || RayleighCrystal2 > 0) {
					event_true = false;
					if (store_scatter && scatter_components[3])
						event_scattered = true;
				}
			}
		}
		if (begin) {
			while (time2 >= time_intervals[pa])
				pa++;
			pa--;
			begin = false;
			aika = time_intervals[pa];
			int_loc[0] = pa;
		}
		const double time = time2;
		uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
		detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
			submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block, cryst_per_block_z, transaxial_multip, rings);
		uint64_t bins = 0;
		if (TOFSize > sinoSize) {
			double timeDif = (time2 - time1) / 2.;
			if (ring_pos2 > ring_pos1)
				timeDif = -timeDif;
			bins = static_cast<uint64_t>(std::floor((std::abs(timeDif) + binSize / 2.) / binSize));
			if (timeDif < 0)
				bins *= 2ULL;
			else if (bins > 2)
				bins = bins * 2ULL - 1ULL;
		}
		if (storeRawData) {
			L1 = ring_number1 * det_per_ring + ring_pos1;
			L2 = ring_number2 * det_per_ring + ring_pos2;
			if (L2 > L1) {
				const uint32_t L3 = L1;
				L1 = L2;
				L2 = L3;
			}
			if (obtain_trues || store_scatter || store_randoms) {
				if ((event_true && obtain_trues) || (event_scattered && store_scatter)) {
					if (outsize2 == 1ULL) {
						if (event_true && obtain_trues) {
							Ltrues[L1 * detectors + L2] = Ltrues[L1 * detectors + L2] + static_cast<octave_uint16>(1);
						}
						else if (event_scattered && store_scatter) {
							Lscatter[L1 * detectors + L2] = Lscatter[L1 * detectors + L2] + static_cast<octave_uint16>(1);
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
						Lrandoms[L1 * detectors + L2] = Lrandoms[L1 * detectors + L2] + static_cast<octave_uint16>(1);
					}
					else
						Lrandoms[kk] = 1u;
				}
			}
			if (outsize2 == 1ULL) {
				LL1[L1 * detectors + L2] = LL1[L1 * detectors + L2] + static_cast<octave_uint16>(1);
			}
			else {
				LL1[kk] = static_cast<octave_uint16>(L1 + 1);
				LL2[kk] = static_cast<octave_uint16>(L2 + 1);
			}
		}
		if (pseudoD) {
			ring_pos1 += ring_pos1 / cryst_per_block;
			ring_pos2 += ring_pos2 / cryst_per_block;
		}
		if (pseudoR) {
			ring_number1 += ring_number1 / gapSize;
			ring_number2 += ring_number2 / gapSize;
		}
		const int64_t sinoIndex = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ringDifference, span, seg, time, NT, TOFSize, 
			vali, alku, detWPseudo, rings, bins, nDistSide);
		if (sinoIndex >= 0) {
			Sino[sinoIndex] = Sino[sinoIndex] + static_cast<octave_uint16>(1);
			if (event_true && obtain_trues)
				SinoT[sinoIndex] = SinoT[sinoIndex] + static_cast<octave_uint16>(1);
			else if (event_scattered && store_scatter)
				SinoC[sinoIndex] = SinoC[sinoIndex] + static_cast<octave_uint16>(1);
			else if (!event_true && store_randoms)
				SinoR[sinoIndex] = SinoR[sinoIndex] + static_cast<octave_uint16>(1);
		}
		if (time2 >= aika && outsize2 > 1ULL) {
			tpoints[ll++] = kk;
			aika = time_intervals[++pa];
			//vali += vali;
		}
		if (source) {
			if (outsize2 == 1ULL || !storeRawData) {
				if (event_true && obtain_trues)
					trues_loc[kk] = true;
				else if (event_scattered && store_scatter)
					scatter_loc[kk] = true;
				else if (!event_true && !event_scattered && store_randoms)
					randoms_loc[kk] = true;
			}
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
	tpoints[ll] = static_cast<uint32_t>(jj);
	if (begin) {
		int_loc[0] = 0;
		int_loc[1] = 0;
	}


	if (randoms_correction) {
		delay->SetBranchAddress("crystalID1", &crystalID1);
		delay->SetBranchAddress("crystalID2", &crystalID2);
		delay->SetBranchAddress("moduleID1", &moduleID1);
		delay->SetBranchAddress("moduleID2", &moduleID2);
		delay->SetBranchAddress("submoduleID1", &submoduleID1);
		delay->SetBranchAddress("submoduleID2", &submoduleID2);
		delay->SetBranchAddress("rsectorID1", &rsectorID1);
		delay->SetBranchAddress("rsectorID2", &rsectorID2);
		if (delay->GetBranchStatus("time2"))
			delay->SetBranchAddress("time2", &time2);

		nbytes = 0;
		uint64_t bins = 0;
		bool begin = false;
		if (outsize2 > 1)
			begin = true;
		int pa = 0;
		double aika = 0.;
		int ll = 0, jj = 0;
		int_loc_delay[0] = 0;
		tpoints_delay[ll] = 0;

		for (int64_t kk = 0; kk < Ndelays; kk++) {
			jj++;
			nbytes += delay->GetEntry(kk);

			if (time2 < alku)
				continue;
			else if (time2 > loppu) {
				break;
			}
			if (begin) {
				while (time2 >= time_intervals[pa])
					pa++;
				//pa--;
				begin = false;
				tpoints_delay[ll++] = 0u;
				aika = time_intervals[pa];
				int_loc_delay[0] = pa;
			}
			if (time2 >= aika && outsize2 > 1ULL) {
				tpoints_delay[ll++] = kk;
				aika = time_intervals[++pa];
			}

			const double time = time2;
			uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
			detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
				submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block, cryst_per_block_z, transaxial_multip, rings);
			if (storeRawData) {
				L1 = ring_number1 * det_per_ring + ring_pos1;
				L2 = ring_number2 * det_per_ring + ring_pos2;
				if (L2 > L1) {
					const uint32_t L3 = L1;
					L1 = L2;
					L2 = L3;
				}
				if (outsize2 == 1ULL) {
					Ldelay1[L1 * detectors + L2] = Ldelay1[L1 * detectors + L2] + static_cast<octave_uint16>(1);
				}
				else {
					Ldelay1[kk] = static_cast<uint16_t>(L1 + 1);
					Ldelay2[kk] = static_cast<uint16_t>(L2 + 1);
				}
			}
			if (pseudoD) {
				ring_pos1 += ring_pos1 / cryst_per_block;
				ring_pos2 += ring_pos2 / cryst_per_block;
			}
			if (pseudoR) {
				ring_number1 += ring_number1 / gapSize;
				ring_number2 += ring_number2 / gapSize;
			}
			const int64_t sinoIndex = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ringDifference, span, seg, time, NT, sinoSize,
				vali, alku, detWPseudo, rings, bins, nDistSide);
			if (sinoIndex >= 0) {
				SinoD[sinoIndex] = SinoD[sinoIndex] + static_cast<octave_uint16>(1);
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


DEFUN_DLD(GATE_root_matlab_oct, prhs, nargout, "GATE ROOT help") {


	double vali = prhs(1).scalar_value();
	double alku = prhs(2).scalar_value();
	double loppu = prhs(3).scalar_value();
	uint32_t detectors = prhs(4).uint32_scalar_value();
	uint32_t blocks_per_ring = prhs(5).uint32_scalar_value();
	uint32_t cryst_per_block = prhs(6).uint32_scalar_value();
	uint32_t det_per_ring = prhs(7).uint32_scalar_value();
	uint32_t linear_multp = prhs(8).uint32_scalar_value();
	bool source = prhs(9).bool_value();
	NDArray time_intervals = prhs(10).array_value();
	bool obtain_trues = prhs(11).bool_value();
	bool store_scatter = prhs(12).bool_value();
	bool store_randoms = prhs(13).bool_value();
	boolNDArray scatter_components = prhs(14).bool_array_value();
	bool randoms_correction = prhs(15).bool_value();
	bool store_coordinates = prhs(16).bool_value();
	uint32_t cryst_per_block_z = prhs(17).uint32_scalar_value();
	uint32_t transaxial_multip = prhs(18).uint32_scalar_value();
	uint32_t rings = prhs(19).uint32_scalar_value();
	uint64_t sinoSize = prhs(20).uint64_scalar_value();
	uint32_t Ndist = prhs(21).uint32_scalar_value();
	uint32_t Nang = prhs(22).uint32_scalar_value();
	uint32_t ringDifference = prhs(23).uint32_scalar_value();
	uint32_t span = prhs(24).uint32_scalar_value();
	uint32NDArray seg = prhs(25).uint32_array_value();
	uint64_t NT = prhs(26).uint64_scalar_value();
	uint64_t TOFSize = prhs(27).uint64_scalar_value();
	int32_t nDistSide = prhs(28).int32_scalar_value();
	bool storeRawData = prhs(29).bool_value();
	uint16NDArray SinoO = prhs(30).uint16_array_value();
	uint16NDArray SinoOT = prhs(31).uint16_array_value();
	uint16NDArray SinoOC = prhs(32).uint16_array_value();
	uint16NDArray SinoOR = prhs(33).uint16_array_value();
	uint16NDArray SinoOD = prhs(34).uint16_array_value();
	const int32_t detWPseudo = prhs(35).int32_scalar_value();
	const int32_t nPseudos = prhs(36).int32_scalar_value();
	double binSize = prhs(37).scalar_value();
	size_t outsize2 = (loppu - alku) / vali;
	const bool* scatter_components_p = scatter_components.fortran_vec();
	const double* time_intervals_p = time_intervals.fortran_vec();
	const octave_uint32* seg_p = seg.fortran_vec();

	// Count inputs and check for char type
	//const char *argv;
	//if (!mxIsChar(prhs[0]))
	//	mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	//argv = mxArrayToString(prhs[0]);
	//charNDArray apu = prhs(0).char_array_value();
	charMatrix apu = prhs(0).char_matrix_value();
	std::string tmp = apu.row_as_string(0);
	//const char* argv = prhs(0).string_value().c_str();
	//const char* argv = apu.fortran_vec();

	//octave_stdout << argv << std::endl;
	//octave_stdout << tmp.c_str() << std::endl;

	bool dynamic = false;

	TChain *Coincidences = new TChain("Coincidences");
	Coincidences->Add(tmp.c_str());

	TChain *delay;
	int64_t Ndelays = 0LL;

	if (randoms_correction) {
		delay = new TChain("delay");
		delay->Add(tmp.c_str());
		Ndelays = delay->GetEntries();
	}


	int64_t Nentries = Coincidences->GetEntries();

	/* Assign pointers to the various parameters */
	uint16NDArray LL1;
	uint16NDArray LL2;
	uint16NDArray Ltrues;
	uint16NDArray Lrandoms;
	uint16NDArray Lscatter;
	uint16NDArray Ldelay1;
	uint16NDArray Ldelay2;
	FloatNDArray x1;
	FloatNDArray x2;
	FloatNDArray y1;
	FloatNDArray y2;
	FloatNDArray z1;
	FloatNDArray z2;
	if (outsize2 == 1) {
		if (storeRawData) {
			LL1.resize(dim_vector(detectors, detectors));
			LL2.resize(dim_vector(1, 1));
		}
		else {
			LL1.resize(dim_vector(1, 1));
			LL2.resize(dim_vector(1, 1));
		}
		if (obtain_trues && storeRawData) {
			Ltrues.resize(dim_vector(detectors, detectors));
		}
		else
			Ltrues.resize(dim_vector(1, 1));
		if (store_randoms && storeRawData)
			Lrandoms.resize(dim_vector(detectors, detectors));
		else
			Lrandoms.resize(dim_vector(1, 1));
		if (store_scatter && storeRawData)
			Lscatter.resize(dim_vector(detectors, detectors));
		else
			Lscatter.resize(dim_vector(1, 1));
		if (randoms_correction && storeRawData) {
			Ldelay1.resize(dim_vector(detectors, detectors));
			Ldelay2.resize(dim_vector(1, 1));
		}
		else {
			Ldelay1.resize(dim_vector(1, 1));
			Ldelay2.resize(dim_vector(1, 1));
		}
	}
	else {
		if (storeRawData) {
			LL1.resize(dim_vector(Nentries, 1));
			LL2.resize(dim_vector(Nentries, 1));
		}
		else {
			LL1.resize(dim_vector(1, 1));
			LL2.resize(dim_vector(1, 1));
		}
		if (obtain_trues && storeRawData) {
			Ltrues.resize(dim_vector(Nentries, 1));
		}
		else
			Ltrues.resize(dim_vector(1, 1));
		if (store_randoms && storeRawData)
			Lrandoms.resize(dim_vector(Nentries, 1));
		else
			Lrandoms.resize(dim_vector(1, 1));
		if (store_scatter && storeRawData)
			Lscatter.resize(dim_vector(Nentries, 1));
		else
			Lscatter.resize(dim_vector(1, 1));
		if (randoms_correction && storeRawData) {
			Ldelay1.resize(dim_vector(Ndelays, 1));
			Ldelay2.resize(dim_vector(Ndelays, 1));
		}
		else {
			Ldelay1.resize(dim_vector(1, 1));
			Ldelay2.resize(dim_vector(1, 1));
		}
		dynamic = true;
	}
	if (store_coordinates) {
		x1.resize(dim_vector(Nentries, 1));
		x2.resize(dim_vector(Nentries, 1));
		y1.resize(dim_vector(Nentries, 1));
		y2.resize(dim_vector(Nentries, 1));
		z1.resize(dim_vector(Nentries, 1));
		z2.resize(dim_vector(Nentries, 1));
	}
	uint32NDArray tpoints, tpoints_delay;
	FloatNDArray S;
	int32NDArray int_loc, int_loc_delay;
	boolNDArray trues_loc, randoms_loc, scatter_loc;
	if (randoms_correction)
		tpoints_delay.resize(dim_vector(outsize2 + 2, 1));
	else
		tpoints_delay.resize(dim_vector(1, 1));

	if (source)
		S.resize(dim_vector(Nentries, 6));
	else
		S.resize(dim_vector(1, 1));
	int_loc.resize(dim_vector(2, 1));
	tpoints.resize(dim_vector(outsize2 + 2, 1));;

	if (obtain_trues && source && outsize2 == 1) {
		trues_loc.resize(dim_vector(Nentries, 1));
	}
	else {
		trues_loc.resize(dim_vector(1, 1));
	}
	if (store_randoms && source && outsize2 == 1) {
		randoms_loc.resize(dim_vector(Nentries, 1));
	}
	else {
		randoms_loc.resize(dim_vector(1, 1));
	}
	if (store_scatter && source && outsize2 == 1) {
		scatter_loc.resize(dim_vector(Nentries, 1));
	}
	else {
		scatter_loc.resize(dim_vector(1, 1));
	}
	octave_uint16* LL1_p = LL1.fortran_vec();
	octave_uint16* LL2_p = LL2.fortran_vec();
	octave_uint32* tpoints_p = tpoints.fortran_vec();
	bool* trues_loc_p, * randoms_loc_p, * scatter_loc_p;
	float* S_p = S.fortran_vec();
	octave_int32* int_loc_p = int_loc.fortran_vec();
	octave_uint16* Ltrues_p = Ltrues.fortran_vec();
	octave_uint16* Lrandoms_p = Lrandoms.fortran_vec();
	octave_uint16* Lscatter_p = Lscatter.fortran_vec();
	trues_loc_p = trues_loc.fortran_vec();
	randoms_loc_p = randoms_loc.fortran_vec();
	scatter_loc_p = scatter_loc.fortran_vec();
	octave_uint16* Ldelay1_p = Ldelay1.fortran_vec();
	octave_uint16* Ldelay2_p = Ldelay2.fortran_vec();
	octave_int32* int_loc_delay_p = int_loc_delay.fortran_vec();
	octave_uint32* tpoints_delay_p = tpoints_delay.fortran_vec();
	float* x1_p = x1.fortran_vec();
	float* x2_p = x2.fortran_vec();
	float* y1_p = y1.fortran_vec();
	float* y2_p = y2.fortran_vec();
	float* z1_p = z1.fortran_vec();
	float* z2_p = z2.fortran_vec();
	octave_uint16* Sino = SinoO.fortran_vec();
	octave_uint16* SinoT = SinoOT.fortran_vec();
	octave_uint16* SinoR = SinoOR.fortran_vec();
	octave_uint16* SinoC = SinoOC.fortran_vec();
	octave_uint16* SinoD = SinoOD.fortran_vec();

	histogram(LL1_p, LL2_p, tpoints_p, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S_p,
		Coincidences, Nentries, time_intervals_p, int_loc_p, obtain_trues, store_scatter, store_randoms, scatter_components_p, Ltrues_p, Lscatter_p,
		Lrandoms_p, trues_loc_p, Ndelays, randoms_correction, delay, Ldelay1_p, Ldelay2_p, int_loc_delay_p, tpoints_delay_p, randoms_loc_p, scatter_loc_p, 
		x1_p, x2_p, y1_p, y1_p, z1_p, z2_p, store_coordinates, dynamic, cryst_per_block_z, transaxial_multip, rings, sinoSize, Ndist, Nang, ringDifference, 
		span, seg_p, NT, TOFSize, nDistSide, storeRawData, Sino, SinoT, SinoC, SinoR, SinoD, detWPseudo, nPseudos, binSize);


	delete Coincidences;
	if (randoms_correction)
		delete delay;
	//mexEvalString("pause(.001);");
	//gROOT->Reset();


	octave_value_list retval(nargout);

	retval(0) = octave_value(LL1);
	retval(1) = octave_value(LL2);
	retval(2) = octave_value(tpoints);
	retval(3) = octave_value(S);
	retval(4) = octave_value(int_loc);
	retval(5) = octave_value(Ltrues);
	retval(6) = octave_value(Lscatter);
	retval(7) = octave_value(Lrandoms);
	retval(8) = octave_value(trues_loc);
	retval(9) = octave_value(Ldelay1);
	retval(10) = octave_value(Ldelay2);
	retval(11) = octave_value(int_loc_delay);
	retval(12) = octave_value(tpoints_delay);
	retval(13) = octave_value(randoms_loc);
	retval(14) = octave_value(scatter_loc);
	retval(15) = octave_value(x1);
	retval(16) = octave_value(x2);
	retval(17) = octave_value(y1);
	retval(18) = octave_value(y2);
	retval(19) = octave_value(z1);
	retval(20) = octave_value(z2);
	retval(21) = octave_value(SinoO);
	retval(22) = octave_value(SinoOT);
	retval(23) = octave_value(SinoOC);
	retval(24) = octave_value(SinoOR);
	retval(25) = octave_value(SinoOD);

	//gROOT->Reset();

	return retval;
}
