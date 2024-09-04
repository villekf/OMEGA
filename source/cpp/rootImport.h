#pragma once
#define NOMINMAX
#include <thread>
#include <TROOT.h>
#include "TTree.h"
#include "TFile.h"
#include "saveSinogram.h"
#include <charconv>
#ifdef MATLABCPP
#include "mex.hpp"
#include "mexAdapter.hpp"
void disp(const char* txt, const std::shared_ptr<matlab::engine::MATLABEngine>& matlabPtr) {
	matlab::data::ArrayFactory factory;
	std::ostringstream stream;
	stream << txt << std::endl;
	// Pass stream content to MATLAB fprintf function
	matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
	// Clear stream buffer
	//stream.str("");
}
#elif defined(MATLABC)
#include "mex.h"
template <typename T>
void disp(const char* txt, const T nullPar = NULL) {
	mexPrintf(txt);
	mexPrintf("\n");
}
template <typename C>
void dispf(const char* txt, const C var) {
	mexPrintf(txt, var);
	mexPrintf("\n");
}
#elif defined(OCTAVE)
#include <octave/oct.h>
template <typename T>
void disp(const char* txt, const T nullPar = NULL) {
	octave_stdout << txt;
	octave_stdout << "\n";
}
#else
template <typename T>
void disp(const char* txt, const T nullPar = NULL) {
	printf(txt);
	printf("\n");
}
#endif

template <typename T>
void formSourceImage(const float bx, const float by, const float bz, const float dx, const float dy, const float dz, const int64_t Nx, const int64_t Ny, const int64_t Nz, 
	const int64_t imDim, const float sourcePosX1, const float sourcePosX2, const float sourcePosY1, const float sourcePosY2, const float sourcePosZ1, const float sourcePosZ2, 
	const int64_t tPoint, T* S) {
	//float xa = bx;
	//float ya = by;
	//float za = bz;
	uint64_t indX = 0, indY = 0, indZ = 0;
	if (sourcePosX1 >= bx && sourcePosX1 <= bx + static_cast<float>(Nx) * dx)
		indX = static_cast<uint64_t>(std::floor((sourcePosX1 - bx) / dx));
	if (sourcePosY1 >= by && sourcePosY1 <= by + static_cast<float>(Ny) * dy)
		indY = static_cast<uint64_t>(std::floor((sourcePosY1 - by) / dy));
	if (sourcePosZ1 >= by && sourcePosZ1 <= bz + static_cast<float>(Nz) * dz)
		indZ = static_cast<uint64_t>(std::floor((sourcePosZ1 - bz) / dz));
	//for (uint64_t xi = 0; xi < Nx; xi++) {
	//	if (sourcePosX1 >= xa && sourcePosX1 < xa + dx) {
	//		indX = xi;
	//		break;
	//	}
	//	xa += dx;
	//}
	//for (uint64_t yi = 0; yi < Ny; yi++) {
	//	if (sourcePosY1 >= ya && sourcePosY1 < ya + dy) {
	//		indY = yi;
	//		break;
	//	}
	//	ya += dy;
	//}
	//for (uint64_t zi = 0; zi < Nz; zi++) {
	//	if (sourcePosZ1 >= za && sourcePosZ1 < za + dz) {
	//		indZ = zi;
	//		break;
	//	}
	//	za += dz;
	//}
//#pragma omp critical 
//	{
//		S[indX + indY * Nx + indZ * Nx * Ny + tPoint * imDim] = S[indX + indY * Nx + indZ * Nx * Ny + tPoint * imDim] + static_cast<T>(1);
//	}
#ifdef _OPENMP
#pragma omp atomic
#endif
	S[indX + indY * Nx + indZ * Nx * Ny + tPoint * imDim]++;
}

template <typename T, typename C, typename K, typename H, typename M, typename D>
void histogram(const char* rootFile, const C* tPoints, const double alku, const double loppu, bool source, const uint32_t linear_multip, const uint32_t* cryst_per_block, const uint32_t blocks_per_ring,
	const uint32_t* det_per_ring, T* S, T* SC, T* RA, T* trIndex, T* axIndex, T* DtrIndex, T* DaxIndex, bool obtain_trues, bool store_scatter, bool store_randoms, K* scatter_components,
	bool randoms_correction, M* coord, M* Dcoord, bool store_coordinates, const bool dynamic, const uint32_t* cryst_per_block_z,
	const uint32_t transaxial_multip, const uint32_t* rings, const uint64_t* sinoSize, const uint32_t Ndist, const uint32_t* Nang, const uint32_t ringDifference, const uint32_t span, 
	const H* seg, const int64_t Nt, const uint64_t TOFSize, const int32_t nDistSide, T* Sino, T* SinoT, T* SinoC, T* SinoR, T* SinoD, 
	const uint32_t* detWPseudo, const int32_t nPseudos, const double binSize, const double FWHM, const bool verbose, const int32_t nLayers, const float dx, const float dy, const float dz,
	const float bx, const float by, const float bz, const int64_t Nx, const int64_t Ny, const int64_t Nz, const bool dualLayerSubmodule, const int64_t imDim, const bool indexBased, T* tIndex, const D mPtr = 0) {

	int nthreads = 1;
	bool scatterTrues[] = {true, true, true, true};

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, FWHM + 1e-20);
	const uint64_t nBins = TOFSize / sinoSize[0];
	const bool TOF = nBins > 1;

#ifdef _OPENMP
	if (omp_get_max_threads() == 1) {
		int n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}
#endif

	Int_t moduleID1F = 0, moduleID2F = 0, submoduleID1F = 0, submoduleID2F = 0;

	TTree* Coincidences;
	TFile* inFile = new TFile(rootFile, "read");
	inFile->GetObject("Coincidences", Coincidences);


	int64_t Nentries;
	Nentries = Coincidences->GetEntries();

	if (Coincidences->GetBranchStatus("moduleID1"))
		Coincidences->SetBranchAddress("moduleID1", &moduleID1F);
	if (Coincidences->GetBranchStatus("moduleID2"))
		Coincidences->SetBranchAddress("moduleID2", &moduleID2F);
	if (Coincidences->GetBranchStatus("submoduleID1"))
		Coincidences->SetBranchAddress("submoduleID1", &submoduleID1F);
	if (Coincidences->GetBranchStatus("submoduleID2"))
		Coincidences->SetBranchAddress("submoduleID2", &submoduleID2F);
	uint64_t summa = 0ULL;
	uint64_t summaS = 0ULL;
	for (uint64_t kk = 0ULL; kk < std::min(static_cast<int64_t>(1000), Nentries); kk++) {
		Coincidences->GetEntry(kk);
		if (summa == 0ULL)
			summa += moduleID1F;
		if (summaS == 0ULL)
			summaS += submoduleID1F;
		if (summa > 0 && summaS > 0)
			break;
	}

	int any = 0;
	int next = 0;
	bool no_time = false;


	if (!Coincidences->GetBranchStatus("crystalID1")) {
		disp("No crystal location information was found from file. Aborting.", mPtr);
		return;
	}
	if (!Coincidences->GetBranchStatus("crystalID2")) {
		disp("No crystal location information was found from file. Aborting.", mPtr);
		return;
	}
	bool no_modules = false;
	bool no_submodules = true;
	if (summa == 0ULL)
		no_modules = true;
	if (summaS > 0ULL)
		no_submodules = false;
	bool layerSubmodule = false;
	if (dualLayerSubmodule && !no_submodules && nLayers > 1) {
		layerSubmodule = true;
		no_submodules = true;
	}
	const bool pseudoD = detWPseudo[0] > det_per_ring[0];
	const bool pseudoR = nPseudos > 0;
	int32_t gapSize = 0;
	if (pseudoR) {
		gapSize = rings[0] / (nPseudos + 1);
	}
	if (source) {
		if (!Coincidences->GetBranchStatus("sourcePosX1")) {
			disp("No X-source coordinates saved for first photon, unable to save source coordinates", mPtr);
			source = false;
		}
		if (!Coincidences->GetBranchStatus("sourcePosX2")) {
			disp("No X-source coordinates saved for second photon, unable to save source coordinates", mPtr);
			source = false;
		}
		if (!Coincidences->GetBranchStatus("sourcePosY1")) {
			disp("No Y-source coordinates saved for first photon, unable to save source coordinates", mPtr);
			source = false;
		}
		if (!Coincidences->GetBranchStatus("sourcePosY2")) {
			disp("No Y-source coordinates saved for second photon, unable to save source coordinates", mPtr);
			source = false;
		}
		if (!Coincidences->GetBranchStatus("sourcePosZ1")) {
			disp("No Z-source coordinates saved for first photon, unable to save source coordinates", mPtr);
			source = false;
		}
		if (!Coincidences->GetBranchStatus("sourcePosZ2")) {
			disp("No Z-source coordinates saved for second photon, unable to save source coordinates", mPtr);
			source = false;
		}
	}
	if (!Coincidences->GetBranchStatus("time1") && TOF) {
		disp("TOF examination selected, but no time information was found from file. Aborting.", mPtr);
		return;
	}
	if (!Coincidences->GetBranchStatus("time2") && (dynamic || TOF)) {
		disp("Dynamic or TOF examination selected, but no time information was found from file. Aborting.", mPtr);
	}
	if (!Coincidences->GetBranchStatus("time1") && !Coincidences->GetBranchStatus("time2"))
		no_time = true;
	if (store_coordinates) {
		if (!Coincidences->GetBranchStatus("globalPosX1")) {
			disp("No X-source coordinates saved for first photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
		if (!Coincidences->GetBranchStatus("globalPosX2")) {
			disp("No X-source coordinates saved for second photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
		if (!Coincidences->GetBranchStatus("globalPosY1")) {
			disp("No Y-source coordinates saved for first photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
		if (!Coincidences->GetBranchStatus("globalPosY2")) {
			disp("No Y-source coordinates saved for second photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
		if (!Coincidences->GetBranchStatus("globalPosZ1")) {
			disp("No Z-source coordinates saved for first photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
		if (!Coincidences->GetBranchStatus("globalPosZ2")) {
			disp("No Z-source coordinates saved for second photon interaction, unable to save interaction coordinates", mPtr);
			store_coordinates = false;
		}
	}
	if (obtain_trues || store_scatter || store_randoms) {
		if (!Coincidences->GetBranchStatus("eventID1")) {
			disp("No event IDs saved for first photon, unable to save trues/scatter/randoms", mPtr);
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
		if (!Coincidences->GetBranchStatus("eventID2")) {
			disp("No event IDs saved for second photon, unable to save trues/scatter/randoms", mPtr);
			obtain_trues = false;
			store_scatter = false;
			store_randoms = false;
		}
	}
	if (obtain_trues || store_scatter || store_randoms) {
		if (!Coincidences->GetBranchStatus("comptonPhantom1"))
			any++;
		if (!Coincidences->GetBranchStatus("comptonPhantom2"))
			any++;

		if (store_scatter && any == 2 && scatter_components[0] >= 1) {
			disp("Compton phantom selected, but no scatter data was found from ROOT-file", mPtr);
			scatter_components[0] = static_cast<K>(0);
		}
		else if (store_scatter && scatter_components[0] >= 1 && verbose) {
			disp("Compton scatter in the phantom will be stored", mPtr);
		}
		if (obtain_trues && any == 2) {
			scatterTrues[0] = false;
		}

		if (any == 2)
			next++;
		if (!Coincidences->GetBranchStatus("comptonCrystal1"))
			any++;
		if (!Coincidences->GetBranchStatus("comptonCrystal2"))
			any++;

		if (store_scatter && ((any == 4 && next == 1) || (any == 2 && next == 0)) && scatter_components[1] >= 1) {
			disp("Compton crystal selected, but no scatter data was found from ROOT-file", mPtr);
			scatter_components[1] = static_cast<K>(0);
		}
		else if (store_scatter && scatter_components[1] >= 1 && verbose) {
			disp("Compton scatter in the detector will be stored", mPtr);
		}
		if (obtain_trues && ((any == 4 && next == 1) || (any == 2 && next == 0))) {
			scatterTrues[1] = false;
		}

		if ((any == 4 && next == 1) || (any == 2 && next == 0))
			next++;
		if (!Coincidences->GetBranchStatus("RayleighPhantom1"))
			any++;
		if (!Coincidences->GetBranchStatus("RayleighPhantom2"))
			any++;

		if (store_scatter && ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1)) && scatter_components[2] >= 1) {
			disp("Rayleigh phantom selected, but no scatter data was found from ROOT-file", mPtr);
			scatter_components[2] = static_cast<K>(0);
		}
		else if (store_scatter && scatter_components[2] >= 1 && verbose) {
			disp("Rayleigh scatter in the phantom will be stored", mPtr);
		}
		if (obtain_trues && ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1))) {
			scatterTrues[2] = false;
		}

		if ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1))
			next++;
		if (!Coincidences->GetBranchStatus("RayleighCrystal1"))
			any++;
		if (!Coincidences->GetBranchStatus("RayleighCrystal2"))
			any++;

		if (store_scatter && ((any == 8 && next == 3) || (any == 2 && next == 0) || (any == 4 && next == 1) || (any == 6 && next == 2)) && scatter_components[3] >= 1) {
			disp("Rayleigh crystal selected, but no scatter data was found from ROOT-file", mPtr);
			scatter_components[3] = static_cast<K>(0);
		}
		else if (store_scatter && scatter_components[3] >= 1 && verbose) {
			disp("Rayleigh scatter in the detector will be stored", mPtr);
		}
		if (obtain_trues && ((any == 8 && next == 3) || (any == 2 && next == 0) || (any == 4 && next == 1) || (any == 6 && next == 2))) {
			scatterTrues[3] = false;
		}

		if (store_scatter && any == 8) {
			disp("Store scatter selected, but no scatter data was found from ROOT-file", mPtr);
			store_scatter = false;
		}

		if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 1 && scatterTrues[2] == 1 && scatterTrues[3] == 1 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the phantom and detector are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 1 && scatterTrues[2] == 1 && scatterTrues[3] == 0 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the phantom are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 1 && scatterTrues[2] == 0 && scatterTrues[3] == 0 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and detector are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 1 && scatterTrues[2] == 0 && scatterTrues[3] == 1 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the detector are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 0 && scatterTrues[2] == 1 && scatterTrues[3] == 1 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the phantom and detector are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 0 && scatterTrues[2] == 0 && scatterTrues[3] == 1 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the detector are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 0 && scatterTrues[2] == 1 && scatterTrues[3] == 0 && verbose) {
			disp("Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the phantom are NOT included in trues", mPtr);
		}
		else if (obtain_trues && scatterTrues[0] == 1 && scatterTrues[1] == 0 && scatterTrues[2] == 0 && scatterTrues[3] == 0 && verbose) {
			disp("Randoms and Compton scattered coincidences in the phantom are NOT included in trues", mPtr);
		}

	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static), num_threads(nthreads), shared(Coincidences)
#endif
	for (int64_t kk = 0; kk < Nentries; kk++) {

		Int_t crystalID1 = 0, crystalID2 = 0, moduleID1 = 0, moduleID2 = 0, submoduleID1 = 0, submoduleID2 = 0, rsectorID1, rsectorID2, eventID1, eventID2, comptonPhantom1 = 0, comptonPhantom2 = 0,
			comptonCrystal1 = 0, comptonCrystal2 = 0, RayleighPhantom1 = 0, RayleighPhantom2 = 0, RayleighCrystal1 = 0, RayleighCrystal2 = 0, layerID1 = 0, layerID2 = 0;
		Float_t sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
		Double_t time1 = alku, time2 = alku;
		int64_t tPoint = 0LL;
#ifdef _OPENMP
#pragma omp critical 
		{
#endif
			Coincidences->SetBranchAddress("rsectorID1", &rsectorID1);
			Coincidences->SetBranchAddress("rsectorID2", &rsectorID2);
			Coincidences->SetBranchAddress("crystalID1", &crystalID1);
			Coincidences->SetBranchAddress("crystalID2", &crystalID2);
			if (nLayers > 1) {
				Coincidences->SetBranchAddress("layerID1", &layerID1);
				Coincidences->SetBranchAddress("layerID2", &layerID2);
			}
			if (!no_modules) {
				Coincidences->SetBranchAddress("moduleID1", &moduleID1);
				Coincidences->SetBranchAddress("moduleID2", &moduleID2);
			}
			if (!no_submodules) {
				Coincidences->SetBranchAddress("submoduleID1", &submoduleID1);
				Coincidences->SetBranchAddress("submoduleID2", &submoduleID2);
			}
			else if (layerSubmodule) {
				Coincidences->SetBranchAddress("submoduleID1", &submoduleID1);
				Coincidences->SetBranchAddress("submoduleID2", &submoduleID2);
			}
			if (source) {
				Coincidences->SetBranchAddress("sourcePosX1", &sourcePosX1);
				Coincidences->SetBranchAddress("sourcePosX2", &sourcePosX2);
				Coincidences->SetBranchAddress("sourcePosY1", &sourcePosY1);
				Coincidences->SetBranchAddress("sourcePosY2", &sourcePosY2);
				Coincidences->SetBranchAddress("sourcePosZ1", &sourcePosZ1);
				Coincidences->SetBranchAddress("sourcePosZ2", &sourcePosZ2);
			}
			if (dynamic || TOF) {
				Coincidences->SetBranchAddress("time1", &time1);
				Coincidences->SetBranchAddress("time2", &time2);
			}
			if (store_coordinates) {
				Coincidences->SetBranchAddress("globalPosX1", &globalPosX1);
				Coincidences->SetBranchAddress("globalPosX2", &globalPosX2);
				Coincidences->SetBranchAddress("globalPosY1", &globalPosY1);
				Coincidences->SetBranchAddress("globalPosY2", &globalPosY2);
				Coincidences->SetBranchAddress("globalPosZ1", &globalPosZ1);
				Coincidences->SetBranchAddress("globalPosZ2", &globalPosZ2);
			}
			if (obtain_trues || store_scatter || store_randoms) {
				Coincidences->SetBranchAddress("eventID1", &eventID1);
				Coincidences->SetBranchAddress("eventID2", &eventID2);
				if (scatter_components[0] || scatterTrues[0])
					Coincidences->SetBranchAddress("comptonPhantom1", &comptonPhantom1);
				if (scatter_components[0] || scatterTrues[0])
					Coincidences->SetBranchAddress("comptonPhantom2", &comptonPhantom2);
				if (scatter_components[1] || scatterTrues[1])
					Coincidences->SetBranchAddress("comptonCrystal1", &comptonCrystal1);
				if (scatter_components[1] || scatterTrues[1])
					Coincidences->SetBranchAddress("comptonCrystal2", &comptonCrystal2);
				if (scatter_components[2] || scatterTrues[2])
					Coincidences->SetBranchAddress("RayleighPhantom1", &RayleighPhantom1);
				if (scatter_components[2] || scatterTrues[2])
					Coincidences->SetBranchAddress("RayleighPhantom2", &RayleighPhantom2);
				if (scatter_components[3] || scatterTrues[3])
					Coincidences->SetBranchAddress("RayleighCrystal1", &RayleighCrystal1);
				if (scatter_components[3] || scatterTrues[3])
					Coincidences->SetBranchAddress("RayleighCrystal2", &RayleighCrystal2);
			}
			Coincidences->GetEntry(kk);
#ifdef _OPENMP
		}
#endif
		///*
		if (!no_time && time2 < alku)
			continue;
		else if (!no_time && time2 > loppu) {
			continue;
		}
		if (nLayers > 1 && layerID1 > 0 && layerSubmodule)
			crystalID1 = submoduleID1;
		if (nLayers > 1 && layerID2 > 0 && layerSubmodule)
			crystalID2 = submoduleID2;
		uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
		detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
			submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block[layerID1], cryst_per_block[layerID2], cryst_per_block_z[layerID1], cryst_per_block_z[layerID2], transaxial_multip, rings[layerID1]);
		uint64_t bins = 0;
		bool event_true = true;
		bool event_scattered = true;
		bool store_scatter_event = false;
		if (obtain_trues || store_scatter || store_randoms) {
			if (eventID1 != eventID2) {
				event_true = false;
				event_scattered = false;
			}
			if (event_true) {
				if (comptonPhantom1 > 0 || comptonPhantom2 > 0) {
					event_true = false;
					if (scatter_components[0] > 0 && (scatter_components[0] <= comptonPhantom1 || scatter_components[0] <= comptonPhantom2))
						store_scatter_event = true;
				}
				else if ((comptonCrystal1 > 0 || comptonCrystal2 > 0)) {
					event_true = false;
					if (scatter_components[1] > 0 && (scatter_components[1] <= comptonCrystal1 || scatter_components[1] <= comptonCrystal2))
						store_scatter_event = true;
				}
				else if ((RayleighPhantom1 > 0 || RayleighPhantom2 > 0)) {
					event_true = false;
					if (scatter_components[2] > 0 && (scatter_components[2] <= RayleighPhantom1 || scatter_components[2] <= RayleighPhantom2))
						store_scatter_event = true;
				}
				else if ((RayleighCrystal1 > 0 || RayleighCrystal2 > 0)) {
					event_true = false;
					if (scatter_components[3] > 0 && (scatter_components[3] <= RayleighCrystal1 || scatter_components[3] <= RayleighCrystal2))
						store_scatter_event = true;
				}
				else
					event_scattered = false;
			}
		}
		if (dynamic) {
			double time = alku;
			for (int64_t ll = 0; ll < Nt; ll++) {
				if (time2 >= time || time2 < tPoints[0]) {
					tPoint = ll;
					break;
				}
				time += tPoints[ll];
			}
		}
		if (TOFSize > sinoSize[0]) {
			double timeDif = (time2 - time1);
			if (ring_pos2 > ring_pos1)
				timeDif = -timeDif;
			if (FWHM > 0.)
				timeDif += distribution(generator);
			if (std::abs(timeDif) > ((binSize / 2.) * static_cast<double>(nBins)))
				continue;
			bins = static_cast<uint64_t>(std::floor((std::abs(timeDif) + binSize / 2.) / binSize));
			const bool tInd = timeDif > 0;
			if (tInd)
				bins *= 2ULL;
			else if (!tInd && bins > 0ULL)
				bins = bins * 2ULL - 1ULL;
		}
		if (pseudoD) {
			ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
			ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
		}
		if (pseudoR) {
			ring_number1 += ring_number1 / gapSize;
			ring_number2 += ring_number2 / gapSize;
		}
		int32_t layer = 0;
		if (nLayers > 1) {
			if (layerID2 == 1 && layerID1 == 1)
				layer = 3;
			else if (layerID2 == 1 && layerID1 == 0)
				layer = 1;
			else if (layerID2 == 0 && layerID1 == 1)
				layer = 2;
			if (nLayers > 2) {
				if (layerID1 == 2 && layerID2 == 2)
					layer = 8;
				else if (layerID1 == 2 && layerID2 == 0)
					layer = 4;
				else if (layerID1 == 0 && layerID2 == 2)
					layer = 5;
				else if (layerID1 == 2 && layerID2 == 1)
					layer = 6;
				else if (layerID1 == 1 && layerID2 == 2)
					layer = 7;
			}
		}
		if (indexBased) {
			trIndex[kk * 2] = static_cast<uint16_t>(ring_pos1) + layerID1 * detWPseudo[0];
			trIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_pos2) + layerID2 * detWPseudo[0];
			axIndex[kk * 2] = static_cast<uint16_t>(ring_number1) + layerID1 * rings[0];
			axIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_number2) + layerID2 * rings[0];
		}
		else {
			if ((layer == 0 || layer == 1) && nLayers > 1) {
				ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
				ring_number1 += moduleID1;
			}
			if ((layer == 0 || layer == 2) && nLayers > 1) {
				ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
				ring_number2 += moduleID2;
			}
			bool swap = false;
			const int64_t sinoIndex = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize[0], Ndist, Nang[0], ringDifference, span, seg, TOFSize,
				detWPseudo[0], rings[0], bins, nDistSide, swap, tPoint, layer, nLayers);
			if (sinoIndex >= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
				Sino[sinoIndex]++;
				if ((event_true && obtain_trues) || (store_scatter_event && store_scatter)) {
					if (event_true && obtain_trues)
#ifdef _OPENMP
#pragma omp atomic
#endif
						SinoT[sinoIndex]++;
					else if (store_scatter_event && store_scatter)
#ifdef _OPENMP
#pragma omp atomic
#endif
						SinoC[sinoIndex]++;
				}
				else if (!event_true && store_randoms && !event_scattered)
#ifdef _OPENMP
#pragma omp atomic
#endif
					SinoR[sinoIndex]++;
			}
			if (source) {
				if (event_true && obtain_trues) {
					formSourceImage(bx, by, bz, dx, dy, dz, Nx, Ny, Nz, imDim, sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, tPoint, S);
				}
				else if (!obtain_trues) {
					if (sourcePosX1 == sourcePosX2 && sourcePosY1 == sourcePosY2 && sourcePosZ1 == sourcePosZ2) {
						formSourceImage(bx, by, bz, dx, dy, dz, Nx, Ny, Nz, imDim, sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, tPoint, S);
					}
				}
				if (store_scatter_event && store_scatter) {
					formSourceImage(bx, by, bz, dx, dy, dz, Nx, Ny, Nz, imDim, sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, tPoint, SC);
				}
				if (!event_true && !event_scattered && store_randoms) {
					formSourceImage(bx, by, bz, dx, dy, dz, Nx, Ny, Nz, imDim, sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2, tPoint, RA);
				}
			}
			if (store_coordinates) {
				coord[kk * 6] = globalPosX1;
				coord[kk * 6 + 1] = globalPosY1;
				coord[kk * 6 + 2] = globalPosZ1;
				coord[kk * 6 + 3] = globalPosX2;
				coord[kk * 6 + 4] = globalPosY2;
				coord[kk * 6 + 5] = globalPosZ2;
				if (dynamic)
					tIndex[kk] = static_cast<uint16_t>(tPoint);
			}
		}
	}


	if (randoms_correction) {

		TTree* delay;
		TFile* inFileD = new TFile(rootFile, "read");
		inFileD->GetObject("delay", delay);

		int64_t Ndelays = delay->GetEntries();


		for (int64_t kk = 0; kk < Ndelays; kk++) {
			Int_t crystalID1 = 0, crystalID2 = 0, moduleID1 = 0, moduleID2 = 0, submoduleID1 = 0, submoduleID2 = 0, rsectorID1, rsectorID2, layerID1 = 0, layerID2 = 0;
			Float_t globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
			Double_t time1 = alku, time2 = alku;
			int64_t tPoint = 0LL;
			delay->SetBranchAddress("crystalID1", &crystalID1);
			delay->SetBranchAddress("crystalID2", &crystalID2);
			if (!no_modules) {
				delay->SetBranchAddress("moduleID1", &moduleID1);
				delay->SetBranchAddress("moduleID2", &moduleID2);
			}
			if (!no_submodules) {
				delay->SetBranchAddress("submoduleID1", &submoduleID1);
				delay->SetBranchAddress("submoduleID2", &submoduleID2);
			}
			else if (layerSubmodule) {
				delay->SetBranchAddress("submoduleID1", &submoduleID1);
				delay->SetBranchAddress("submoduleID2", &submoduleID2);
			}
			delay->SetBranchAddress("rsectorID1", &rsectorID1);
			delay->SetBranchAddress("rsectorID2", &rsectorID2);
			if (dynamic) {
				if (delay->GetBranchStatus("time2"))
					delay->SetBranchAddress("time2", &time2);
			}
			if (nLayers > 1) {
				delay->SetBranchAddress("layerID1", &layerID1);
				delay->SetBranchAddress("layerID2", &layerID2);
			}
			if (store_coordinates) {
				delay->SetBranchAddress("globalPosX1", &globalPosX1);
				delay->SetBranchAddress("globalPosX2", &globalPosX2);
				delay->SetBranchAddress("globalPosY1", &globalPosY1);
				delay->SetBranchAddress("globalPosY2", &globalPosY2);
				delay->SetBranchAddress("globalPosZ1", &globalPosZ1);
				delay->SetBranchAddress("globalPosZ2", &globalPosZ2);
			}
			delay->GetEntry(kk);

			uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
			if (nLayers > 1 && layerID1 > 0 && layerSubmodule)
				crystalID1 = submoduleID1;
			if (nLayers > 1 && layerID2 > 0 && layerSubmodule)
				crystalID2 = submoduleID2;
			detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
				submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block[layerID1], cryst_per_block[layerID2], cryst_per_block_z[layerID1], cryst_per_block_z[layerID2], transaxial_multip, rings[layerID1]);
			uint64_t bins = 0;
			uint64_t L1 = static_cast<uint64_t>(ring_number1) * static_cast<uint64_t>(det_per_ring[layerID1]) + static_cast<uint64_t>(ring_pos1);
			uint64_t L2 = static_cast<uint64_t>(ring_number2) * static_cast<uint64_t>(det_per_ring[layerID1]) + static_cast<uint64_t>(ring_pos2);
			if (dynamic) {
				double time = alku;
				for (int64_t ll = 0; ll < Nt; ll++) {
					if (time2 >= time || time2 < tPoints[0]) {
						tPoint = ll;
						break;
					}
					time += tPoints[ll];
				}
			}
			if (pseudoD) {
				ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
				ring_pos2 += ring_pos2 / cryst_per_block[layerID1];
			}
			if (pseudoR) {
				ring_number1 += ring_number1 / gapSize;
				ring_number2 += ring_number2 / gapSize;
			}
			int32_t layer = 0;
			if (nLayers > 1) {
				if (layerID2 == 1 && layerID1 == 1)
					layer = 3;
				else if (layerID2 == 1 && layerID1 == 0)
					layer = 1;
				else if (layerID2 == 0 && layerID1 == 1)
					layer = 2;
				if (nLayers > 2) {
					if (layerID1 == 2 && layerID2 == 2)
						layer = 8;
					else if (layerID1 == 2 && layerID2 == 0)
						layer = 4;
					else if (layerID1 == 0 && layerID2 == 2)
						layer = 5;
					else if (layerID1 == 2 && layerID2 == 1)
						layer = 6;
					else if (layerID1 == 1 && layerID2 == 2)
						layer = 7;
				}
			}
			if (indexBased) {
				DtrIndex[kk * 2] = static_cast<uint16_t>(ring_pos1) + layerID1 * detWPseudo[0];
				DtrIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_pos2) + layerID2 * detWPseudo[0];
				DaxIndex[kk * 2] = static_cast<uint16_t>(ring_number1) + layerID1 * rings[0];
				DaxIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_number2) + layerID2 * rings[0];
			}
			else {
				if ((layer == 0 || layer == 1) && nLayers > 1) {
					ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
					ring_number1 += moduleID1;
				}
				if ((layer == 0 || layer == 2) && nLayers > 1) {
					ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
					ring_number2 += moduleID2;
			}
				bool swap = false;
				const int64_t sinoIndex = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize[0], Ndist, Nang[0], ringDifference, span, seg, TOFSize,
					detWPseudo[0], rings[0], bins, nDistSide, swap, tPoint, layer, nLayers);
				if (sinoIndex >= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
					SinoD[sinoIndex]++;
				}
				if (store_coordinates) {
					Dcoord[kk * 6] = globalPosX1;
					Dcoord[kk * 6 + 1] = globalPosY1;
					Dcoord[kk * 6 + 2] = globalPosZ1;
					Dcoord[kk * 6 + 3] = globalPosX2;
					Dcoord[kk * 6 + 4] = globalPosY2;
					Dcoord[kk * 6 + 5] = globalPosZ2;
				}
			}
		}
		delete inFileD;
	}
	delete inFile;
	return;
}