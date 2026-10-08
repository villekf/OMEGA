#pragma once
#define NOMINMAX
#include <thread>
#include <TROOT.h>
#include "TTree.h"
#include "TFile.h"
#include "saveSinogram.h"
#include <charconv>
#include <atomic>
#include <memory>
#include <vector>
#include <algorithm>
#include <random>
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
	printf("%s\n", txt);
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
	if (sourcePosZ1 >= bz && sourcePosZ1 <= bz + static_cast<float>(Nz) * dz)
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

#ifndef ROOT_IMPORT_MAX_THREADS
#define ROOT_IMPORT_MAX_THREADS 16
#endif
#ifndef ROOT_IMPORT_BLOCK_SIZE
#define ROOT_IMPORT_BLOCK_SIZE 2097152LL
#endif

enum RootIntCol : int { cRsector1, cRsector2, cCrystal1, cCrystal2, cModule1, cModule2, cSubmodule1, cSubmodule2, cLayer1, cLayer2, cEvent1, cEvent2,
	cComptonPhantom1, cComptonPhantom2, cComptonCrystal1, cComptonCrystal2, cRayleighPhantom1, cRayleighPhantom2, cRayleighCrystal1, cRayleighCrystal2, nRootIntCols };
enum RootFloatCol : int { cSourceX1, cSourceX2, cSourceY1, cSourceY2, cSourceZ1, cSourceZ2, cGlobalX1, cGlobalX2, cGlobalY1, cGlobalY2, cGlobalZ1, cGlobalZ2, nRootFloatCols };
enum RootDoubleCol : int { cTime1, cTime2, nRootDoubleCols };

static const char* const rootIntNames[nRootIntCols] = { "rsectorID1", "rsectorID2", "crystalID1", "crystalID2", "moduleID1", "moduleID2", "submoduleID1", "submoduleID2",
	"layerID1", "layerID2", "eventID1", "eventID2", "comptonPhantom1", "comptonPhantom2", "comptonCrystal1", "comptonCrystal2", "RayleighPhantom1", "RayleighPhantom2",
	"RayleighCrystal1", "RayleighCrystal2" };
static const char* const rootFloatNames[nRootFloatCols] = { "sourcePosX1", "sourcePosX2", "sourcePosY1", "sourcePosY2", "sourcePosZ1", "sourcePosZ2",
	"globalPosX1", "globalPosX2", "globalPosY1", "globalPosY2", "globalPosZ1", "globalPosZ2" };
static const char* const rootDoubleNames[nRootDoubleCols] = { "time1", "time2" };

// Reads a ROOT tree in blocks of entries into column vectors. Only the active branches are enabled (and thus decompressed).
// Each worker thread owns its own TFile/TTree, since TTree::GetEntry is not thread-safe on a shared tree.
class RootBlockReader {
public:
	bool intActive[nRootIntCols];
	bool floatActive[nRootFloatCols];
	bool doubleActive[nRootDoubleCols];
	std::vector<Int_t> ints[nRootIntCols];
	std::vector<Float_t> floats[nRootFloatCols];
	std::vector<Double_t> doubles[nRootDoubleCols];
	int64_t entries = 0;

	RootBlockReader() {
		std::fill(intActive, intActive + nRootIntCols, false);
		std::fill(floatActive, floatActive + nRootFloatCols, false);
		std::fill(doubleActive, doubleActive + nRootDoubleCols, false);
	}

	~RootBlockReader() {
		for (int w = 0; w < nWorkers; w++) {
			delete workers[w].file;
			workers[w].file = nullptr;
		}
	}

	bool open(const char* fileName, const char* treeName, int nThreads, int64_t blockSize) {
		nWorkers = std::max(1, nThreads);
		// Allocated once, branch addresses point into the workers
		workers.reset(new Worker[nWorkers]);
		for (int w = 0; w < nWorkers; w++) {
			Worker& wk = workers[w];
			wk.file = TFile::Open(fileName, "READ");
			if (wk.file == nullptr || wk.file->IsZombie())
				return false;
			wk.file->GetObject(treeName, wk.tree);
			if (wk.tree == nullptr)
				return false;
			wk.tree->SetBranchStatus("*", 0);
			for (int c = 0; c < nRootIntCols; c++) {
				if (!intActive[c])
					continue;
				if (wk.tree->GetBranch(rootIntNames[c]) == nullptr) {
					intActive[c] = false;
					continue;
				}
				wk.tree->SetBranchStatus(rootIntNames[c], 1);
				if (wk.tree->SetBranchAddress(rootIntNames[c], &wk.iBuf[c]) < 0)
					return false;
			}
			for (int c = 0; c < nRootFloatCols; c++) {
				if (!floatActive[c])
					continue;
				if (wk.tree->GetBranch(rootFloatNames[c]) == nullptr) {
					floatActive[c] = false;
					continue;
				}
				wk.tree->SetBranchStatus(rootFloatNames[c], 1);
				if (wk.tree->SetBranchAddress(rootFloatNames[c], &wk.fBuf[c]) < 0)
					return false;
			}
			for (int c = 0; c < nRootDoubleCols; c++) {
				if (!doubleActive[c])
					continue;
				if (wk.tree->GetBranch(rootDoubleNames[c]) == nullptr) {
					doubleActive[c] = false;
					continue;
				}
				wk.tree->SetBranchStatus(rootDoubleNames[c], 1);
				if (wk.tree->SetBranchAddress(rootDoubleNames[c], &wk.dBuf[c]) < 0)
					return false;
			}
		}
		entries = workers[0].tree->GetEntries();
		for (int c = 0; c < nRootIntCols; c++)
			if (intActive[c])
				ints[c].resize(blockSize);
		for (int c = 0; c < nRootFloatCols; c++)
			if (floatActive[c])
				floats[c].resize(blockSize);
		for (int c = 0; c < nRootDoubleCols; c++)
			if (doubleActive[c])
				doubles[c].resize(blockSize);
		return true;
	}

	// Reads entries [start, start + n) into the column vectors (index 0 corresponds to entry start)
	bool readBlock(const int64_t start, const int64_t n) {
		std::atomic<bool> fail(false);
		const int64_t chunk = (n + nWorkers - 1) / nWorkers;
		std::vector<std::thread> threads;
		for (int w = 1; w < nWorkers; w++) {
			const int64_t b = std::min<int64_t>(n, chunk * w);
			const int64_t e = std::min<int64_t>(n, chunk * (w + 1));
			if (e > b)
				threads.emplace_back(&RootBlockReader::readChunk, this, w, start, b, e, &fail);
		}
		readChunk(0, start, 0, std::min<int64_t>(n, chunk), &fail);
		for (auto& t : threads)
			t.join();
		return !fail.load();
	}

private:
	struct Worker {
		TFile* file = nullptr;
		TTree* tree = nullptr;
		Int_t iBuf[nRootIntCols] = {};
		Float_t fBuf[nRootFloatCols] = {};
		Double_t dBuf[nRootDoubleCols] = {};
	};
	std::unique_ptr<Worker[]> workers;
	int nWorkers = 0;

	void readChunk(const int w, const int64_t start, const int64_t b, const int64_t e, std::atomic<bool>* fail) {
		Worker& wk = workers[w];
		for (int64_t i = b; i < e; i++) {
			if (wk.tree->GetEntry(start + i) <= 0) {
				fail->store(true);
				return;
			}
			for (int c = 0; c < nRootIntCols; c++)
				if (intActive[c])
					ints[c][i] = wk.iBuf[c];
			for (int c = 0; c < nRootFloatCols; c++)
				if (floatActive[c])
					floats[c][i] = wk.fBuf[c];
			for (int c = 0; c < nRootDoubleCols; c++)
				if (doubleActive[c])
					doubles[c][i] = wk.dBuf[c];
		}
	}
};

template <typename T, typename C, typename K, typename H, typename M, typename D>
void histogram(const char* rootFile, const C* tPoints, const double alku, const double loppu, bool source, const uint32_t linear_multip, const uint32_t* cryst_per_block, const uint32_t blocks_per_ring,
	const uint32_t* det_per_ring, T* S, T* SC, T* RA, T* trIndex, T* axIndex, T* DtrIndex, T* DaxIndex, bool obtain_trues, bool store_scatter, bool store_randoms, K* scatter_components,
	bool randoms_correction, M* coord, M* Dcoord, bool store_coordinates, const bool dynamic, const uint32_t* cryst_per_block_z,
	const uint32_t transaxial_multip, const uint32_t* rings, const uint64_t* sinoSize, const uint32_t Ndist, const uint32_t* Nang, const uint32_t ringDifference, const uint32_t span, 
	const H* seg, const int64_t Nt, const uint64_t TOFSize, const int32_t nDistSide, T* Sino, T* SinoT, T* SinoC, T* SinoR, T* SinoD, 
	const uint32_t* detWPseudo, const int32_t nPseudos, const double binSize, const double FWHM, const bool verbose, const int32_t nLayers, const float dx, const float dy, const float dz,
	const float bx, const float by, const float bz, const int64_t Nx, const int64_t Ny, const int64_t Nz, const bool dualLayerSubmodule, const int64_t imDim, const bool indexBased, T* tIndex, 
	uint8_t* TOFIndex, const D mPtr) {

	int nThreads = std::max(1, std::min<int>(static_cast<int>(std::thread::hardware_concurrency()), ROOT_IMPORT_MAX_THREADS));
	if (nThreads > 1)
		ROOT::EnableThreadSafety();
	bool scatterTrues[] = {true, true, true, true};

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, FWHM + 1e-20);
	const uint64_t nBins = TOFSize / sinoSize[0];
	const bool TOF = nBins > 1;
	// Custom time window in a static (non-dynamic) examination
	const bool customWindow = alku > 0. || loppu < 1e9;

#ifdef _OPENMP
	if (omp_get_max_threads() == 1) {
		int n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}
#endif

	Int_t moduleID1F = 0, submoduleID1F = 0;

	TTree* Coincidences;
	TFile* inFile = new TFile(rootFile, "read");
	inFile->GetObject("Coincidences", Coincidences);


	int64_t Nentries;
	Nentries = Coincidences->GetEntries();

	TBranch* bMod = nullptr;
	TBranch* bSub = nullptr;
	Coincidences->SetBranchAddress("moduleID1", &moduleID1F, &bMod);
	Coincidences->SetBranchAddress("submoduleID1", &submoduleID1F, &bSub);
	uint64_t summa = 0ULL;
	uint64_t summaS = 0ULL;
	for (uint64_t kk = 0ULL; kk < std::min(static_cast<int64_t>(1000), Nentries); kk++) {
		if (bMod != nullptr)
			bMod->GetEntry(kk);
		if (bSub != nullptr)
			bSub->GetEntry(kk);
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
		delete inFile;
		return;
	}
	if (!Coincidences->GetBranchStatus("crystalID2")) {
		disp("No crystal location information was found from file. Aborting.", mPtr);
		delete inFile;
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
		delete inFile;
		return;
	}
	if (!Coincidences->GetBranchStatus("time2") && (dynamic || TOF)) {
		disp("Dynamic or TOF examination selected, but no time information was found from file. Aborting.", mPtr);
		delete inFile;
		return;
	}
	if (!Coincidences->GetBranchStatus("time1") && !Coincidences->GetBranchStatus("time2"))
		no_time = true;
	const bool timeWindow = customWindow && !dynamic && Coincidences->GetBranchStatus("time2");
	if (customWindow && !dynamic) {
		char windowMsg[256];
		if (timeWindow)
			snprintf(windowMsg, sizeof(windowMsg), "Using time window from %g s to %g s", alku, loppu);
		else
			snprintf(windowMsg, sizeof(windowMsg), "A custom time window was selected, but no time information was found from file. The time window cannot be applied.");
		disp(windowMsg, mPtr);
	}
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

	nThreads = std::max<int64_t>(1, std::min<int64_t>(nThreads, Nentries / 100000));
	{
		RootBlockReader reader;
		reader.intActive[cRsector1] = true;
		reader.intActive[cRsector2] = true;
		reader.intActive[cCrystal1] = true;
		reader.intActive[cCrystal2] = true;
		if (nLayers > 1) {
			reader.intActive[cLayer1] = true;
			reader.intActive[cLayer2] = true;
		}
		if (!no_modules) {
			reader.intActive[cModule1] = true;
			reader.intActive[cModule2] = true;
		}
		if (!no_submodules || layerSubmodule) {
			reader.intActive[cSubmodule1] = true;
			reader.intActive[cSubmodule2] = true;
		}
		if (source) {
			for (int c = cSourceX1; c <= cSourceZ2; c++)
				reader.floatActive[c] = true;
		}
		if (dynamic || TOF)
			reader.doubleActive[cTime1] = true;
		if (dynamic || TOF || timeWindow)
			reader.doubleActive[cTime2] = true;
		if (store_coordinates) {
			for (int c = cGlobalX1; c <= cGlobalZ2; c++)
				reader.floatActive[c] = true;
		}
		if (obtain_trues || store_scatter || store_randoms) {
			reader.intActive[cEvent1] = true;
			reader.intActive[cEvent2] = true;
			if (scatter_components[0] || scatterTrues[0]) {
				reader.intActive[cComptonPhantom1] = true;
				reader.intActive[cComptonPhantom2] = true;
			}
			if (scatter_components[1] || scatterTrues[1]) {
				reader.intActive[cComptonCrystal1] = true;
				reader.intActive[cComptonCrystal2] = true;
			}
			if (scatter_components[2] || scatterTrues[2]) {
				reader.intActive[cRayleighPhantom1] = true;
				reader.intActive[cRayleighPhantom2] = true;
			}
			if (scatter_components[3] || scatterTrues[3]) {
				reader.intActive[cRayleighCrystal1] = true;
				reader.intActive[cRayleighCrystal2] = true;
			}
		}
		const int64_t blockSize = std::max<int64_t>(1, std::min<int64_t>(ROOT_IMPORT_BLOCK_SIZE, Nentries));
		if (Nentries > 0 && !reader.open(rootFile, "Coincidences", nThreads, blockSize)) {
			disp("Error opening the ROOT file for reading", mPtr);
			delete inFile;
			return;
		}
		for (int64_t blockStart = 0; blockStart < Nentries; blockStart += blockSize) {
			const int64_t nBlock = std::min<int64_t>(blockSize, Nentries - blockStart);
			if (!reader.readBlock(blockStart, nBlock)) {
				disp("Error reading the ROOT file", mPtr);
				break;
			}
			for (int64_t ll = 0; ll < nBlock; ll++) {
				const int64_t kk = blockStart + ll;
				const Int_t rsectorID1 = reader.intActive[cRsector1] ? reader.ints[cRsector1][ll] : 0;
				const Int_t rsectorID2 = reader.intActive[cRsector2] ? reader.ints[cRsector2][ll] : 0;
				Int_t crystalID1 = reader.intActive[cCrystal1] ? reader.ints[cCrystal1][ll] : 0;
				Int_t crystalID2 = reader.intActive[cCrystal2] ? reader.ints[cCrystal2][ll] : 0;
				const Int_t moduleID1 = reader.intActive[cModule1] ? reader.ints[cModule1][ll] : 0;
				const Int_t moduleID2 = reader.intActive[cModule2] ? reader.ints[cModule2][ll] : 0;
				const Int_t submoduleID1 = reader.intActive[cSubmodule1] ? reader.ints[cSubmodule1][ll] : 0;
				const Int_t submoduleID2 = reader.intActive[cSubmodule2] ? reader.ints[cSubmodule2][ll] : 0;
				const Int_t layerID1 = reader.intActive[cLayer1] ? reader.ints[cLayer1][ll] : 0;
				const Int_t layerID2 = reader.intActive[cLayer2] ? reader.ints[cLayer2][ll] : 0;
				const Int_t eventID1 = reader.intActive[cEvent1] ? reader.ints[cEvent1][ll] : 0;
				const Int_t eventID2 = reader.intActive[cEvent2] ? reader.ints[cEvent2][ll] : 0;
				const Int_t comptonPhantom1 = reader.intActive[cComptonPhantom1] ? reader.ints[cComptonPhantom1][ll] : 0;
				const Int_t comptonPhantom2 = reader.intActive[cComptonPhantom2] ? reader.ints[cComptonPhantom2][ll] : 0;
				const Int_t comptonCrystal1 = reader.intActive[cComptonCrystal1] ? reader.ints[cComptonCrystal1][ll] : 0;
				const Int_t comptonCrystal2 = reader.intActive[cComptonCrystal2] ? reader.ints[cComptonCrystal2][ll] : 0;
				const Int_t RayleighPhantom1 = reader.intActive[cRayleighPhantom1] ? reader.ints[cRayleighPhantom1][ll] : 0;
				const Int_t RayleighPhantom2 = reader.intActive[cRayleighPhantom2] ? reader.ints[cRayleighPhantom2][ll] : 0;
				const Int_t RayleighCrystal1 = reader.intActive[cRayleighCrystal1] ? reader.ints[cRayleighCrystal1][ll] : 0;
				const Int_t RayleighCrystal2 = reader.intActive[cRayleighCrystal2] ? reader.ints[cRayleighCrystal2][ll] : 0;
				const Float_t sourcePosX1 = reader.floatActive[cSourceX1] ? reader.floats[cSourceX1][ll] : 0.f;
				const Float_t sourcePosX2 = reader.floatActive[cSourceX2] ? reader.floats[cSourceX2][ll] : 0.f;
				const Float_t sourcePosY1 = reader.floatActive[cSourceY1] ? reader.floats[cSourceY1][ll] : 0.f;
				const Float_t sourcePosY2 = reader.floatActive[cSourceY2] ? reader.floats[cSourceY2][ll] : 0.f;
				const Float_t sourcePosZ1 = reader.floatActive[cSourceZ1] ? reader.floats[cSourceZ1][ll] : 0.f;
				const Float_t sourcePosZ2 = reader.floatActive[cSourceZ2] ? reader.floats[cSourceZ2][ll] : 0.f;
				const Float_t globalPosX1 = reader.floatActive[cGlobalX1] ? reader.floats[cGlobalX1][ll] : 0.f;
				const Float_t globalPosX2 = reader.floatActive[cGlobalX2] ? reader.floats[cGlobalX2][ll] : 0.f;
				const Float_t globalPosY1 = reader.floatActive[cGlobalY1] ? reader.floats[cGlobalY1][ll] : 0.f;
				const Float_t globalPosY2 = reader.floatActive[cGlobalY2] ? reader.floats[cGlobalY2][ll] : 0.f;
				const Float_t globalPosZ1 = reader.floatActive[cGlobalZ1] ? reader.floats[cGlobalZ1][ll] : 0.f;
				const Float_t globalPosZ2 = reader.floatActive[cGlobalZ2] ? reader.floats[cGlobalZ2][ll] : 0.f;
				const Double_t time1 = reader.doubleActive[cTime1] ? reader.doubles[cTime1][ll] : alku;
				const Double_t time2 = reader.doubleActive[cTime2] ? reader.doubles[cTime2][ll] : alku;
				int64_t tPoint = 0LL;
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
					tPoint = Nt - 1;
					for (int64_t ll = 0; ll < Nt; ll++) {
						time += tPoints[ll];
						if (time2 < time) {
							tPoint = ll;
							break;
						}
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
					// Index-based TOF indexing also swaps the TOF directions when needed
					// The behavior should be the same to the sinogram version
					if (ring_pos2 < ring_pos1) {
						trIndex[kk * 2] = static_cast<uint16_t>(ring_pos2) + layerID2 * detWPseudo[0];
						trIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_pos1) + layerID1 * detWPseudo[0];
						axIndex[kk * 2] = static_cast<uint16_t>(ring_number2) + layerID2 * rings[0];
						axIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_number1) + layerID1 * rings[0];
					}
					else {
						trIndex[kk * 2] = static_cast<uint16_t>(ring_pos1) + layerID1 * detWPseudo[0];
						trIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_pos2) + layerID2 * detWPseudo[0];
						axIndex[kk * 2] = static_cast<uint16_t>(ring_number1) + layerID1 * rings[0];
						axIndex[kk * 2 + 1] = static_cast<uint16_t>(ring_number2) + layerID2 * rings[0];
					}
					if (dynamic)
						tIndex[kk] = static_cast<uint16_t>(tPoint);
					if (TOFSize > sinoSize[0])
						TOFIndex[kk] = static_cast<uint8_t>(bins);
				}
				else {
					if (pseudoD) {
						ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
						ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
					}
					if (pseudoR) {
						ring_number1 += ring_number1 / gapSize;
						ring_number2 += ring_number2 / gapSize;
					}
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
						if (TOFSize > sinoSize[0])
							TOFIndex[kk] = static_cast<uint8_t>(bins);
					}
				}
			}
		}
	}


	if (randoms_correction) {

		RootBlockReader reader;
		reader.intActive[cCrystal1] = true;
		reader.intActive[cCrystal2] = true;
		reader.intActive[cRsector1] = true;
		reader.intActive[cRsector2] = true;
		if (!no_modules) {
			reader.intActive[cModule1] = true;
			reader.intActive[cModule2] = true;
		}
		if (!no_submodules || layerSubmodule) {
			reader.intActive[cSubmodule1] = true;
			reader.intActive[cSubmodule2] = true;
		}
		TTree* delay = nullptr;
		inFile->GetObject("delay", delay);
		const bool delayWindow = (dynamic || customWindow) && delay != nullptr && delay->GetBranch("time2") != nullptr;
		if (delayWindow)
			reader.doubleActive[cTime2] = true;
		if (nLayers > 1) {
			reader.intActive[cLayer1] = true;
			reader.intActive[cLayer2] = true;
		}
		if (store_coordinates) {
			for (int c = cGlobalX1; c <= cGlobalZ2; c++)
				reader.floatActive[c] = true;
		}
		const int64_t Ndelays = (delay != nullptr) ? delay->GetEntries() : 0;
		int nThreadsD = std::max(1, std::min<int>(static_cast<int>(std::thread::hardware_concurrency()), ROOT_IMPORT_MAX_THREADS));
		nThreadsD = std::max<int64_t>(1, std::min<int64_t>(nThreadsD, Ndelays / 100000));
		const int64_t blockSizeD = std::max<int64_t>(1, std::min<int64_t>(ROOT_IMPORT_BLOCK_SIZE, Ndelays));
		if (!reader.open(rootFile, "delay", nThreadsD, blockSizeD)) {
			disp("Error opening the delayed coincidences from the ROOT file", mPtr);
			delete inFile;
			return;
		}
		for (int64_t blockStart = 0; blockStart < Ndelays; blockStart += blockSizeD) {
			const int64_t nBlock = std::min<int64_t>(blockSizeD, Ndelays - blockStart);
			if (!reader.readBlock(blockStart, nBlock)) {
				disp("Error reading the delayed coincidences from the ROOT file", mPtr);
				break;
			}
			for (int64_t ll = 0; ll < nBlock; ll++) {
				const int64_t kk = blockStart + ll;
				Int_t crystalID1 = reader.intActive[cCrystal1] ? reader.ints[cCrystal1][ll] : 0;
				Int_t crystalID2 = reader.intActive[cCrystal2] ? reader.ints[cCrystal2][ll] : 0;
				const Int_t moduleID1 = reader.intActive[cModule1] ? reader.ints[cModule1][ll] : 0;
				const Int_t moduleID2 = reader.intActive[cModule2] ? reader.ints[cModule2][ll] : 0;
				const Int_t submoduleID1 = reader.intActive[cSubmodule1] ? reader.ints[cSubmodule1][ll] : 0;
				const Int_t submoduleID2 = reader.intActive[cSubmodule2] ? reader.ints[cSubmodule2][ll] : 0;
				const Int_t rsectorID1 = reader.intActive[cRsector1] ? reader.ints[cRsector1][ll] : 0;
				const Int_t rsectorID2 = reader.intActive[cRsector2] ? reader.ints[cRsector2][ll] : 0;
				const Int_t layerID1 = reader.intActive[cLayer1] ? reader.ints[cLayer1][ll] : 0;
				const Int_t layerID2 = reader.intActive[cLayer2] ? reader.ints[cLayer2][ll] : 0;
				const Float_t globalPosX1 = reader.floatActive[cGlobalX1] ? reader.floats[cGlobalX1][ll] : 0.f;
				const Float_t globalPosX2 = reader.floatActive[cGlobalX2] ? reader.floats[cGlobalX2][ll] : 0.f;
				const Float_t globalPosY1 = reader.floatActive[cGlobalY1] ? reader.floats[cGlobalY1][ll] : 0.f;
				const Float_t globalPosY2 = reader.floatActive[cGlobalY2] ? reader.floats[cGlobalY2][ll] : 0.f;
				const Float_t globalPosZ1 = reader.floatActive[cGlobalZ1] ? reader.floats[cGlobalZ1][ll] : 0.f;
				const Float_t globalPosZ2 = reader.floatActive[cGlobalZ2] ? reader.floats[cGlobalZ2][ll] : 0.f;
				const Double_t time1 = alku;
				const Double_t time2 = reader.doubleActive[cTime2] ? reader.doubles[cTime2][ll] : alku;
				int64_t tPoint = 0LL;
				if (delayWindow && (time2 < alku || time2 > loppu))
					continue;
				uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
				if (nLayers > 1 && layerID1 > 0 && layerSubmodule)
					crystalID1 = submoduleID1;
				if (nLayers > 1 && layerID2 > 0 && layerSubmodule)
					crystalID2 = submoduleID2;
				detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
					submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block[layerID1], cryst_per_block[layerID2], cryst_per_block_z[layerID1], cryst_per_block_z[layerID2], transaxial_multip, rings[layerID1]);
				uint64_t bins = 0;
				if (dynamic) {
					double time = alku;
					tPoint = Nt - 1;
					for (int64_t ll = 0; ll < Nt; ll++) {
						time += tPoints[ll];
						if (time2 < time) {
							tPoint = ll;
							break;
						}
					}
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
					if (pseudoD) {
						ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
						ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
					}
					if (pseudoR) {
						ring_number1 += ring_number1 / gapSize;
						ring_number2 += ring_number2 / gapSize;
					}
					if ((layer == 0 || layer == 1) && nLayers > 1) {
						ring_pos1 += ring_pos1 / cryst_per_block[layerID1];
						ring_number1 += moduleID1;
					}
					if ((layer == 0 || layer == 2) && nLayers > 1) {
						ring_pos2 += ring_pos2 / cryst_per_block[layerID2];
						ring_number2 += moduleID2;
				}
					bool swap = false;
					const int64_t sinoIndex = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize[0], Ndist, Nang[0], ringDifference, span, seg, sinoSize[0],
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
		}
	}
	delete inFile;
	return;
}