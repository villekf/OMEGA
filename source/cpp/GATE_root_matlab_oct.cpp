// Load ROOT data in Octave
#include <octave/oct.h>
#include "TChain.h"
#include "rootImport.h"

DEFUN_DLD(GATE_root_matlab_oct, prhs, nargout, "GATE ROOT help") {


	NDArray tPoints = prhs(1).array_value();
	const double alku = prhs(2).scalar_value();
	const double loppu = prhs(3).scalar_value();
	uint32NDArray detectorsP = prhs(4).uint32_array_value();
	const uint32_t blocks_per_ring = prhs(5).uint32_scalar_value();
	uint32NDArray cryst_per_blockP = prhs(6).uint32_array_value();
	uint32NDArray det_per_ringP = prhs(7).uint32_array_value();
	const uint32_t linear_multp = prhs(8).uint32_scalar_value();
	bool source = prhs(9).bool_value();
	const int64_t Nt = prhs(10).int64_scalar_value();
	bool obtain_trues = prhs(11).bool_value();
	bool store_scatter = prhs(12).bool_value();
	bool store_randoms = prhs(13).bool_value();
	uint8NDArray scatter_components = prhs(14).uint8_array_value();
	bool randoms_correction = prhs(15).bool_value();
	bool store_coordinates = prhs(16).bool_value();
	const uint32_t transaxial_multip = prhs(17).uint32_scalar_value();
	uint32NDArray cryst_per_block_zP = prhs(18).uint32_array_value();
	uint32NDArray ringsP = prhs(19).uint32_array_value();
	uint64NDArray sinoSizeP = prhs(20).uint64_array_value();
	const bool TOF = prhs(21).bool_value();
	const bool verbose = prhs(22).bool_value();
	const int32_t nLayers = prhs(23).int32_scalar_value();
	const uint32_t Ndist = prhs(24).uint32_scalar_value();
	uint32NDArray NangP = prhs(25).uint32_array_value();
	const uint32_t ringDifference = prhs(26).uint32_scalar_value();
	const uint32_t span = prhs(27).uint32_scalar_value();
	uint32NDArray seg = prhs(28).uint32_array_value();
	const uint64_t TOFSize = prhs(29).uint64_scalar_value();
	const int32_t nDistSide = prhs(30).int32_scalar_value();
	uint16NDArray SinoO = prhs(31).uint16_array_value();
	uint16NDArray SinoOT = prhs(32).uint16_array_value();
	uint16NDArray SinoOC = prhs(33).uint16_array_value();
	uint16NDArray SinoOR = prhs(34).uint16_array_value();
	uint16NDArray SinoOD = prhs(35).uint16_array_value();
	uint32NDArray detWPseudoP = prhs(36).uint32_array_value();
	const int32_t nPseudos = prhs(37).int32_scalar_value();
	const double binSize = prhs(38).scalar_value();
	const double FWHM = prhs(39).scalar_value();
	const float dx = prhs(40).float_value();
	const float dy = prhs(41).float_value();
	const float dz = prhs(42).float_value();
	const float bx = prhs(43).float_value();
	const float by = prhs(44).float_value();
	const float bz = prhs(45).float_value();
	const int64_t Nx = prhs(46).int64_scalar_value();
	const int64_t Ny = prhs(47).int64_scalar_value();
	const int64_t Nz = prhs(48).int64_scalar_value();
	const bool dualLayerSubmodule = prhs(49).bool_value();
	const bool indexBased = prhs(50).bool_value();

	uint8_t* scatter_components_p = reinterpret_cast<uint8_t*>(scatter_components.fortran_vec());
	const double* tPoints_p = tPoints.fortran_vec();
	const uint32_t* seg_p = reinterpret_cast<uint32_t*>(seg.fortran_vec());
	const uint32_t* detectors = reinterpret_cast<uint32_t*>(detectorsP.fortran_vec());
	const uint32_t* cryst_per_block = reinterpret_cast<uint32_t*>(cryst_per_blockP.fortran_vec());
	const uint32_t* det_per_ring = reinterpret_cast<uint32_t*>(det_per_ringP.fortran_vec());
	const uint32_t* cryst_per_block_z = reinterpret_cast<uint32_t*>(cryst_per_block_zP.fortran_vec());
	const uint32_t* rings = reinterpret_cast<uint32_t*>(ringsP.fortran_vec());
	const uint64_t* sinoSize = reinterpret_cast<uint64_t*>(sinoSizeP.fortran_vec());
	const uint32_t* Nang = reinterpret_cast<uint32_t*>(NangP.fortran_vec());
	const uint32_t* detWPseudo = reinterpret_cast<uint32_t*>(detWPseudoP.fortran_vec());

	const int64_t imDim = Nx * Ny * Nz;

	const bool dynamic = Nt > 1;

	/* Pointer to character array */
	charMatrix apu = prhs(0).char_matrix_value();
	std::string tmp = apu.row_as_string(0);

	TChain *Coincidences = new TChain("Coincidences");
	Coincidences->Add(tmp.c_str());

	int64_t Nentries = Coincidences->GetEntries();
	TChain* delay = nullptr;
	int64_t Ndelays = 0LL;

	if (randoms_correction) {
		delay = new TChain("delay");
		delay->Add(tmp.c_str());
		Ndelays = delay->GetEntries();
	}

	if (randoms_correction)
		delete delay;
	delete Coincidences;

	/* Assign pointers to the various parameters */
	uint16NDArray trIndex;
	uint16NDArray axIndex;
	uint16NDArray DtrIndex;
	uint16NDArray DaxIndex;
	FloatNDArray coord;
	FloatNDArray Dcoord;
	if (store_coordinates) {
		coord.resize(dim_vector(6, Nentries));
		if (randoms_correction)
			Dcoord.resize(dim_vector(6, Nentries));
		else
			Dcoord.resize(dim_vector(1, 1));
	}
	else {
		coord.resize(dim_vector(1, 1));
		Dcoord.resize(dim_vector(1, 1));
	}
	if (indexBased) {
		trIndex.resize(dim_vector(2, Nentries));
		axIndex.resize(dim_vector(2, Nentries));
		if (randoms_correction) {
			DtrIndex.resize(dim_vector(2, Nentries));
			DaxIndex.resize(dim_vector(2, Nentries));
		}
		else {
			DtrIndex.resize(dim_vector(1, 1));
			DaxIndex.resize(dim_vector(1, 1));
		}
	}
	else {
		trIndex.resize(dim_vector(1, 1));
		axIndex.resize(dim_vector(1, 1));
		DtrIndex.resize(dim_vector(1, 1));
		DaxIndex.resize(dim_vector(1, 1));
	}
	uint16NDArray tIndex;
	uint16NDArray S, RA, SC;

	if (source) {
		S.resize(dim_vector(imDim * Nt, 1));
		if (store_randoms)
			RA.resize(dim_vector(imDim * Nt, 1));
		else
			RA.resize(dim_vector(1, 1));
		if (store_scatter)
			SC.resize(dim_vector(imDim * Nt, 1));
		else
			SC.resize(dim_vector(1, 1));
	}
	else {
		S.resize(dim_vector(1, 1));
		RA.resize(dim_vector(1, 1));
		SC.resize(dim_vector(1, 1));
	}
	if (store_coordinates)
		tIndex.resize(dim_vector(Nentries, 1));
	else
		tIndex.resize(dim_vector(1, 1));

	uint16_t* tIndex_p = reinterpret_cast<uint16_t*>(tIndex.fortran_vec());
	uint16_t* S_p = reinterpret_cast<uint16_t*>(S.fortran_vec());
	uint16_t* RA_p = reinterpret_cast<uint16_t*>(RA.fortran_vec());
	uint16_t* SC_p = reinterpret_cast<uint16_t*>(SC.fortran_vec());
	float* coordP = coord.fortran_vec();
	float* DcoordP = Dcoord.fortran_vec();
	uint16_t* trIndexP = reinterpret_cast<uint16_t*>(trIndex.fortran_vec());
	uint16_t* axIndexP = reinterpret_cast<uint16_t*>(axIndex.fortran_vec());
	uint16_t* DtrIndexP = reinterpret_cast<uint16_t*>(DtrIndex.fortran_vec());
	uint16_t* DaxIndexP = reinterpret_cast<uint16_t*>(DaxIndex.fortran_vec());
	uint16_t* Sino = reinterpret_cast<uint16_t*>(SinoO.fortran_vec());
	uint16_t* SinoT = reinterpret_cast<uint16_t*>(SinoOT.fortran_vec());
	uint16_t* SinoR = reinterpret_cast<uint16_t*>(SinoOR.fortran_vec());
	uint16_t* SinoC = reinterpret_cast<uint16_t*>(SinoOC.fortran_vec());
	uint16_t* SinoD = reinterpret_cast<uint16_t*>(SinoOD.fortran_vec());

	const float matlabPtr = 0;

	histogram(tmp.c_str(), tPoints_p, alku, loppu, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S_p, SC_p, RA_p, trIndexP, axIndexP, DtrIndexP, DaxIndexP, obtain_trues, store_scatter, store_randoms,
		scatter_components_p, randoms_correction, coordP, DcoordP, store_coordinates, dynamic, cryst_per_block_z, transaxial_multip, rings, sinoSize, Ndist, Nang, ringDifference, span,
		seg_p, Nt, TOFSize, nDistSide, Sino, SinoT, SinoC, SinoR, SinoD, detWPseudo, nPseudos, binSize, FWHM, verbose, nLayers, dx, dy, dz, bx, by, bz, Nx, Ny, Nz, dualLayerSubmodule, imDim, indexBased, tIndex_p, matlabPtr);



	octave_value_list retval(nargout);

	retval(0) = octave_value(SinoO);
	retval(1) = octave_value(SinoOT);
	retval(2) = octave_value(SinoOC);
	retval(3) = octave_value(SinoOR);
	retval(4) = octave_value(SinoOD);
	retval(5) = octave_value(S);
	retval(6) = octave_value(SC);
	retval(7) = octave_value(RA);
	retval(8) = octave_value(tIndex);
	retval(9) = octave_value(coord);
	retval(10) = octave_value(Dcoord);
	retval(11) = octave_value(trIndex);
	retval(12) = octave_value(axIndex);
	retval(13) = octave_value(DtrIndex);
	retval(14) = octave_value(DaxIndex);

	gROOT->Reset();

	return retval;
}
