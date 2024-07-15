/**************************************************************************
* Constructs a 3D/4D/5D sinogram from the input ring number and positions.
* This code uses the Octave API.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
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
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#include <octave/oct.h>
#include "saveSinogram.h"
#include <thread>


DEFUN_DLD(createSinogramASCIIOct, prhs, nargout, "ASCII to sinogram help") {


	/* Check for proper number of arguments */

	uint16NDArray ringPos1 = prhs(0).uint16_array_value();
	uint16NDArray ringPos2 = prhs(1).uint16_array_value();
	uint16NDArray ringNumber1 = prhs(2).uint16_array_value();
	uint16NDArray ringNumber2 = prhs(3).uint16_array_value();
	const boolNDArray truesIndex = prhs(4).bool_array_value();
	const boolNDArray scatterIndex = prhs(5).bool_array_value();
	const boolNDArray randomsIndex = prhs(6).bool_array_value();
	const uint64_t sinoSize = prhs(7).uint64_scalar_value();
	const uint32_t Ndist = prhs(8).uint32_scalar_value();
	const uint32_t Nang = prhs(9).uint32_scalar_value();
	const uint32_t ringDifference = prhs(10).uint32_scalar_value();
	const uint32_t span = prhs(11).uint32_scalar_value();
	uint32NDArray seg = prhs(12).uint32_array_value();
	const size_t pituus = seg.numel();
	const uint64_t TOFSize = prhs(13).uint64_scalar_value();
	uint16NDArray time = prhs(14).array_value();
	const uint64_t NT = prhs(15).uint64_scalar_value();
	const int32_t detPerRing = prhs(16).int32_scalar_value();
	const int32_t rings = prhs(17).int32_scalar_value();
	uint16NDArray bins = prhs(18).uint16_array_value();
	const int32_t nDistSide = prhs(23).int32_scalar_value();
	const int32_t detWPseudo = prhs(24).int32_scalar_value();
	const int32_t nPseudos = prhs(25).int32_scalar_value();
	const int32_t crystPerBlock = prhs(26).int32_scalar_value();
	const int32_t nLayers = prhs(27).int32_scalar_value();
	uint8NDArray layers1 = prhs(28).uint8_array_value();
	uint8NDArray layers2 = prhs(29).uint8_array_value();

	bool storeTrues = false;
	bool storeScatter = false;
	bool storeRandoms = false;
	const int64_t koko = ringPos1.numel();
	if (truesIndex.numel() > 1) {
		storeTrues = true;
	}
	if (scatterIndex.numel() > 1) {
		storeScatter = true;
	}
	if (randomsIndex.numel() > 1) {
		storeRandoms = true;
	}

	const uint16_t* ring_pos1 = reinterpret_cast<uint16_t*>(ringPos1.fortran_vec());
	const uint16_t* ring_pos2 = reinterpret_cast<uint16_t*>(ringPos2.fortran_vec());
	const uint16_t* ring_number1 = reinterpret_cast<uint16_t*>(ringNumber1.fortran_vec());
	const uint16_t* ring_number2 = reinterpret_cast<uint16_t*>(ringNumber2.fortran_vec());
	const uint8_t* layer1 = reinterpret_cast<uint8_t*>(layers1.fortran_vec());
	const uint8_t* layer2 = reinterpret_cast<uint8_t*>(layers2.fortran_vec());

	const uint32_t* seg_p = reinterpret_cast<uint32_t*>(seg.fortran_vec());
	const uint16_t* bins_p = reinterpret_cast<uint16_t*>(bins.fortran_vec());
	const uint16_t* tPoint = reinterpret_cast<uint16_t*>(time.fortran_vec());
	const bool* tIndex = truesIndex.fortran_vec();
	const bool* cIndex = scatterIndex.fortran_vec();
	const bool* rIndex = randomsIndex.fortran_vec();


	/* Assign pointers to the various parameters */
	uint16NDArray SinoO = prhs(19).uint16_array_value();
	uint16NDArray SinoTO = prhs(20).uint16_array_value();
	uint16NDArray SinoCO = prhs(21).uint16_array_value();
	uint16NDArray SinoRO = prhs(22).uint16_array_value();
	uint16_t* Sino = reinterpret_cast<uint16_t*>(SinoO.fortran_vec());
	uint16_t* SinoT = reinterpret_cast<uint16_t*>(SinoTO.fortran_vec());
	uint16_t* SinoC = reinterpret_cast<uint16_t*>(SinoCO.fortran_vec());
	uint16_t* SinoR = reinterpret_cast<uint16_t*>(SinoRO.fortran_vec());

	openMPSino(ring_pos1, ring_pos2, ring_number1, ring_number2, tIndex, cIndex, rIndex, sinoSize, Ndist, Nang, ringDifference,
		span, seg_p, tPoint, NT, TOFSize, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko,
		bins_p, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock, layer1, layer2, nLayers);

	octave_value_list retval(nargout);

	retval(0) = octave_value(SinoO);
	retval(1) = octave_value(SinoTO);
	retval(2) = octave_value(SinoCO);
	retval(3) = octave_value(SinoRO);

	return retval;

}
