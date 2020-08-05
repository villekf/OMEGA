/**************************************************************************
* Constructs a 3D/4D/5D sinogram from the input ring number and positions.
* This code uses the new C++ MEX API. For MATLAB R2018a and newer.
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

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "saveSinogram.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <thread>

class MexFunction : public matlab::mex::Function {

	matlab::data::ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {


		checkArguments(outputs, inputs);

		const double vali = inputs[0][0];
		const matlab::data::TypedArray<uint16_t> ringPos1 = std::move(inputs[1]);
		const matlab::data::TypedArray<uint16_t> ringPos2 = std::move(inputs[2]);
		const matlab::data::TypedArray<uint16_t> ringNumber1 = std::move(inputs[3]);
		const matlab::data::TypedArray<uint16_t> ringNumber2 = std::move(inputs[4]);
		const matlab::data::TypedArray<bool> truesIndex = std::move(inputs[5]);
		const matlab::data::TypedArray<bool> scatterIndex = std::move(inputs[6]);
		const matlab::data::TypedArray<bool> randomsIndex = std::move(inputs[7]);
		const uint64_t sinoSize = inputs[8][0];
		const uint32_t Ndist = inputs[9][0];
		const uint32_t Nang = inputs[10][0];
		const uint32_t ringDifference = inputs[11][0];
		const uint32_t span = inputs[12][0];
		const matlab::data::TypedArray<uint32_t> seg = std::move(inputs[13]);
		const uint64_t TOFSize = inputs[14][0];
		const matlab::data::TypedArray<double> time = std::move(inputs[15]);
		const uint64_t NT = inputs[16][0];
		const double alku = inputs[17][0];
		const int32_t detPerRing = inputs[18][0];
		const int32_t rings = inputs[19][0];
		const matlab::data::TypedArray<uint64_t> bins = std::move(inputs[20]);
		const int32_t nDistSide = inputs[25][0];
		const int32_t detWPseudo = inputs[26][0];
		const int32_t nPseudos = inputs[27][0];
		const int32_t crystPerBlock = inputs[28][0];

		bool storeTrues = false;
		bool storeScatter = false;
		bool storeRandoms = false;
		const int64_t koko = ringPos1.getNumberOfElements();
		const size_t pituus = seg.getNumberOfElements();
		if (truesIndex.getNumberOfElements() > 1) {
			storeTrues = true;
		}
		if (scatterIndex.getNumberOfElements() > 1) {
			storeScatter = true;
		}
		if (randomsIndex.getNumberOfElements() > 1) {
			storeRandoms = true;
		}


		/* Assign pointers to the various parameters */
		matlab::data::TypedArray<uint16_t> Sino = std::move(inputs[21]);
		matlab::data::TypedArray<uint16_t> SinoT = std::move(inputs[22]);
		matlab::data::TypedArray<uint16_t> SinoC = std::move(inputs[23]);
		matlab::data::TypedArray<uint16_t> SinoR = std::move(inputs[24]);

		openMPSino(ringPos1, ringPos2, ringNumber1, ringNumber2, truesIndex, scatterIndex, randomsIndex, sinoSize, Ndist, Nang, ringDifference,
			span, seg, time, NT, TOFSize, vali, alku, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko,
			bins, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock);

		outputs[0] = Sino;
		outputs[1] = SinoT;
		outputs[2] = SinoC;
		outputs[3] = SinoR;

	}

	void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {
		if (inputs.size() != 29) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("29 inputs required") }));
		}

		if (outputs.size() > 4) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Too many output arguments") }));
		}
	}

	void openMPSino(const matlab::data::TypedArray<uint16_t>& ringPos1, const matlab::data::TypedArray<uint16_t>& ringPos2, const matlab::data::TypedArray<uint16_t>& ringNumber1, 
		const matlab::data::TypedArray<uint16_t>& ringNumber2, const matlab::data::TypedArray<bool>& trues_index, const matlab::data::TypedArray<bool>& scatter_index, 
		const matlab::data::TypedArray<bool>& randoms_index, const uint64_t sinoSize, const uint32_t Ndist,	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, 
		const matlab::data::TypedArray<uint32_t>& seg, const matlab::data::TypedArray<double>& time, const uint64_t NT, const uint64_t TOFSize, const double vali, const double alku,
		matlab::data::TypedArray<uint16_t>& Sino, matlab::data::TypedArray<uint16_t>& SinoT, matlab::data::TypedArray<uint16_t>& SinoC, matlab::data::TypedArray<uint16_t>& SinoR,
		const bool store_trues, const bool store_scatter, const bool store_randoms, const int32_t det_per_ring, const int32_t rings,
		const int64_t koko, const matlab::data::TypedArray<uint64_t>& bins, const int32_t nDistSide, const size_t pituus, const int32_t detWPseudo, const int32_t nPseudos,
		const int32_t cryst_per_block) {

#ifdef _OPENMP
		if (omp_get_max_threads() == 1) {
			int n_threads = std::thread::hardware_concurrency();
			omp_set_num_threads(n_threads);
		}
#endif
		const bool pseudoD = detWPseudo > det_per_ring;
		const bool pseudoR = nPseudos > 0;
		int32_t gapSize = 0;
		if (pseudoR) {
			gapSize = rings / (nPseudos + 1);
		}
		uint32_t* seg_p = new uint32_t[pituus];
		for (int ll = 0; ll < pituus; ll++)
			seg_p[ll] = static_cast<uint32_t>(seg[ll]);

#pragma omp parallel for
		for (int64_t kk = 0; kk < koko; kk++) {
			double aika = 0.;
			if (NT > 1)
				double aika = time[kk];
			uint64_t binN = 0ULL;
			if (TOFSize > sinoSize)
				binN = bins[kk];
			int32_t ring_pos1 = static_cast<int32_t>(ringPos1[kk]);
			int32_t ring_pos2 = static_cast<int32_t>(ringPos2[kk]);
			int32_t ring_number1 = static_cast<int32_t>(ringNumber1[kk]);
			int32_t ring_number2 = static_cast<int32_t>(ringNumber2[kk]);
			if (pseudoD) {
				ring_pos1 += ring_pos1 / cryst_per_block;
				ring_pos2 += ring_pos2 / cryst_per_block;
			}
			if (pseudoR) {
				ring_number1 += ring_number1 / gapSize;
				ring_number2 += ring_number2 / gapSize;
			}
			const int64_t indeksi = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ring_difference, span, seg_p, aika, NT, TOFSize,
				vali, alku, detWPseudo, rings, binN, nDistSide);
			if (indeksi >= 0) {
				if (store_trues && trues_index[kk]) {
#pragma omp atomic
					SinoT[indeksi]++;
				}
				else if (store_scatter && scatter_index[kk]) {
#pragma omp atomic
					SinoC[indeksi]++;
				}
				else if (store_randoms && randoms_index[kk]) {
#pragma omp atomic
					SinoR[indeksi]++;
				}
#pragma omp atomic
				Sino[indeksi]++;
			}
		}
		delete[] seg_p;
	}

	void displayOnMATLAB(std::ostringstream& stream) {
		// Pass stream content to MATLAB fprintf function
		matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
		// Clear stream buffer
		stream.str("");
	}
};
