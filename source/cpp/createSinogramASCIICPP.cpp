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
#define MATLAB
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "saveSinogram.h"
#include <type_traits>

//! Extracts the pointer to underlying data from the non-const iterator (`TypedIterator<T>`).
/*! This function does not throw any exceptions. */
template <typename T>
inline T* toPointer(const matlab::data::TypedIterator<T>& it) MW_NOEXCEPT {
	static_assert(std::is_arithmetic<T>::value && !std::is_const<T>::value,
		"Template argument T must be a std::is_arithmetic and non-const type.");
	return it.operator->();
}
/*! Extracts pointer to the first element in the array.
 *  Example usage:
 *  \code
 *  ArrayFactory factory;
 *  TypedArray<double> A = factory.createArray<double>({ 2,2 }, { 1.0, 3.0, 2.0, 4.0 });
 *  auto ptr = getPointer(A);
 *  \endcode
 *  \note Do not call `getPointer` with temporary object. e.g., the following code is ill-formed.
 *        auto ptr=getPointer(factory.createArray<double>({ 2,2 },{ 1.0, 3.0, 2.0, 4.0 }));
 */
template <typename T>
inline T* getPointer(matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
	static_assert(std::is_arithmetic<T>::value, "Template argument T must be a std::is_arithmetic type.");
	return toPointer(arr.begin());
}
template <typename T>
inline const T* getPointer(const matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
	return getPointer(const_cast<matlab::data::TypedArray<T>&>(arr));
}

class MexFunction : public matlab::mex::Function {

	matlab::data::ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {


		checkArguments(outputs, inputs);

		const matlab::data::TypedArray<uint16_t> ringPos1 = inputs[0];
		const matlab::data::TypedArray<uint16_t> ringPos2 = inputs[1];
		const matlab::data::TypedArray<uint16_t> ringNumber1 = inputs[2];
		const matlab::data::TypedArray<uint16_t> ringNumber2 = inputs[3];
		const matlab::data::TypedArray<bool> truesIndex = inputs[4];
		const matlab::data::TypedArray<bool> scatterIndex = inputs[5];
		const matlab::data::TypedArray<bool> randomsIndex = inputs[6];
		const uint64_t sinoSize = inputs[7][0];
		const uint32_t Ndist = inputs[8][0];
		const uint32_t Nang = inputs[9][0];
		const uint32_t ringDifference = inputs[10][0];
		const uint32_t span = inputs[11][0];
		const matlab::data::TypedArray<uint32_t> seg = inputs[12];
		const uint64_t TOFSize = inputs[13][0];
		const matlab::data::TypedArray<uint16_t> tIndex = inputs[14];
		const uint64_t NT = inputs[15][0];
		const int32_t detPerRing = inputs[16][0];
		const int32_t rings = inputs[17][0];
		const matlab::data::TypedArray<uint16_t> bins = inputs[18];
		const int32_t nDistSide = inputs[23][0];
		const int32_t detWPseudo = inputs[24][0];
		const int32_t nPseudos = inputs[25][0];
		const int32_t crystPerBlock = inputs[26][0];
		const int32_t nLayers = inputs[27][0];
		const matlab::data::TypedArray<uint8_t> layer1 = inputs[28];
		const matlab::data::TypedArray<uint8_t> layer2 = inputs[29];

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
		matlab::data::TypedArray<uint16_t> Sino = std::move(inputs[19]);
		matlab::data::TypedArray<uint16_t> SinoT = std::move(inputs[20]);
		matlab::data::TypedArray<uint16_t> SinoC = std::move(inputs[21]);
		matlab::data::TypedArray<uint16_t> SinoR = std::move(inputs[22]);

		uint16_t* SinoP = getPointer(Sino);
		uint16_t* SinoTP = getPointer(SinoT);
		uint16_t* SinoCP = getPointer(SinoC);
		uint16_t* SinoRP = getPointer(SinoR);
		const uint16_t* tIndP = getPointer(tIndex);
		const uint16_t* ringPos1P = getPointer(ringPos1);
		const uint16_t* ringPos2P = getPointer(ringPos2);
		const uint16_t* ringNumber1P = getPointer(ringNumber1);
		const uint16_t* ringNumber2P = getPointer(ringNumber2);
		const uint16_t* binsP = getPointer(bins);
		const uint32_t* segP = getPointer(seg);
		const bool* truesIndexP = getPointer(truesIndex);
		const bool* scatterIndexP = getPointer(scatterIndex);
		const bool* randomsIndexP = getPointer(randomsIndex);
		const uint8_t* layer1P = getPointer(layer1);
		const uint8_t* layer2P = getPointer(layer2);

		openMPSino(ringPos1P, ringPos2P, ringNumber1P, ringNumber2P, truesIndexP, scatterIndexP, randomsIndexP, sinoSize, Ndist, Nang, ringDifference,
			span, segP, tIndP, NT, TOFSize, SinoP, SinoTP, SinoCP, SinoRP, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko,
			binsP, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock, layer1P, layer2P, nLayers);

		outputs[0] = Sino;
		outputs[1] = SinoT;
		outputs[2] = SinoC;
		outputs[3] = SinoR;

	}

	void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {
		if (inputs.size() != 30) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("30 inputs required") }));
		}

		if (outputs.size() > 4) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Too many output arguments") }));
		}
	}

	void displayOnMATLAB(std::ostringstream& stream) {
		// Pass stream content to MATLAB fprintf function
		matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
		// Clear stream buffer
		stream.str("");
	}
};
