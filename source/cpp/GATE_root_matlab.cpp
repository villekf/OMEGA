/**************************************************************************
* ROOT file import into MATLAB. This file contains the C++ implementation.
* Requires MATLAB 2019a or later.
*
* Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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
#include "rootImport.h"
#include "TChain.h"
#include <type_traits>

// source: https://se.mathworks.com/matlabcentral/answers/436916-how-to-access-raw-data-via-the-mex-c-api#answer_415645
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
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = matlab::mex::Function::getEngine();
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		checkArguments(outputs, inputs);

		int ind = 1;
		const matlab::data::TypedArray<double> tPoints = inputs[ind++];
		const double alku = inputs[ind++][0];
		const double loppu = inputs[ind++][0];
		const matlab::data::TypedArray<uint32_t> detectorsP = std::move(inputs[ind++]);
		const uint32_t blocks_per_ring = inputs[ind++][0];
		const matlab::data::TypedArray<uint32_t> cryst_per_blockP = std::move(inputs[ind++]);
		const matlab::data::TypedArray<uint32_t> det_per_ringP = std::move(inputs[ind++]);
		const uint32_t linear_multp = inputs[ind++][0];
		bool source = std::move(inputs[ind++][0]);
		const int64_t Nt = inputs[ind++][0];
		bool obtain_trues = std::move(inputs[ind++][0]);
		bool store_scatter = std::move(inputs[ind++][0]);
		bool store_randoms = std::move(inputs[ind++][0]);
		matlab::data::TypedArray<uint8_t> scatter_components = std::move(inputs[ind++]);
		const bool randoms_correction = inputs[ind++][0];
		bool store_coordinates = std::move(inputs[ind++][0]);
		const uint32_t transaxial_multip = inputs[ind++][0];
		const matlab::data::TypedArray<uint32_t> cryst_per_block_zP = std::move(inputs[ind++]);
		const matlab::data::TypedArray<uint32_t> ringsP = std::move(inputs[ind++]);
		const matlab::data::TypedArray<uint64_t> sinoSizeP = std::move(inputs[ind++]);
		const bool TOF = inputs[ind++][0];
		const bool verbose = inputs[ind++][0];
		const int32_t nLayers = inputs[ind++][0];
		const uint32_t Ndist = inputs[ind++][0];
		const matlab::data::TypedArray<uint32_t> NangP = std::move(inputs[ind++]);
		const uint32_t ringDifference = inputs[ind++][0];
		const uint32_t span = inputs[ind++][0];
		const matlab::data::TypedArray<uint32_t> seg = inputs[ind++];
		const uint64_t TOFSize = inputs[ind++][0];
		const int32_t nDistSide = inputs[ind++][0];
		matlab::data::TypedArray<uint16_t> Sino = std::move(inputs[ind++]);
		matlab::data::TypedArray<uint16_t> SinoT = std::move(inputs[ind++]);
		matlab::data::TypedArray<uint16_t> SinoC = std::move(inputs[ind++]);
		matlab::data::TypedArray<uint16_t> SinoR = std::move(inputs[ind++]);
		matlab::data::TypedArray<uint16_t> SinoD = std::move(inputs[ind++]);
		const matlab::data::TypedArray<uint32_t> detWPseudoP = std::move(inputs[ind++]);
		const int32_t nPseudos = inputs[ind++][0];
		const double binSize = inputs[ind++][0];
		const double FWHM = inputs[ind++][0];
		const float dx = inputs[ind++][0];
		const float dy = inputs[ind++][0];
		const float dz = inputs[ind++][0];
		const float bx = inputs[ind++][0];
		const float by = inputs[ind++][0];
		const float bz = inputs[ind++][0];
		const int64_t Nx = inputs[ind++][0];
		const int64_t Ny = inputs[ind++][0];
		const int64_t Nz = inputs[ind++][0];
		const bool dualLayerSubmodule = inputs[ind++][0];
		const bool indexBased = inputs[ind++][0];

		const int64_t imDim = Nx * Ny * Nz;

		const bool dynamic = Nt > 1;


		matlab::data::CharArray newName(inputs[0]);
		std::string newPropertyValue = newName.toAscii();
		const char* argv = newPropertyValue.c_str();

		TChain* Coincidences = new TChain("Coincidences");
		Coincidences->Add(argv);
		size_t Nentries = Coincidences->GetEntries();

		// Output data
		matlab::data::TypedArray<float> coord = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> Dcoord = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> trIndex = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> axIndex = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> DtrIndex = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> DaxIndex = factory.createArray<uint16_t>({ 1, 1 });

		TChain* delay = nullptr;
		int64_t Ndelays = 0LL;

		if (randoms_correction) {
			delay = new TChain("delay");
			delay->Add(argv);
			Ndelays = delay->GetEntries();
		}

		delete Coincidences;
		if (randoms_correction)
			delete delay;

		if (store_coordinates) {
			coord = factory.createArray<float>({ 6, Nentries });
			std::fill(coord.begin(), coord.end(), 0.f);
			if (randoms_correction) {
				Dcoord = factory.createArray<float>({ 6, Nentries });
				std::fill(Dcoord.begin(), Dcoord.end(), 0.f);
			}
		}
		if (indexBased) {
			trIndex = factory.createArray<uint16_t>({ 2, Nentries });
			std::fill(trIndex.begin(), trIndex.end(), 0.f);
			axIndex = factory.createArray<uint16_t>({ 2, Nentries });
			std::fill(axIndex.begin(), axIndex.end(), 0.f);
			if (randoms_correction) {
				DtrIndex = factory.createArray<uint16_t>({ 2, Nentries });
				std::fill(DtrIndex.begin(), DtrIndex.end(), 0.f);
				DaxIndex = factory.createArray<uint16_t>({ 2, Nentries });
				std::fill(DaxIndex.begin(), DaxIndex.end(), 0.f);
			}
		}

		matlab::data::TypedArray<uint16_t> tIndex = factory.createArray<uint16_t>({ 1, 1 });
		if (dynamic) {
			tIndex = factory.createArray<uint16_t>({ Nentries, 1 });
			std::fill(tIndex.begin(), tIndex.end(), static_cast<uint16_t>(32768));
		}
		matlab::data::TypedArray<uint16_t> S = factory.createArray<uint16_t>({ 1, 1 }, { 0 });
		matlab::data::TypedArray<uint16_t> SC = factory.createArray<uint16_t>({ 1, 1 }, { 0 });
		matlab::data::TypedArray<uint16_t> RA = factory.createArray<uint16_t>({ 1, 1 }, { 0 });

		if (source) {
			S = factory.createArray<uint16_t>({ static_cast<uint64_t>(imDim * Nt), 1 });
			std::fill(S.begin(), S.end(), static_cast<uint16_t>(0));
			if (store_scatter) {
				SC = factory.createArray<uint16_t>({ static_cast<uint64_t>(imDim * Nt), 1 });
				std::fill(SC.begin(), SC.end(), static_cast<uint16_t>(0));
			}
			if (store_randoms) {
				RA = factory.createArray<uint16_t>({ static_cast<uint64_t>(imDim * Nt), 1 });
				std::fill(RA.begin(), RA.end(), static_cast<uint16_t>(0));
			}
		}

		const double* tPointer = getPointer(tPoints);
		uint16_t* SPointer = getPointer(S);
		uint16_t* SCPointer = getPointer(SC);
		uint16_t* RAPointer = getPointer(RA);
		uint8_t* scatterPointer = getPointer(scatter_components);
		float* coordP = getPointer(coord);
		float* DcoordP = getPointer(Dcoord);
		uint16_t* trIndexPtr = getPointer(trIndex);
		uint16_t* axIndexPtr = getPointer(axIndex);
		uint16_t* DtrIndexPtr = getPointer(DtrIndex);
		uint16_t* DaxIndexPtr = getPointer(DaxIndex);
		const uint32_t* segP = getPointer(seg);
		const uint32_t* detectors = getPointer(detectorsP);
		const uint32_t* cryst_per_block = getPointer(cryst_per_blockP);
		const uint32_t* det_per_ring = getPointer(det_per_ringP);
		const uint32_t* cryst_per_block_z = getPointer(cryst_per_block_zP);
		const uint32_t* rings = getPointer(ringsP);
		const uint64_t* sinoSize = getPointer(sinoSizeP);
		const uint32_t* Nang = getPointer(NangP);
		const uint32_t* detWPseudo = getPointer(detWPseudoP);
		uint16_t* SinoP = getPointer(Sino);
		uint16_t* SinoTP = getPointer(SinoT);
		uint16_t* SinoCP = getPointer(SinoC);
		uint16_t* SinoRP = getPointer(SinoR);
		uint16_t* SinoDP = getPointer(SinoD);
		uint16_t* tIndP = getPointer(tIndex);


		histogram(argv, tPointer, alku, loppu, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, SPointer, SCPointer, RAPointer, trIndexPtr, axIndexPtr, DtrIndexPtr, DaxIndexPtr, obtain_trues, store_scatter, store_randoms,
			scatterPointer, randoms_correction, coordP, DcoordP, store_coordinates, dynamic, cryst_per_block_z, transaxial_multip, rings, sinoSize, Ndist, Nang, ringDifference, span,
			segP, Nt, TOFSize, nDistSide, SinoP, SinoTP, SinoCP, SinoRP, SinoDP, detWPseudo, nPseudos, binSize, FWHM, verbose, nLayers, dx, dy, dz, bx, by, bz, Nx, Ny, Nz, dualLayerSubmodule, imDim, indexBased, tIndP, matlabPtr);


		outputs[0] = std::move(Sino);
		outputs[1] = std::move(SinoT);
		outputs[2] = std::move(SinoC);
		outputs[3] = std::move(SinoR);
		outputs[4] = std::move(SinoD);
		outputs[5] = std::move(S);
		outputs[6] = std::move(SC);
		outputs[7] = std::move(RA);
		outputs[8] = std::move(tIndex);
		outputs[9] = std::move(coord);
		outputs[10] = std::move(Dcoord);
		outputs[11] = std::move(trIndex);
		outputs[12] = std::move(axIndex);
		outputs[13] = std::move(DtrIndex);
		outputs[14] = std::move(DaxIndex);

	}

	void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {
		if (inputs.size() != 51) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("51 inputs required") }));
		}

		if (outputs.size() > 15) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Too many output arguments") }));
		}

		if (inputs[0].getType() != matlab::data::ArrayType::CHAR) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("First input must be char") }));
		}
	}

};
