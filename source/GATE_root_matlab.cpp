/**************************************************************************
* ROOT file import into MATLAB. This file contains the C++ implementation.
* Requires MATLAB 2019a or later.
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
#include <stdint.h>
#include <cmath>
#include <algorithm>
#include <TROOT.h>
#include "TChain.h"


class MexFunction : public matlab::mex::Function {

	matlab::data::ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		checkArguments(outputs, inputs);

		const double vali = inputs[1][0];
		const double alku = inputs[2][0];
		const double loppu = inputs[3][0];
		const uint32_t detectors = inputs[4][0];
		const uint32_t blocks_per_ring = inputs[5][0];
		const uint32_t cryst_per_block = inputs[6][0];
		const uint32_t det_per_ring = inputs[7][0];
		const uint32_t linear_multp = inputs[8][0];
		const bool source = inputs[9][0];
		const matlab::data::TypedArray<double> time_intervals = std::move(inputs[10]);
		const bool obtain_trues = inputs[11][0];
		const bool store_scatter = inputs[12][0];
		const bool store_randoms = inputs[13][0];
		const matlab::data::TypedArray<bool> scatter_components = std::move(inputs[14]);
		const bool randoms_correction = inputs[15][0];
		const bool store_coordinates = inputs[16][0];
		size_t outsize2 = (loppu - alku) / vali;
		bool large_case = false;
		if (static_cast<uint64_t>(detectors) * static_cast<uint64_t>(detectors) * 2ULL >= 1024ULL * 1024ULL * 1024ULL * 2ULL) {
			large_case = true;
		}

		bool dynamic = false;


		matlab::data::CharArray newName(inputs[0]);
		std::string newPropertyValue = newName.toAscii();
		const char* argv = newPropertyValue.c_str();


		TChain* Coincidences = new TChain("Coincidences");
		Coincidences->Add(argv);

		TChain* delay;
		size_t Ndelays = 0ULL;

		if (randoms_correction) {
			delay = new TChain("delay");
			delay->Add(argv);
			Ndelays = delay->GetEntries();
		}


		size_t Nentries = Coincidences->GetEntries();

		// Output data
		matlab::data::TypedArray<uint16_t> LL1 = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> LL2 = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> Ltrues = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> Lrandoms = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> Lscatter = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> Ldelay1 = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<uint16_t> Ldelay2 = factory.createArray<uint16_t>({ 1, 1 });
		matlab::data::TypedArray<float> x1 = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> x2 = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> y1 = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> y2 = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> z1 = factory.createArray<float>({ 1, 1 });
		matlab::data::TypedArray<float> z2 = factory.createArray<float>({ 1, 1 });

		if (outsize2 == 1ULL && !large_case) {
			if (obtain_trues) {
				Ltrues = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Ltrues.begin(), Ltrues.end(), static_cast<uint16_t>(0));
			}
			else {
				//plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
				//matlab::data::TypedArray<uint16_t> Ltrues = factory.createArray<uint16_t>({ 1, 1 });
				//matlab::data::TypedArray<uint16_t> Ltrues = factory.createArrayFromBuffer<uint16_t>({ 1, 1 }, true_buf);
			}
			if (store_randoms) {
				//plhs[6] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
				Lrandoms = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Lrandoms.begin(), Lrandoms.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	Lrandoms = factory.createArray<uint16_t>({ 1, 1 });
			//}
			if (store_scatter) {
				//plhs[7] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
				Lscatter = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Lscatter.begin(), Lscatter.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	Lscatter = factory.createArray<uint16_t>({ 1, 1 });
			//}
			if (randoms_correction) {
				//plhs[9] = mxCreateNumericMatrix(detectors, detectors, mxUINT16_CLASS, mxREAL);
				//plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
				Ldelay1 = factory.createArray<uint16_t>({ detectors, detectors });
				Ldelay2 = factory.createArray<uint16_t>({ 1, 1 }, { static_cast<uint16_t>(0) });
				std::fill(Ldelay1.begin(), Ldelay1.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[9] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	//plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	Ldelay1 = factory.createArray<uint16_t>({ 1, 1 });
			//	Ldelay2 = factory.createArray<uint16_t>({ 1, 1 });
			//}
			LL1 = factory.createArray<uint16_t>({ detectors, detectors });
			std::fill(LL1.begin(), LL1.end(), static_cast<uint16_t>(0));
		}
		else {
			//plhs[0] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
			//plhs[1] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
			LL1 = factory.createArray<uint16_t>({ Nentries, 1 });
			LL2 = factory.createArray<uint16_t>({ Nentries, 1 });
			std::fill(LL1.begin(), LL1.end(), static_cast<uint16_t>(0));
			std::fill(LL2.begin(), LL2.end(), static_cast<uint16_t>(0));

			if (obtain_trues) {
			//	//plhs[5] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
				Ltrues = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Ltrues.begin(), Ltrues.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	matlab::data::TypedArray<uint16_t> Ltrues = factory.createArray<uint16_t>({ 1, 1 });
			//}
			if (store_randoms) {
				//plhs[6] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
				Lrandoms = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Lrandoms.begin(), Lrandoms.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	matlab::data::TypedArray<uint16_t> Lrandoms = factory.createArray<uint16_t>({ 1, 1 });
			//}
			if (store_scatter) {
				//plhs[7] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
				Lscatter = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Lscatter.begin(), Lscatter.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	matlab::data::TypedArray<uint16_t> Lscatter = factory.createArray<uint16_t>({ 1, 1 });
			//}
			if (randoms_correction) {
				//plhs[9] = mxCreateNumericMatrix(Ndelays, 1, mxUINT16_CLASS, mxREAL);
				//plhs[10] = mxCreateNumericMatrix(Ndelays, 1, mxUINT16_CLASS, mxREAL);
				Ldelay1 = factory.createArray<uint16_t>({ Ndelays, 1 });
				Ldelay2 = factory.createArray<uint16_t>({ Ndelays, 1 });
				std::fill(Ldelay1.begin(), Ldelay1.end(), static_cast<uint16_t>(0));
				std::fill(Ldelay2.begin(), Ldelay2.end(), static_cast<uint16_t>(0));
			}
			//else {
			//	//plhs[9] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	//plhs[10] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			//	matlab::data::TypedArray<uint16_t> Ldelay1 = factory.createArray<uint16_t>({ 1, 1 });
			//	matlab::data::TypedArray<uint16_t> Ldelay2 = factory.createArray<uint16_t>({ 1, 1 });
			//}
			dynamic = true;
		}
		if (store_coordinates) {
			x1 = factory.createArray<float>({ Nentries, 1 });
			std::fill(x1.begin(), x1.end(), 0.f);
			x2 = factory.createArray<float>({ Nentries, 1 });
			std::fill(x2.begin(), x2.end(), 0.f);
			y1 = factory.createArray<float>({ Nentries, 1 });
			std::fill(y1.begin(), y1.end(), 0.f);
			y2 = factory.createArray<float>({ Nentries, 1 });
			std::fill(y2.begin(), y2.end(), 0.f);
			z1 = factory.createArray<float>({ Nentries, 1 });
			std::fill(z1.begin(), z1.end(), 0.f);
			z2 = factory.createArray<float>({ Nentries, 1 });
			std::fill(z2.begin(), z2.end(), 0.f);
		}
		//plhs[2] = mxCreateNumericMatrix(outsize2 + 2, 1, mxUINT32_CLASS, mxREAL);
		matlab::data::TypedArray<uint32_t> tpoints = factory.createArray<uint32_t>({ outsize2 + 2, 1 });
		std::fill(tpoints.begin(), tpoints.end(), 0U);
		matlab::data::TypedArray<uint32_t> tpoints_delay = factory.createArray<uint32_t>({ 1, 1 }, { 0 });
		matlab::data::TypedArray<float> S = factory.createArray<float>({ 1, 1 }, { 0.f });
		matlab::data::TypedArray<uint8_t> trues_loc = factory.createArray<uint8_t>({ 1, 1 }, { static_cast<uint8_t>(0) });
		matlab::data::TypedArray<uint8_t> randoms_loc = factory.createArray<uint8_t>({ 1, 1 }, { static_cast<uint8_t>(0) });
		matlab::data::TypedArray<uint8_t> scatter_loc = factory.createArray<uint8_t>({ 1, 1 }, { static_cast<uint8_t>(0) });
		if (randoms_correction) {
			//plhs[12] = mxCreateNumericMatrix(outsize2 + 2, 1, mxUINT32_CLASS, mxREAL);
			tpoints_delay = factory.createArray<uint32_t>({ outsize2 + 2, 1 });
			std::fill(tpoints_delay.begin(), tpoints_delay.end(), 0U);
		}
		//else {
		//	//plhs[12] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
		//	matlab::data::TypedArray<uint32_t> tpoints_delay = factory.createArray<uint32_t>({ 1, 1 });
		//}

		if (source) {
			//plhs[3] = mxCreateNumericMatrix(Nentries, 6, mxSINGLE_CLASS, mxREAL);
			S = factory.createArray<float>({ Nentries, 6 });
			std::fill(S.begin(), S.end(), 0.f);
		}
		//else {
		//	//plhs[3] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//	matlab::data::TypedArray<float> S = factory.createArray<float>({ 1, 1 });
		//}
		//plhs[4] = mxCreateNumericMatrix(2, 1, mxINT32_CLASS, mxREAL);
		matlab::data::TypedArray<int32_t> int_loc = factory.createArray<int32_t>({ 2, 1 }, { 0, 0 });

		//plhs[11] = mxCreateNumericMatrix(2, 1, mxINT32_CLASS, mxREAL);
		matlab::data::TypedArray<int32_t> int_loc_delay = factory.createArray<int32_t>({ 2, 1 }, { 0, 0 });

		if (obtain_trues && source && outsize2 == 1ULL) {
			//plhs[8] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
			//trues_loc = (bool*)mxGetData(plhs[8]);
			trues_loc = factory.createArray<uint8_t>({ Nentries, 1 });
			std::fill(trues_loc.begin(), trues_loc.end(), static_cast<uint8_t>(0));
		}
		//else {
		//	//trues_loc = 0;
		//	//plhs[8] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
		//	matlab::data::TypedArray<bool> trues_loc = factory.createArray<bool>({ 1, 1 });
		//}
		if (store_randoms && source && outsize2 == 1ULL) {
			//plhs[13] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
			//randoms_loc = (bool*)mxGetData(plhs[13]);
			randoms_loc = factory.createArray<uint8_t>({ Nentries, 1 });
			std::fill(randoms_loc.begin(), randoms_loc.end(), static_cast<uint8_t>(0));
		}
		//else {
		//	//randoms_loc = 0;
		//	//plhs[13] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
		//	matlab::data::TypedArray<bool> randoms_loc = factory.createArray<bool>({ 1, 1 });
		//}
		if (store_scatter && source && outsize2 == 1) {
			//plhs[14] = mxCreateNumericMatrix(Nentries, 1, mxLOGICAL_CLASS, mxREAL);
			//scatter_loc = (bool*)mxGetData(plhs[14]);
			scatter_loc = factory.createArray<uint8_t>({ Nentries, 1 });
			std::fill(scatter_loc.begin(), scatter_loc.end(), static_cast<uint8_t>(0));
		}
		//else {
		//	//scatter_loc = 0;
		//	//plhs[14] = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
		//	matlab::data::TypedArray<bool> scatter_loc = factory.createArray<bool>({ 1, 1 });
		//}

		histogram(LL1, LL2, tpoints, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S,
			Coincidences, Nentries, time_intervals, int_loc, obtain_trues, store_scatter, store_randoms, scatter_components, Ltrues, Lscatter,
			Lrandoms, trues_loc, Ndelays, randoms_correction, delay, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_loc, scatter_loc, x1, x2, 
			y1, y2, z1, z2, store_coordinates, dynamic, large_case);


		//std::ostringstream stream;
		//stream << "Data loaded" << std::endl;
		//displayOnMATLAB(stream);

		delete Coincidences;
		if (randoms_correction)
			delete delay;

		outputs[0] = std::move(LL1);
		outputs[1] = std::move(LL2);
		outputs[2] = std::move(tpoints);
		outputs[3] = std::move(S);
		outputs[4] = std::move(int_loc);
		outputs[5] = std::move(Ltrues);
		outputs[6] = std::move(Lscatter);
		outputs[7] = std::move(Lrandoms);
		outputs[8] = std::move(trues_loc);
		outputs[9] = std::move(Ldelay1);
		outputs[10] = std::move(Ldelay2);
		outputs[11] = std::move(int_loc_delay);
		outputs[12] = std::move(tpoints_delay);
		outputs[13] = std::move(randoms_loc);
		outputs[14] = std::move(scatter_loc);
		outputs[15] = std::move(x1);
		outputs[16] = std::move(x2);
		outputs[17] = std::move(y1);
		outputs[18] = std::move(y2);
		outputs[19] = std::move(z1);
		outputs[20] = std::move(z2);

		//mexEvalString("pause(.001);");

		//std::ostringstream stream1;
		//stream1 << "Outputs formed" << std::endl;
		//displayOnMATLAB(stream1);
		//gROOT->Reset();

	}

	void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {
		if (inputs.size() != 17) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("17 inputs required") }));
		}

		if (outputs.size() > 21) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Too many output arguments") }));
		}

		if (inputs[0].getType() != matlab::data::ArrayType::CHAR) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("First input must be char") }));
		}
	}

	void histogram(matlab::data::TypedArray<uint16_t>& LL1, matlab::data::TypedArray<uint16_t>& LL2, matlab::data::TypedArray<uint32_t>& tpoints, double vali, const double alku, const double loppu, const size_t outsize2,
		const uint32_t detectors, bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring,
		const uint32_t det_per_ring, matlab::data::TypedArray<float>& S, TTree* Coincidences, const size_t Nentries, const matlab::data::TypedArray<double>& time_intervals, matlab::data::TypedArray<int32_t>& int_loc,
		bool obtain_trues, bool store_scatter, bool store_randoms, const matlab::data::TypedArray<bool> scatter_components, matlab::data::TypedArray<uint16_t>& Ltrues, matlab::data::TypedArray<uint16_t>& Lscatter,
		matlab::data::TypedArray<uint16_t>& Lrandoms, matlab::data::TypedArray<uint8_t>& trues_loc, const size_t Ndelays, bool randoms_correction, TTree* delay, matlab::data::TypedArray<uint16_t>& Ldelay1, matlab::data::TypedArray<uint16_t>& Ldelay2,
		matlab::data::TypedArray<int32_t>& int_loc_delay, matlab::data::TypedArray<uint32_t>& tpoints_delay, matlab::data::TypedArray<uint8_t>& randoms_loc, matlab::data::TypedArray<uint8_t>& scatter_loc,
		matlab::data::TypedArray<float>& x1, matlab::data::TypedArray<float>& x2, matlab::data::TypedArray<float>& y1, matlab::data::TypedArray<float>& y2,
		matlab::data::TypedArray<float>& z1, matlab::data::TypedArray<float>& z2, bool store_coordinates, const bool dynamic, const bool large_case)
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
			std::ostringstream stream;
			stream << "No crystal location information was found from file. Aborting." << std::endl;
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
				std::ostringstream stream;
				stream << "No X-source coordinates saved for first photon, unable to save source coordinates" << std::endl;
				source = false;
			}
			if (Coincidences->GetBranchStatus("sourcePosX2"))
				Coincidences->SetBranchAddress("sourcePosX2", &sourcePosX2);
			else {
				std::ostringstream stream;
				stream << "No X-source coordinates saved for second photon, unable to save source coordinates" << std::endl;
				source = false;
			}
			if (Coincidences->GetBranchStatus("sourcePosY1"))
				Coincidences->SetBranchAddress("sourcePosY1", &sourcePosY1);
			else {
				std::ostringstream stream;
				stream << "No Y-source coordinates saved for first photon, unable to save source coordinates" << std::endl;
				source = false;
			}
			if (Coincidences->GetBranchStatus("sourcePosY2"))
				Coincidences->SetBranchAddress("sourcePosY2", &sourcePosY2);
			else {
				std::ostringstream stream;
				stream << "No Y-source coordinates saved for second photon, unable to save source coordinates" << std::endl;
				source = false;
			}
			if (Coincidences->GetBranchStatus("sourcePosZ1"))
				Coincidences->SetBranchAddress("sourcePosZ1", &sourcePosZ1);
			else {
				std::ostringstream stream;
				stream << "No Z-source coordinates saved for first photon, unable to save source coordinates" << std::endl;
				source = false;
			}
			if (Coincidences->GetBranchStatus("sourcePosZ2"))
				Coincidences->SetBranchAddress("sourcePosZ2", &sourcePosZ2);
			else {
				std::ostringstream stream;
				stream << "No Z-source coordinates saved for second photon, unable to save source coordinates" << std::endl;
				source = false;
			}
		}
		//Coincidences->SetBranchAddress("time1", &time1);
		if (Coincidences->GetBranchStatus("time2"))
			Coincidences->SetBranchAddress("time2", &time2);
		else {
			if (dynamic) {
				std::ostringstream stream;
				stream << "Dynamic examination selected, but no time information was found from file. Aborting." << std::endl;
				return;
			}
			no_time = true;
		}
		if (store_coordinates) {
			if (Coincidences->GetBranchStatus("globalPosX1"))
				Coincidences->SetBranchAddress("globalPosX1", &globalPosX1);
			else {
				std::ostringstream stream;
				stream << "No X-source coordinates saved for first photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
			if (Coincidences->GetBranchStatus("globalPosX2"))
				Coincidences->SetBranchAddress("globalPosX2", &globalPosX2);
			else {
				std::ostringstream stream;
				stream << "No X-source coordinates saved for second photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
			if (Coincidences->GetBranchStatus("globalPosY1"))
				Coincidences->SetBranchAddress("globalPosY1", &globalPosY1);
			else {
				std::ostringstream stream;
				stream << "No Y-source coordinates saved for first photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
			if (Coincidences->GetBranchStatus("globalPosY2"))
				Coincidences->SetBranchAddress("globalPosY2", &globalPosY2);
			else {
				std::ostringstream stream;
				stream << "No Y-source coordinates saved for second photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
			if (Coincidences->GetBranchStatus("globalPosZ1"))
				Coincidences->SetBranchAddress("globalPosZ1", &globalPosZ1);
			else {
				std::ostringstream stream;
				stream << "No Z-source coordinates saved for first photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
			if (Coincidences->GetBranchStatus("globalPosZ2"))
				Coincidences->SetBranchAddress("globalPosZ2", &globalPosZ2);
			else {
				std::ostringstream stream;
				stream << "No Z-source coordinates saved for second photon interaction, unable to save interaction coordinates" << std::endl;
				store_coordinates = false;
			}
		}
		if (obtain_trues || store_scatter || store_randoms) {
			if (Coincidences->GetBranchStatus("eventID1"))
				Coincidences->SetBranchAddress("eventID1", &eventID1);
			else {
				std::ostringstream stream;
				stream << "No event IDs saved for first photon, unable to save trues/scatter/randoms" << std::endl;
				obtain_trues = false;
				store_scatter = false;
				store_randoms = false;
			}
			if (Coincidences->GetBranchStatus("eventID2"))
				Coincidences->SetBranchAddress("eventID2", &eventID2);
			else {
				std::ostringstream stream;
				stream << "No event IDs saved for second photon, unable to save trues/scatter/randoms" << std::endl;
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

			if (store_scatter && any == 2 && scatter_components[0]) {
				std::ostringstream stream;
				stream << "Compton phantom selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
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

			if (store_scatter && ((any == 4 && next == 1) || (any == 2 && next == 0)) && scatter_components[1]) {
				std::ostringstream stream;
				stream << "Compton crystal selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
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

			if (store_scatter && ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1)) && scatter_components[2]) {
				std::ostringstream stream;
				stream << "Rayleigh phantom selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
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

			if (store_scatter && ((any == 8 && next == 3) || (any == 2 && next == 0) || (any == 4 && next == 1) || (any == 6 && next == 2)) && scatter_components[3]) {
				std::ostringstream stream;
				stream << "Rayleigh crystal selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
			}


			if (store_scatter && any == 8) {
				std::ostringstream stream;
				stream << "Store scatter selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
				store_scatter = false;
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
		float ff = 0.01f;
		uint64_t NN = Nentries * ff;

		for (uint64_t kk = 0; kk < Nentries; kk++) {

			jj++;
			nbytes += Coincidences->GetEntry(kk);

			//if (kk == NN) {
			//	ff += 0.01f;
			//	NN = Nentries * ff; 
			//	std::ostringstream stream;
			//	stream << static_cast<uint32_t>(ff * 100.f) << "%" << std::endl;
			//	displayOnMATLAB(stream);
			//}
			//if (kk == 0) {
			//	std::ostringstream stream;
			//	stream << "ROOT file load started" << std::endl;
			//	displayOnMATLAB(stream);
			//}

			if (time2 < alku)
				continue;
			else if (time2 > loppu) {
				//int_loc[1] = pa;
				//tpoints[ll] = static_cast<uint32_t>(jj);
				//if (begin) {
				//	int_loc[0] = -1;
				//	int_loc[1] = -1;
				//}
				break;
			}
			uint32_t ring_number1, ring_number2;

			// Detector numbers axially
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
			// Detector number transaxially
			const uint32_t ring_pos1 = (rsectorID1 % blocks_per_ring) * cryst_per_block + (crystalID1 % cryst_per_block);
			const uint32_t ring_pos2 = (rsectorID2 % blocks_per_ring) * cryst_per_block + (crystalID2 % cryst_per_block);
			// Detector numbers
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
							if (large_case)
								Ltrues[kk] = 1u;
							else
								Ltrues[L1][L2]++;
							if (source)
								trues_loc[kk] = 1u;
						}
						else if (event_scattered && store_scatter) {
							if (large_case)
								Lscatter[kk] = 1u;
							else
								Lscatter[L1][L2]++;
							if (source)
								scatter_loc[kk] = 1u;
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
						if (large_case)
							Lrandoms[kk] = 1u;
						else
							Lrandoms[L1][L2]++;
						if (source)
							randoms_loc[kk] = 1u;
					}
					else
						Lrandoms[kk] = 1u;
				}
				//if (!event_true && obtain_trues && source && outsize2 == 1ULL)
				//	trues_loc[kk] = 0u;
				//if (!event_scattered && store_scatter && source && outsize2 == 1ULL)
				//	scatter_loc[kk] = 0u;
				//if (event_true && store_randoms && source && outsize2 == 1ULL)
				//	randoms_loc[kk] = 0u;
			}
			if (outsize2 == 1) {
				if (large_case) {
					LL1[kk] = static_cast<uint16_t>(L1 + 1);
					LL2[kk] = static_cast<uint16_t>(L2 + 1);
				}
				else
					LL1[L1][L2]++;
			}
			else {
				LL1[kk] = static_cast<uint16_t>(L1 + 1);
				LL2[kk] = static_cast<uint16_t>(L2 + 1);
			}
			if (outsize2 > 1ULL && time2 >= aika) {
				tpoints[ll++] = static_cast<uint32_t>(jj);
				aika = time_intervals[++pa];
			}
			if (source) {
				S[kk][0] = sourcePosX1;
				S[kk][1] = sourcePosY1;
				S[kk][2] = sourcePosZ1;
				S[kk][3] = sourcePosX2;
				S[kk][4] = sourcePosY2;
				S[kk][5] = sourcePosZ2;
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
		//std::ostringstream stream5;
		//stream5 << "ROOT file load finished" << std::endl;
		//displayOnMATLAB(stream5);


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
			int_loc_delay[0] = 0;
			tpoints_delay[ll] = 0;

			for (uint64_t kk = 0; kk < Ndelays; kk++) {
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
					//pa--;
					begin = false;
					tpoints_delay[ll++] = 0u;
					aika = time_intervals[pa];
					int_loc_delay[0] = pa;
				}
				if (outsize2 == 1) {
					if (large_case) {
						Ldelay1[kk] = static_cast<uint16_t>(L1 + 1);
						Ldelay2[kk] = static_cast<uint16_t>(L2 + 1);
					}
					else
						Ldelay1[L1][L2]++;
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
	}

	void displayOnMATLAB(std::ostringstream& stream) {
		// Pass stream content to MATLAB fprintf function
		matlabPtr->feval(u"fprintf", 0,	std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
		// Clear stream buffer
		stream.str("");
	}

};
