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
#include <TROOT.h>
#include "TChain.h"
#include "saveSinogram.h"


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
		const matlab::data::TypedArray<uint8_t> scatter_components = std::move(inputs[14]);
		const bool randoms_correction = inputs[15][0];
		const bool store_coordinates = inputs[16][0];
		const uint32_t transaxial_multip = inputs[17][0];
		const uint32_t cryst_per_block_z = inputs[18][0];
		const uint32_t rings = inputs[19][0];
		const bool large_case = inputs[20][0];
		const bool TOF = inputs[21][0];
		const bool verbose = inputs[22][0];
		const int32_t nLayers = inputs[23][0];
		size_t outsize2 = (loppu - alku) / vali;

		bool dynamic = outsize2 > 1 ? true : false;


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
		matlab::data::TypedArray<uint8_t> layer1 = factory.createArray<uint8_t>({ 1, 1 });
		matlab::data::TypedArray<uint8_t> layer2 = factory.createArray<uint8_t>({ 1, 1 });
		matlab::data::TypedArray<uint8_t> layerD1 = factory.createArray<uint8_t>({ 1, 1 });
		matlab::data::TypedArray<uint8_t> layerD2 = factory.createArray<uint8_t>({ 1, 1 });

		if (outsize2 == 1ULL && !large_case) {
			if (obtain_trues) {
				Ltrues = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Ltrues.begin(), Ltrues.end(), static_cast<uint16_t>(0));
			}
			if (store_randoms) {
				Lrandoms = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Lrandoms.begin(), Lrandoms.end(), static_cast<uint16_t>(0));
			}
			if (store_scatter) {
				Lscatter = factory.createArray<uint16_t>({ detectors, detectors });
				std::fill(Lscatter.begin(), Lscatter.end(), static_cast<uint16_t>(0));
			}
			if (randoms_correction) {
				Ldelay1 = factory.createArray<uint16_t>({ detectors, detectors });
				Ldelay2 = factory.createArray<uint16_t>({ 1, 1 }, { static_cast<uint16_t>(0) });
				std::fill(Ldelay1.begin(), Ldelay1.end(), static_cast<uint16_t>(0));
			}
			LL1 = factory.createArray<uint16_t>({ detectors, detectors });
			std::fill(LL1.begin(), LL1.end(), static_cast<uint16_t>(0));
		}
		else {
			LL1 = factory.createArray<uint16_t>({ Nentries, 1 });
			LL2 = factory.createArray<uint16_t>({ Nentries, 1 });
			std::fill(LL1.begin(), LL1.end(), static_cast<uint16_t>(0));
			std::fill(LL2.begin(), LL2.end(), static_cast<uint16_t>(0));

			if (obtain_trues) {
				Ltrues = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Ltrues.begin(), Ltrues.end(), static_cast<uint16_t>(0));
			}
			if (store_randoms) {
				Lrandoms = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Lrandoms.begin(), Lrandoms.end(), static_cast<uint16_t>(0));
			}
			if (store_scatter) {
				Lscatter = factory.createArray<uint16_t>({ Nentries, 1 });
				std::fill(Lscatter.begin(), Lscatter.end(), static_cast<uint16_t>(0));
			}
			if (randoms_correction) {
				Ldelay1 = factory.createArray<uint16_t>({ Ndelays, 1 });
				Ldelay2 = factory.createArray<uint16_t>({ Ndelays, 1 });
				std::fill(Ldelay1.begin(), Ldelay1.end(), static_cast<uint16_t>(0));
				std::fill(Ldelay2.begin(), Ldelay2.end(), static_cast<uint16_t>(0));
			}
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
		if (nLayers > 1) {
			layer1 = factory.createArray<uint8_t>({ Nentries, 1 });
			layer2 = factory.createArray<uint8_t>({ Nentries, 1 });
			if (randoms_correction) {
				layerD1 = factory.createArray<uint8_t>({ Ndelays, 1 });
				layerD2 = factory.createArray<uint8_t>({ Ndelays, 1 });
			}
		}
		std::fill(layer1.begin(), layer1.end(), static_cast<uint8_t>(0));
		std::fill(layer2.begin(), layer2.end(), static_cast<uint8_t>(0));
		std::fill(layerD1.begin(), layerD1.end(), static_cast<uint8_t>(0));
		std::fill(layerD2.begin(), layerD2.end(), static_cast<uint8_t>(0));
		matlab::data::TypedArray<uint32_t> tpoints = factory.createArray<uint32_t>({ outsize2 + 2, 1 });
		std::fill(tpoints.begin(), tpoints.end(), 0U);
		matlab::data::TypedArray<uint32_t> tpoints_delay = factory.createArray<uint32_t>({ 1, 1 }, { 0 });
		matlab::data::TypedArray<float> S = factory.createArray<float>({ 1, 1 }, { 0.f });
		matlab::data::TypedArray<bool> trues_loc = factory.createArray<bool>({ 1, 1 }, { false });
		matlab::data::TypedArray<bool> randoms_loc = factory.createArray<bool>({ 1, 1 }, { false });
		matlab::data::TypedArray<bool> scatter_loc = factory.createArray<bool>({ 1, 1 }, { false });
		if (randoms_correction) {
			tpoints_delay = factory.createArray<uint32_t>({ outsize2 + 2, 1 });
			std::fill(tpoints_delay.begin(), tpoints_delay.end(), 0U);
		}

		if (source) {
			S = factory.createArray<float>({ Nentries, 6 });
			std::fill(S.begin(), S.end(), 0.f);
		}
		matlab::data::TypedArray<int32_t> int_loc = factory.createArray<int32_t>({ 2, 1 }, { 0, 0 });

		matlab::data::TypedArray<int32_t> int_loc_delay = factory.createArray<int32_t>({ 2, 1 }, { 0, 0 });

		if (obtain_trues && outsize2 == 1ULL) {
			trues_loc = factory.createArray<bool>({ Nentries, 1 });
			std::fill(trues_loc.begin(), trues_loc.end(), false);
		}
		if (store_randoms && outsize2 == 1ULL) {
			randoms_loc = factory.createArray<bool>({ Nentries, 1 });
			std::fill(randoms_loc.begin(), randoms_loc.end(), false);
		}
		if (store_scatter && outsize2 == 1ULL) {
			scatter_loc = factory.createArray<bool>({ Nentries, 1 });
			std::fill(scatter_loc.begin(), scatter_loc.end(), false);
		}
		matlab::data::TypedArray<double> TP = factory.createArray<double>({ 1, 1 });
		if ((dynamic && !TOF) || (!dynamic && TOF))
			TP = factory.createArray<double>({ Nentries, 1 });
		else if (dynamic && TOF)
			TP = factory.createArray<double>({ Nentries, 2 });
		std::fill(TP.begin(), TP.end(), 0.);
		matlab::data::TypedArray<double> TPd = factory.createArray<double>({ 1, 1 });
		if (dynamic && randoms_correction)
			TPd = factory.createArray<double>({ Ndelays, 1 });
		std::fill(TPd.begin(), TPd.end(), 0.);

		histogram(LL1, LL2, tpoints, vali, alku, loppu, outsize2, detectors, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S,
			Coincidences, Nentries, time_intervals, int_loc, obtain_trues, store_scatter, store_randoms, scatter_components, Ltrues, Lscatter,
			Lrandoms, trues_loc, Ndelays, randoms_correction, delay, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_loc, scatter_loc, x1, x2,
			y1, y2, z1, z2, store_coordinates, dynamic, large_case, rings, cryst_per_block_z, transaxial_multip, TP, TPd, TOF, verbose, layer1, layer2,
			layerD1, layerD2, nLayers);

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
		outputs[21] = std::move(TP);
		outputs[22] = std::move(TPd);
		outputs[23] = std::move(layer1);
		outputs[24] = std::move(layer2);
		outputs[25] = std::move(layerD1);
		outputs[26] = std::move(layerD2);

		//mexEvalString("pause(.001);");

		//std::ostringstream stream1;
		//stream1 << "Outputs formed" << std::endl;
		//displayOnMATLAB(stream1);
		//gROOT->Reset();

	}

	void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {
		if (inputs.size() != 24) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("24 inputs required") }));
		}

		if (outputs.size() > 27) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Too many output arguments") }));
		}

		if (inputs[0].getType() != matlab::data::ArrayType::CHAR) {
			matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("First input must be char") }));
		}
	}

	void histogram(matlab::data::TypedArray<uint16_t>& LL1, matlab::data::TypedArray<uint16_t>& LL2, matlab::data::TypedArray<uint32_t>& tpoints, double vali, const double alku, const double loppu, 
		const size_t outsize2, const uint32_t detectors, bool source, const uint32_t linear_multip, const uint32_t cryst_per_block, const uint32_t blocks_per_ring,	const uint32_t det_per_ring, 
		matlab::data::TypedArray<float>& S, TTree* Coincidences, const size_t Nentries, const matlab::data::TypedArray<double>& time_intervals, matlab::data::TypedArray<int32_t>& int_loc,
		bool obtain_trues, bool store_scatter, bool store_randoms, const matlab::data::TypedArray<uint8_t> scatter_components, matlab::data::TypedArray<uint16_t>& Ltrues, 
		matlab::data::TypedArray<uint16_t>& Lscatter, matlab::data::TypedArray<uint16_t>& Lrandoms, matlab::data::TypedArray<bool>& trues_loc, const size_t Ndelays, bool randoms_correction, 
		TTree* delay, matlab::data::TypedArray<uint16_t>& Ldelay1, matlab::data::TypedArray<uint16_t>& Ldelay2, matlab::data::TypedArray<int32_t>& int_loc_delay, 
		matlab::data::TypedArray<uint32_t>& tpoints_delay, matlab::data::TypedArray<bool>& randoms_loc, matlab::data::TypedArray<bool>& scatter_loc,
		matlab::data::TypedArray<float>& x1, matlab::data::TypedArray<float>& x2, matlab::data::TypedArray<float>& y1, matlab::data::TypedArray<float>& y2,
		matlab::data::TypedArray<float>& z1, matlab::data::TypedArray<float>& z2, bool store_coordinates, const bool dynamic, const bool large_case, const uint32_t rings, 
		const uint32_t cryst_per_block_z, const uint32_t transaxial_multip, matlab::data::TypedArray<double>& TP, matlab::data::TypedArray<double>& TPd, const bool TOF, 
		const bool verbose, matlab::data::TypedArray<uint8_t>& layer1, matlab::data::TypedArray<uint8_t>& layer2, matlab::data::TypedArray<uint8_t>& layerD1,
		matlab::data::TypedArray<uint8_t>& layerD2, const int32_t nLayers)
	{

		Int_t crystalID1 = 0, crystalID2 = 0, moduleID1 = 0, moduleID2 = 0, submoduleID1 = 0, submoduleID2 = 0, rsectorID1, rsectorID2, eventID1, eventID2, comptonPhantom1 = 0, comptonPhantom2 = 0,
			comptonCrystal1 = 0, comptonCrystal2 = 0, RayleighPhantom1 = 0, RayleighPhantom2 = 0, RayleighCrystal1 = 0, RayleighCrystal2 = 0, layerID1 = 0, layerID2 = 0;
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
		if (Coincidences->GetBranchStatus("crystalID2"))
			Coincidences->SetBranchAddress("crystalID2", &crystalID2);
		else {
			std::ostringstream stream;
			stream << "No crystal location information was found from file. Aborting." << std::endl;
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
		Coincidences->SetBranchAddress("layerID1", &layerID1);
		Coincidences->SetBranchAddress("layerID2", &layerID2);
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
		if (Coincidences->GetBranchStatus("time1"))
			Coincidences->SetBranchAddress("time1", &time1);
		else {
			if (TOF) {
				std::ostringstream stream;
				stream << "TOF examination selected, but no time information was found from file. Aborting." << std::endl;
				return;
			}
		}
		if (Coincidences->GetBranchStatus("time2"))
			Coincidences->SetBranchAddress("time2", &time2);
		else {
			if (dynamic || TOF) {
				std::ostringstream stream;
				stream << "Dynamic or TOF examination selected, but no time information was found from file. Aborting." << std::endl;
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
		if (obtain_trues || store_scatter || store_randoms) {
			//if (scatter_components[0]) {
			if (Coincidences->GetBranchStatus("comptonPhantom1"))
				Coincidences->SetBranchAddress("comptonPhantom1", &comptonPhantom1);
			else
				any++;
			if (Coincidences->GetBranchStatus("comptonPhantom2"))
				Coincidences->SetBranchAddress("comptonPhantom2", &comptonPhantom2);
			else
				any++;

			if (store_scatter && any == 2 && scatter_components[0] >= 1) {
				std::ostringstream stream;
				stream << "Compton phantom selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (store_scatter && scatter_components[0] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Compton scatter in the phantom will be stored" << std::endl;
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

			if (store_scatter && ((any == 4 && next == 1) || (any == 2 && next == 0)) && scatter_components[1] >= 1) {
				std::ostringstream stream;
				stream << "Compton crystal selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (store_scatter && scatter_components[1] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Compton scatter in the detector will be stored" << std::endl;
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

			if (store_scatter && ((any == 6 && next == 2) || (any == 2 && next == 0) || (any == 4 && next == 1)) && scatter_components[2] >= 1) {
				std::ostringstream stream;
				stream << "Rayleigh phantom selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (store_scatter && scatter_components[2] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Rayleigh scatter in the phantom will be stored" << std::endl;
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

			if (store_scatter && ((any == 8 && next == 3) || (any == 2 && next == 0) || (any == 4 && next == 1) || (any == 6 && next == 2)) && scatter_components[3] >= 1) {
				std::ostringstream stream;
				stream << "Rayleigh crystal selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (store_scatter && scatter_components[3] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Rayleigh scatter in the detector will be stored" << std::endl;
				displayOnMATLAB(stream);
			}


			if (store_scatter && any == 8) {
				std::ostringstream stream;
				stream << "Store scatter selected, but no scatter data was found from ROOT-file" << std::endl;
				displayOnMATLAB(stream);
				store_scatter = false;
			}


			if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] >= 1 && scatter_components[2] >= 1 && scatter_components[3] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the phantom and detector are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] >= 1 && scatter_components[2] >= 1 && scatter_components[3] == 0 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the phantom are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] >= 1 && scatter_components[2] == 0 && scatter_components[3] == 0 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and detector are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] >= 1 && scatter_components[2] == 0 && scatter_components[3] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and detector and Rayleigh scattered coincidences in the detector are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] == 0 && scatter_components[2] >= 1 && scatter_components[3] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the phantom and detector are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] == 0 && scatter_components[2] == 0 && scatter_components[3] >= 1 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the detector are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] == 0 && scatter_components[2] >= 1 && scatter_components[3] == 0 && verbose) {
				std::ostringstream stream;
				stream << "Randoms, Compton scattered coincidences in the phantom and Rayleigh scattered coincidences in the phantom are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
			}
			else if (obtain_trues && scatter_components[0] >= 1 && scatter_components[1] == 0 && scatter_components[2] == 0 && scatter_components[3] == 0 && verbose) {
				std::ostringstream stream;
				stream << "Randoms and Compton scattered coincidences in the phantom are NOT included in trues" << std::endl;
				displayOnMATLAB(stream);
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
				break;
			}
			if (dynamic && !TOF)
				TP[kk] = time2;
			else if (dynamic && TOF)
				TP[kk][0] = time2;
			uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
			detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1, 
				submoduleID2, rsectorID1, rsectorID2,crystalID1, crystalID2, cryst_per_block, cryst_per_block_z, transaxial_multip, rings);

			if (!dynamic && TOF) {
				double aika = (time2 - time1);
				if (ring_pos2 > ring_pos1)
					aika = -aika;
				TP[kk] = aika;
			}
			else if (dynamic && TOF) {
				double aika = (time2 - time1);
				if (ring_pos2 > ring_pos1)
					aika = -aika;
				TP[kk][1] = aika;
			}
			if (begin) {
				while (time2 >= time_intervals[pa])
					pa++;
				begin = false;
				tpoints[ll++] = 0u;
				aika = time_intervals[pa];
				int_loc[0] = pa;
			}
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
					else if ((comptonCrystal1 > 0 || comptonCrystal2 > 0) && scatter_components[1] > 0) {
						event_true = false;
						if ((scatter_components[1] <= comptonCrystal1 || scatter_components[1] <= comptonCrystal2))
							store_scatter_event = true;
					}
					else if ((RayleighPhantom1 > 0 || RayleighPhantom2 > 0) && scatter_components[2] > 0) {
						event_true = false;
						if ((scatter_components[2] <= RayleighPhantom1 || scatter_components[2] <= RayleighPhantom2))
							store_scatter_event = true;
					}
					else if ((RayleighCrystal1 > 0 || RayleighCrystal2 > 0) && scatter_components[3] > 0) {
						event_true = false;
						if ((scatter_components[3] <= RayleighCrystal1 || scatter_components[3] <= RayleighCrystal2))
							store_scatter_event = true;
					}
					else
						event_scattered = false;
				}
			}
			// Detector numbers
			uint32_t L1 = static_cast<uint32_t>(ring_number1) * det_per_ring + static_cast<uint32_t>(ring_pos1);
			uint32_t L2 = static_cast<uint32_t>(ring_number2) * det_per_ring + static_cast<uint32_t>(ring_pos2);
			if (L2 > L1 && !large_case) {
				const uint32_t L3 = L1;
				L1 = L2;
				L2 = L3;
			}
			if (nLayers > 1) {
				layer1[kk] = layerID2;
				layer2[kk] = layerID1;
			}
			if (obtain_trues || store_scatter || store_randoms) {
				if ((event_true && obtain_trues) || (store_scatter_event && store_scatter)) {
					if (outsize2 == 1ULL) {
						if (event_true && obtain_trues) {
							if (large_case)
								Ltrues[kk] = 1u;
							else
								Ltrues[L1][L2]++;
						}
						else if (store_scatter_event && store_scatter) {
							if (large_case)
								Lscatter[kk] = 1u;
							else
								Lscatter[L1][L2]++;
						}
					}
					else {
						if (event_true && obtain_trues)
							Ltrues[kk] = 1u;
						else if (store_scatter_event && store_scatter)
							Lscatter[kk] = 1u;
					}
				}
				else if (!event_true && store_randoms && !event_scattered) {
					if (outsize2 == 1ULL) {
						if (large_case)
							Lrandoms[kk] = 1u;
						else
							Lrandoms[L1][L2]++;
					}
					else
						Lrandoms[kk] = 1u;
				}
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
			if (outsize2 == 1ULL) {
				if (event_true && obtain_trues)
					trues_loc[kk] = true;
				else if (store_scatter_event && store_scatter)
					scatter_loc[kk] = true;
				else if (!event_true && store_randoms && !event_scattered)
					randoms_loc[kk] = true;
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
			delay->SetBranchAddress("submoduleID1", &submoduleID1);
			delay->SetBranchAddress("submoduleID2", &submoduleID2);
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

				if (dynamic)
					TPd[kk] = time2;
				uint32_t ring_number1 = 0, ring_number2 = 0, ring_pos1 = 0, ring_pos2 = 0;
				detectorIndices(ring_number1, ring_number2, ring_pos1, ring_pos2, blocks_per_ring, linear_multip, no_modules, no_submodules, moduleID1, moduleID2, submoduleID1,
					submoduleID2, rsectorID1, rsectorID2, crystalID1, crystalID2, cryst_per_block, cryst_per_block_z, transaxial_multip, rings);
				uint32_t L1 = ring_number1 * det_per_ring + ring_pos1;
				uint32_t L2 = ring_number2 * det_per_ring + ring_pos2;
				if (begin) {
					while (time2 >= time_intervals[pa])
						pa++;
					//pa--;
					begin = false;
					tpoints_delay[ll++] = 0u;
					aika = time_intervals[pa];
					int_loc_delay[0] = pa;
				}
				if (nLayers > 1) {
					layerD1[kk] = layerID2;
					layerD2[kk] = layerID1;
				}
				if (outsize2 == 1) {
					if (large_case) {
						Ldelay1[kk] = static_cast<uint16_t>(L1 + 1);
						Ldelay2[kk] = static_cast<uint16_t>(L2 + 1);
					}
					else {
						if (L2 > L1) {
							const uint32_t L3 = L1;
							L1 = L2;
							L2 = L3;
						}
						Ldelay1[L1][L2]++;
					}
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
		matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
		// Clear stream buffer
		stream.str("");
	}

};
