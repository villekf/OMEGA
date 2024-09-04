// Sinogram creation for, e.g., Python
#if defined(_MSC_VER)
#define DLL_FUNCTION __declspec(dllexport)
#endif
#include "sinoLib.h"

int sinoMain(const uint16_t* ringPos1, const uint16_t* ringPos2, const uint16_t* ringNumber1, const uint16_t* ringNumber2, const bool* truesIndex, const bool* scatterIndex,
	const bool* randomsIndex, const uint64_t sinoSize, const uint32_t Ndist, const uint32_t Nang, const uint32_t ringDifference, const uint32_t span, const uint32_t* seg, const uint64_t pituus, 
	const uint64_t TOFSize, const uint16_t* time, const uint64_t NT, const int32_t detPerRing, const int32_t rings, const uint16_t* bins, const int32_t nDistSide, const int32_t detWPseudo,
	const int32_t nPseudos, const int32_t crystPerBlock, const int32_t nLayers, const uint8_t* layer1, const uint8_t* layer2, const int64_t koko, const int64_t tSize, const int64_t sSize, const int64_t rSize, 
	uint16_t* Sino, uint16_t* SinoT, uint16_t* SinoC, uint16_t* SinoR) {

	bool storeTrues = false;
	bool storeScatter = false;
	bool storeRandoms = false;
	if (tSize > 1) {
		storeTrues = true;
	}
	if (sSize > 1) {
		storeScatter = true;
	}
	if (rSize > 1) {
		storeRandoms = true;
	}

	openMPSino(ringPos1, ringPos2, ringNumber1, ringNumber2, truesIndex, scatterIndex, randomsIndex, sinoSize, Ndist, Nang, ringDifference,
		span, seg, time, NT, TOFSize, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko,
		bins, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock, layer1, layer2, nLayers);
	return 0;
}