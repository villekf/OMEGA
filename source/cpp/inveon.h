#pragma once

#include <cstdint>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cstdio>
#include <thread>
#include <atomic>
#ifdef MATLAB
#include "mex.h"
#endif

#define DET_PER_RING 320
#define RINGS 80

void saveSinogram(uint16_t L1, uint16_t L2, uint16_t* Sino, const uint32_t Ndist, const uint32_t Nang, const uint32_t ring_difference, const uint32_t span,
	const uint64_t sinoSize, const uint32_t* seg, const int32_t nDistSide, const int tPoint = 0) {
	int32_t ring_pos1 = L1 % DET_PER_RING;
	int32_t ring_pos2 = L2 % DET_PER_RING;
	int32_t ring_number1 = L1 / DET_PER_RING;
	int32_t ring_number2 = L2 / DET_PER_RING;
	const int32_t xa = std::max(ring_pos1, ring_pos2);
	const int32_t ya = std::min(ring_pos1, ring_pos2);
	int32_t j = ((xa + ya + DET_PER_RING / 2 + 1) % DET_PER_RING) / 2;
	const int32_t b = j + DET_PER_RING / 2;
	int32_t i = std::abs(xa - ya - DET_PER_RING / 2);
	const bool ind = ya < j || b < xa;
	if (ind)
		i = -i;
	const bool swap = (j * 2) < -i || i <= ((j - DET_PER_RING / 2) * 2);
	bool accepted_lors;
	if (Ndist % 2U == 0)
		accepted_lors = (i <= (static_cast<int32_t>(Ndist) / 2 + std::min(0, nDistSide)) && i >= (-static_cast<int32_t>(Ndist) / 2 + std::max(0, nDistSide)));
	else
		accepted_lors = (i <= static_cast<int32_t>(Ndist) / 2 && i >= (-static_cast<int32_t>(Ndist) / 2));
	accepted_lors = accepted_lors && (std::abs(ring_number1 - ring_number2) <= ring_difference);
	int32_t sinoIndex = 0;
	if (accepted_lors) {
		j = j / (DET_PER_RING / 2 / Nang);
		if (swap) {
			const int32_t ring_number3 = ring_number1;
			ring_number1 = ring_number2;
			ring_number2 = ring_number3;
		}
		i = i + Ndist / 2 - std::max(0, nDistSide);
		const bool swappi = ring_pos2 > ring_pos1;
		if (swappi) {
			const int32_t ring_number3 = ring_number1;
			ring_number1 = ring_number2;
			ring_number2 = ring_number3;
		}
		if (span <= 1) {
			sinoIndex = ring_number2 * RINGS + ring_number1;
		}
		else {
			const int32_t erotus = ring_number1 - ring_number2;
			const int32_t summa = ring_number1 + ring_number2;
			if (std::abs(erotus) <= span / 2) {
				sinoIndex = summa;
			}
			else {
				sinoIndex = ((std::abs(erotus) + (span / 2)) / span);
				if (erotus < 0) {
					sinoIndex = (summa - ((span / 2) * (sinoIndex * 2 - 1) + sinoIndex)) + static_cast<uint32_t>(seg[(sinoIndex - 1) * 2]);
				}
				else {
					sinoIndex = (summa - ((span / 2) * (sinoIndex * 2 - 1) + sinoIndex)) + static_cast<uint32_t>(seg[(sinoIndex - 1) * 2 + 1]);
				}
			}
		}
		const uint64_t indeksi = static_cast<uint64_t>(i) + static_cast<uint64_t>(j) * static_cast<uint64_t>(Ndist) +
			static_cast<uint64_t>(sinoIndex) * static_cast<uint64_t>(Ndist) * static_cast<uint64_t>(Nang) + sinoSize * static_cast<uint64_t>(tPoint);
		Sino[indeksi]++;
	}
	else
		return;
}

// Parameters that stay constant during the data loading
struct InveonParams {
	const double* vali;
	double alku;
	double loppu;
	uint32_t detectors;
	bool randoms_correction;
	uint32_t Ndist;
	uint32_t Nang;
	uint32_t ringDifference;
	uint32_t span;
	uint64_t sinoSize;
	const uint32_t* seg;
	int32_t nDistSide;
	bool storeCoordinates;
	uint64_t Nt;
};

// Output arrays. In the counting mode (LL1 == nullptr) nothing is written, only the counters are increased.
struct InveonOutputs {
	uint16_t* LL1;
	uint16_t* LL2;
	uint16_t* tpoints;
	uint16_t* DD1;
	uint16_t* DD2;
	uint16_t* Sino;
	uint16_t* SinoD;
};

// Decoding state shared between the packet blocks (and, in the parallel version, between the chunks)
struct InveonState {
	double ms = 0.;		// seconds
	double aika = 0.;
	uint64_t tPoint = 0;
	uint64_t ll = 0;
	uint64_t dd = 0;
	uint64_t tBase = 0;	// the first time point stored in the sinogram buffer (parallel version)
};

// Processes nPackets 6-byte packets starting at buf (buf must be readable for 2 bytes beyond the last packet).
// Returns true if the loop has to stop (Nt time points reached).
static inline bool processInveonBlock(const uint8_t* buf, const size_t nPackets, InveonState& st, const InveonParams& P, const InveonOutputs& O)
{
	const double* vali = P.vali;
	const double alku = P.alku;
	const double loppu = P.loppu;
	const uint32_t detectors = P.detectors;
	const bool randoms_correction = P.randoms_correction;
	const bool storeCoordinates = P.storeCoordinates;
	const uint64_t Nt = P.Nt;
	for (size_t k = 0; k < nPackets; k++) {
		// ms is checked before decoding each packet (using the value after the previous packet)
		if (st.ms > loppu)
			return false;
		uint64_t ew1;
		std::memcpy(&ew1, buf + k * 6, 8);
		ew1 &= 0xFFFFFFFFFFFFULL; // 48-bit packet, the two top bytes are zero (little-endian)

		const int tag = static_cast<int>((ew1 >> 43) & 1);

		if (!tag) {
			if (st.ms >= alku) {
				const int prompt = static_cast<int>((ew1 >> 42) & 1);
				if (prompt) {
					uint32_t L1 = (ew1 >> 19) & 0x1ffff;
					uint32_t L2 = ew1 & 0x1ffff;
					if (L1 >= detectors || L2 >= detectors)
						continue;
					if (!storeCoordinates)
						saveSinogram(L1, L2, O.Sino, P.Ndist, P.Nang, P.ringDifference, P.span, P.sinoSize, P.seg, P.nDistSide, static_cast<int>(st.tPoint - st.tBase));
					else {
						if (O.LL1) {
							if (L2 > L1) {
								const uint32_t L3 = L1;
								L1 = L2;
								L2 = L3;
							}
							O.LL1[st.ll] = static_cast<uint16_t>(L1 + 1);
							O.LL2[st.ll] = static_cast<uint16_t>(L2 + 1);
							if (Nt > 1)
								O.tpoints[st.ll] = static_cast<uint16_t>(st.tPoint + 1);
						}
						st.ll++;
					}
				}
				else if (randoms_correction) {
					uint32_t L1 = (ew1 >> 19) & 0x1ffff;
					uint32_t L2 = ew1 & 0x1ffff;
					if (L1 >= detectors || L2 >= detectors)
						continue;
					if (!storeCoordinates)
						saveSinogram(L1, L2, O.SinoD, P.Ndist, P.Nang, P.ringDifference, P.span, P.sinoSize, P.seg, P.nDistSide, static_cast<int>(st.tPoint - st.tBase));
					else {
						if (O.DD1) {
							if (L2 > L1) {
								const uint32_t L3 = L1;
								L1 = L2;
								L2 = L3;
							}
							O.DD1[st.dd] = static_cast<uint16_t>(L1 + 1);
							O.DD2[st.dd] = static_cast<uint16_t>(L2 + 1);
						}
						st.dd++;
					}
				}
			}
		}
		else {
			if (((ew1 >> 36) & 0xff) == 160) { // Elapsed Time Tag Packet
				st.ms += 200e-6; // 200 microsecond increments

				if (Nt > 1 && st.tPoint < Nt && st.ms >= st.aika + vali[st.tPoint]) {
					st.aika += vali[st.tPoint];
					st.tPoint++;
					if (st.tPoint == Nt)
						return true;
				}
			}
		}
	}
	return false;
}

static inline FILE* inveonOpen(const char* name) {
	FILE* f = nullptr;
#if (defined(WIN32) || defined(_WIN32) || (defined(__WIN32) && !defined(__CYGWIN__)) || defined(_WIN64)) && defined(_MSC_VER)
	if (fopen_s(&f, name, "rb") != 0)
		f = nullptr;
#else
	f = fopen(name, "rb");
#endif
	return f;
}

static inline int inveonSeek(FILE* f, const uint64_t offset) {
#if defined(_WIN32) || defined(_WIN64) || defined(WIN32)
	return _fseeki64(f, static_cast<long long>(offset), SEEK_SET);
#else
	return fseeko(f, static_cast<off_t>(offset), SEEK_SET);
#endif
}

static inline int64_t inveonFileSize(FILE* f) {
#if defined(_WIN32) || defined(_WIN64) || defined(WIN32)
	if (_fseeki64(f, 0, SEEK_END) != 0)
		return -1;
	const int64_t sz = static_cast<int64_t>(_ftelli64(f));
	_fseeki64(f, 0, SEEK_SET);
#else
	if (fseeko(f, 0, SEEK_END) != 0)
		return -1;
	const int64_t sz = static_cast<int64_t>(ftello(f));
	fseeko(f, 0, SEEK_SET);
#endif
	return sz;
}

// Runs fun(0), ..., fun(nThreads - 1), each in its own thread (the first one in the calling thread).
// If a thread cannot be created the work of that thread is done by the calling thread.
template <typename F>
static void inveonRunThreads(const uint32_t nThreads, F&& fun) {
	std::vector<std::thread> threads;
	std::vector<uint32_t> notStarted;
	threads.reserve(nThreads);
	for (uint32_t t = 1; t < nThreads; t++) {
		try {
			threads.emplace_back([&fun, t]() { fun(t); });
		}
		catch (...) {
			notStarted.push_back(t);
		}
	}
	fun(0);
	for (auto t : notStarted)
		fun(t);
	for (auto& th : threads)
		th.join();
}

// Bit-identical parallel version of the loading loop. Pass 1 counts the elapsed time tags of fixed-size blocks, the block start times (and
// time points) are then evaluated sequentially with exactly the same repeated additions as in the serial loop, after which the blocks are
// processed in parallel. With stored coordinates an additional counting pass is needed for the output offsets.
// Returns false if the parallel version could not be used (the outputs are unchanged then, assuming zero-initialized sinograms).
static bool histogramParallel(const char* fname, const uint64_t nPacketsTotal, const InveonParams& P, const InveonOutputs& O, uint32_t nThreads,
	double& msEnd, uint64_t& nPrompts, uint64_t& nDelays)
{
	const uint64_t BP = static_cast<uint64_t>(1) << 20; // packets per block
	const uint64_t nBlocks = (nPacketsTotal + BP - 1) / BP;
	if (nThreads < 2 || nBlocks < 2 * static_cast<uint64_t>(nThreads))
		return false;
	std::atomic<bool> failed(false);
	std::atomic<uint64_t> nextBlock(0);
	std::vector<uint64_t> nTags(nBlocks, 0);

	auto readBlock = [&](FILE* f, const uint64_t b, std::vector<uint8_t>& buf) -> uint64_t {
		const uint64_t cnt = std::min(BP, nPacketsTotal - b * BP);
		if (inveonSeek(f, b * BP * 6) != 0 || fread(buf.data(), 6, static_cast<size_t>(cnt), f) != static_cast<size_t>(cnt)) {
			failed = true;
			return 0;
		}
		return cnt;
	};

	// Pass 1: count the elapsed time tags
	inveonRunThreads(nThreads, [&](uint32_t) {
		FILE* f = nullptr;
		try {
			f = inveonOpen(fname);
			if (f == nullptr) {
				failed = true;
				return;
			}
			std::vector<uint8_t> buf(static_cast<size_t>(BP) * 6 + 8, 0);
			while (!failed) {
				const uint64_t b = nextBlock++;
				if (b >= nBlocks)
					break;
				const uint64_t cnt = readBlock(f, b, buf);
				if (cnt == 0)
					break;
				uint64_t n = 0;
				for (uint64_t k = 0; k < cnt; k++) {
					uint64_t ew1;
					std::memcpy(&ew1, buf.data() + k * 6, 8);
					n += static_cast<uint64_t>((((ew1 >> 36) & 0xff) == 160) && (((ew1 >> 43) & 1) == 1));
				}
				nTags[b] = n;
			}
		}
		catch (...) {
			failed = true;
		}
		if (f)
			fclose(f);
	});
	if (failed)
		return false;

	// Sequential evaluation of the block start/end states (identical arithmetic to the serial loop)
	struct BlockInfo {
		InveonState start;
		InveonState end;
		bool relevant = false;
		uint64_t nP = 0;
		uint64_t nD = 0;
		uint64_t offP = 0;
		uint64_t offD = 0;
	};
	std::vector<BlockInfo> blk(nBlocks);
	InveonState cur;
	cur.aika = P.alku;
	int64_t lastBlock = static_cast<int64_t>(nBlocks) - 1;
	for (uint64_t b = 0; b < nBlocks; b++) {
		blk[b].start = cur;
		if (cur.ms > P.loppu) {
			lastBlock = static_cast<int64_t>(b) - 1;
			break;
		}
		bool stopped = false;
		for (uint64_t t = 0; t < nTags[b]; t++) {
			cur.ms += 200e-6;
			if (P.Nt > 1 && cur.tPoint < P.Nt && cur.ms >= cur.aika + P.vali[cur.tPoint]) {
				cur.aika += P.vali[cur.tPoint];
				cur.tPoint++;
				if (cur.tPoint == P.Nt) {
					stopped = true;
					break;
				}
			}
			if (cur.ms > P.loppu) {
				stopped = true;
				break;
			}
		}
		blk[b].end = cur;
		blk[b].relevant = cur.ms >= P.alku; // otherwise no events are stored from this block
		if (stopped) {
			lastBlock = static_cast<int64_t>(b);
			break;
		}
	}
	msEnd = cur.ms;
	std::vector<uint64_t> rel;
	for (int64_t b = 0; b <= lastBlock; b++)
		if (blk[b].relevant)
			rel.push_back(static_cast<uint64_t>(b));
	if (rel.empty()) {
		nPrompts = 0;
		nDelays = 0;
		return true;
	}
	if (rel.size() < 2 * static_cast<size_t>(nThreads))
		nThreads = static_cast<uint32_t>(rel.size() / 2);
	if (nThreads < 2)
		return false;

	// Contiguous groups of relevant blocks, one per thread. Each group (except the first one, which writes directly to
	// the output) has its own sinogram copy that covers only the time points the group touches.
	struct Group {
		size_t first;
		size_t last;	// inclusive
		uint64_t fLo;
		uint64_t fHi;
	};
	std::vector<Group> groups;
	const uint64_t maxBytes = static_cast<uint64_t>(4) * 1024 * 1024 * 1024; // memory limit for the sinogram copies
	const uint64_t nSino = P.randoms_correction ? 2 : 1;
	for (;;) {
		groups.clear();
		uint64_t extra = 0;
		for (uint32_t g = 0; g < nThreads; g++) {
			Group gr;
			gr.first = static_cast<size_t>(static_cast<uint64_t>(rel.size()) * g / nThreads);
			gr.last = static_cast<size_t>(static_cast<uint64_t>(rel.size()) * (g + 1) / nThreads) - 1;
			gr.fLo = blk[rel[gr.first]].start.tPoint;
			gr.fHi = std::min(P.Nt - 1, blk[rel[gr.last]].end.tPoint);
			if (g > 0 && !P.storeCoordinates)
				extra += (gr.fHi - gr.fLo + 1) * P.sinoSize * 2 * nSino;
			groups.push_back(gr);
		}
		if (extra <= maxBytes)
			break;
		if (nThreads <= 2)
			return false;
		nThreads--;
	}

	// Counting pass for the stored coordinates (output offsets)
	if (P.storeCoordinates) {
		inveonRunThreads(nThreads, [&](uint32_t g) {
			FILE* f = nullptr;
			try {
				f = inveonOpen(fname);
				if (f == nullptr) {
					failed = true;
					return;
				}
				std::vector<uint8_t> buf(static_cast<size_t>(BP) * 6 + 8, 0);
				InveonOutputs Oc = {};
				for (size_t r = groups[g].first; r <= groups[g].last && !failed; r++) {
					const uint64_t b = rel[r];
					const uint64_t cnt = readBlock(f, b, buf);
					if (cnt == 0)
						break;
					InveonState st = blk[b].start;
					st.ll = 0;
					st.dd = 0;
					processInveonBlock(buf.data(), static_cast<size_t>(cnt), st, P, Oc);
					blk[b].nP = st.ll;
					blk[b].nD = st.dd;
				}
			}
			catch (...) {
				failed = true;
			}
			if (f)
				fclose(f);
		});
		if (failed)
			return false;
		uint64_t offP = 0, offD = 0;
		for (size_t r = 0; r < rel.size(); r++) {
			blk[rel[r]].offP = offP;
			blk[rel[r]].offD = offD;
			offP += blk[rel[r]].nP;
			offD += blk[rel[r]].nD;
		}
		nPrompts = offP;
		nDelays = offD;
	}
	else {
		nPrompts = 0;
		nDelays = 0;
	}

	// Main pass
	std::vector<std::vector<uint16_t>> sinoBuf(nThreads), sinoDBuf(nThreads);
	inveonRunThreads(nThreads, [&](uint32_t g) {
		FILE* f = nullptr;
		try {
			f = inveonOpen(fname);
			if (f == nullptr) {
				failed = true;
				return;
			}
			std::vector<uint8_t> buf(static_cast<size_t>(BP) * 6 + 8, 0);
			InveonOutputs Oc = O;
			uint64_t tBase = 0;
			if (!P.storeCoordinates && g > 0) {
				const size_t nEl = static_cast<size_t>((groups[g].fHi - groups[g].fLo + 1) * P.sinoSize);
				sinoBuf[g].assign(nEl, 0);
				Oc.Sino = sinoBuf[g].data();
				if (P.randoms_correction) {
					sinoDBuf[g].assign(nEl, 0);
					Oc.SinoD = sinoDBuf[g].data();
				}
				tBase = groups[g].fLo;
			}
			for (size_t r = groups[g].first; r <= groups[g].last && !failed; r++) {
				const uint64_t b = rel[r];
				const uint64_t cnt = readBlock(f, b, buf);
				if (cnt == 0)
					break;
				InveonState st = blk[b].start;
				st.tBase = tBase;
				st.ll = blk[b].offP;
				st.dd = blk[b].offD;
				processInveonBlock(buf.data(), static_cast<size_t>(cnt), st, P, Oc);
			}
		}
		catch (...) {
			failed = true;
		}
		if (f)
			fclose(f);
	});
	if (!failed && !P.storeCoordinates) {
		// Sum the sinogram copies to the output (uint16 wrap-around gives the same result as the serial increments)
		const uint64_t total = P.sinoSize * P.Nt;
		inveonRunThreads(nThreads, [&](uint32_t t) {
			const uint64_t lo = total * t / nThreads;
			const uint64_t hi = total * (t + 1) / nThreads;
			for (uint32_t g = 1; g < nThreads; g++) {
				const uint64_t gLo = groups[g].fLo * P.sinoSize;
				const uint64_t gHi = gLo + sinoBuf[g].size();
				const uint64_t a = std::max(lo, gLo);
				const uint64_t e = std::min(hi, gHi);
				for (uint64_t i = a; i < e; i++) {
					O.Sino[i] = static_cast<uint16_t>(O.Sino[i] + sinoBuf[g][i - gLo]);
					if (P.randoms_correction)
						O.SinoD[i] = static_cast<uint16_t>(O.SinoD[i] + sinoDBuf[g][i - gLo]);
				}
			}
		});
	}
	if (failed) {
		// Restore the zero-initialized outputs so that the serial version can be used
		if (!P.storeCoordinates) {
			std::fill(O.Sino, O.Sino + P.sinoSize * P.Nt, static_cast<uint16_t>(0));
			if (P.randoms_correction)
				std::fill(O.SinoD, O.SinoD + P.sinoSize * P.Nt, static_cast<uint16_t>(0));
		}
		return false;
	}
	return true;
}

// Returns 1 on success, 0 if the file could not be opened. The numbers of stored prompts and delays (only when storeCoordinates is true) are returned in nPrompts and nDelays.
// maxThreads: maximum number of threads (1 = serial). The result is bit-identical regardless of the number of threads.
int histogram(uint16_t* LL1, uint16_t* LL2, uint16_t* tpoints, const char* argv, const double* vali, const double alku, const double loppu,
	const uint32_t detectors, const size_t pituus, const bool randoms_correction, uint16_t* DD1, uint16_t* DD2, uint16_t* Sino, uint16_t* SinoD, const bool saveRawData,
	const uint32_t Ndist, const uint32_t Nang, const uint32_t ringDifference, const uint32_t span, const uint64_t sinoSize, const uint32_t* seg, const int32_t nDistSide,
	const bool storeCoordinates, const uint64_t Nt, uint64_t* nPrompts = nullptr, uint64_t* nDelays = nullptr, const uint32_t maxThreads = 16)
{
	if (nPrompts)
		*nPrompts = 0;
	if (nDelays)
		*nDelays = 0;

	FILE* streami = inveonOpen(argv);
	if (streami == NULL) {
#ifdef MATLAB
		mexErrMsgIdAndTxt("MATLAB:inveon_list2matlab:invalidFile",
			"Error opening file or no file opened");
#else
		fprintf(stdout, "Error opening file or no file opened\n");
		fflush(stdout);
#endif
		return 0;
	}

#ifdef MATLAB
	mexPrintf("File opened \n");
#else
	fprintf(stdout, "File opened \n");
	fflush(stdout);
#endif

	InveonParams P;
	P.vali = vali;
	P.alku = alku;
	P.loppu = loppu;
	P.detectors = detectors;
	P.randoms_correction = randoms_correction;
	P.Ndist = Ndist;
	P.Nang = Nang;
	P.ringDifference = ringDifference;
	P.span = span;
	P.sinoSize = sinoSize;
	P.seg = seg;
	P.nDistSide = nDistSide;
	P.storeCoordinates = storeCoordinates;
	P.Nt = Nt;
	InveonOutputs O;
	O.LL1 = LL1;
	O.LL2 = LL2;
	O.tpoints = tpoints;
	O.DD1 = DD1;
	O.DD2 = DD2;
	O.Sino = Sino;
	O.SinoD = SinoD;

	InveonState st;
	st.aika = alku;

	// Use the multithreaded version for large files
	uint32_t nThreads = std::min(maxThreads, static_cast<uint32_t>(16));
	const uint32_t hw = std::thread::hardware_concurrency();
	nThreads = std::min(nThreads, hw == 0 ? static_cast<uint32_t>(4) : hw);
	const int64_t fileSize = inveonFileSize(streami);
	bool done = false;
	if (nThreads > 1 && fileSize > 0) {
		double msEnd = 0.;
		uint64_t nP = 0, nD = 0;
		if (histogramParallel(argv, static_cast<uint64_t>(fileSize) / 6, P, O, nThreads, msEnd, nP, nD)) {
			st.ms = msEnd;
			st.ll = nP;
			st.dd = nD;
			done = true;
		}
	}

	// Serial version: read the file in large blocks instead of one packet at a time
	if (!done) {
		const size_t blockPackets = static_cast<size_t>(1) << 22;
		std::vector<uint8_t> buffer(blockPackets * 6 + 8, 0);
		bool stop = false;
		while (!stop) {
			const size_t nRead = fread(buffer.data(), 6, blockPackets, streami); // complete packets only
			if (nRead == 0)
				break;
			if (st.ms > loppu)
				break;
			stop = processInveonBlock(buffer.data(), nRead, st, P, O);
			if (st.ms > loppu)
				break;
			if (nRead < blockPackets)
				break;
		}
	}
	if (nPrompts)
		*nPrompts = st.ll;
	if (nDelays)
		*nDelays = st.dd;
#ifdef MATLAB
	mexPrintf("End time %f\n", st.ms);
	mexEvalString("pause(.0001);");
#else
	printf("End time %f\n", st.ms);
#endif
	fclose(streami);
	return 1;
}
