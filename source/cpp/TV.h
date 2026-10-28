#pragma once
#include "projector_functions.h"

template <typename T>
T sqrtVal(const float3<T> input, const T epps, const int type, const float3<T> w = float3<T>{}) {
if (type == 6)
	return std::sqrt(w.x * input.x * input.x + w.y * input.y * input.y + w.z * input.z * input.z + epps);
else
	return std::sqrt(input.x * input.x + input.y * input.y + input.z * input.z + epps);
}

// Image/reference image read, each coordinate is clamped to the volume (clamp-to-edge, i.e. Neumann boundary), same as the GPU kernels
template <typename T>
inline T tvRead(const T* buf, const int64_t x, const int64_t y, const int64_t z, const int64_t Nx, const int64_t Ny, const int64_t Nz) {
	const int64_t xc = std::max<int64_t>(0, std::min<int64_t>(x, Nx - 1));
	const int64_t yc = std::max<int64_t>(0, std::min<int64_t>(y, Ny - 1));
	const int64_t zc = std::max<int64_t>(0, std::min<int64_t>(z, Nz - 1));
	return buf[xc + yc * Nx + zc * Nx * Ny];
}

// APLS vector q = (g - (xi.g) xi) / sqrt(|g|^2 - (xi.g)^2 + epps), xi = r / sqrt(|r|^2 + eta^2)
// g = forward differences of the image, r = forward differences of the reference image
template <typename T>
inline float3<T> aplsQ(const float3<T> g, const float3<T> r, const T eta, const T epps) {
	const T rNorm = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z + eta * eta);
	const T xix = r.x / rNorm;
	const T xiy = r.y / rNorm;
	const T xiz = r.z / rNorm;
	const T gxi = g.x * xix + g.y * xiy + g.z * xiz;
	const T p = std::sqrt(std::max<T>(g.x * g.x + g.y * g.y + g.z * g.z - gxi * gxi, static_cast<T>(0.)) + epps);
	float3<T> q;
	MFLOAT3((g.x - gxi * xix) / p, (g.y - gxi * xiy) / p, (g.z - gxi * xiz) / p, q);
	return q;
}

// Anatomical TV type 1: q = 0.5 * grad_g Q(g) / sqrt(Q(g) + epps), Q(g) = g^T A g with the 3x3 weight matrix A of the voxel
// The 9 weights S[n + NN * kk] are read at the (clamped) voxel (x, y, z), g = forward differences of the image at that voxel
template <typename T>
inline float3<T> anatTVQ(const T* S, const int64_t x, const int64_t y, const int64_t z, const int64_t Nx, const int64_t Ny, const int64_t Nz, const float3<T> g, const T epps) {
	const int64_t xc = std::max<int64_t>(0, std::min<int64_t>(x, Nx - 1));
	const int64_t yc = std::max<int64_t>(0, std::min<int64_t>(y, Ny - 1));
	const int64_t zc = std::max<int64_t>(0, std::min<int64_t>(z, Nz - 1));
	const int64_t NN = Nx * Ny * Nz;
	const int64_t nc = xc + yc * Nx + zc * Nx * Ny;
	T s[9];
	for (int kk = 0; kk < 9; kk++)
		s[kk] = S[nc + NN * kk];
	const T sxy = s[1] + s[3];
	const T sxz = s[2] + s[6];
	const T syz = s[5] + s[7];
	const T p = std::sqrt(s[0] * g.x * g.x + s[4] * g.y * g.y + s[8] * g.z * g.z + sxy * g.x * g.y + sxz * g.x * g.z + syz * g.y * g.z + epps);
	float3<T> q;
	MFLOAT3((T).5 * ((T)2. * s[0] * g.x + sxy * g.y + sxz * g.z) / p, (T).5 * ((T)2. * s[4] * g.y + sxy * g.x + syz * g.z) / p, (T).5 * ((T)2. * s[8] * g.z + sxz * g.x + syz * g.y) / p, q);
	return q;
}

template <typename T>
void TVKernel(T* grad, const T* u, const int Nx, const int Ny, const int Nz, const int NxOrig, const int NyOrig, const int NzOrig, const T sigma, const T epps, const T beta, 
	const int type, const bool anatomical, const T C = 0.f, const T* S = nullptr, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr) {

	setThreads();
	int64_t start = 0;
	int64_t end = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny) * static_cast<int64_t>(Nz);
	int64_t Nxy = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny);
	const int64_t Nx64 = static_cast<int64_t>(Nx);
	const int64_t Ny64 = static_cast<int64_t>(Ny);
	const int64_t Nz64 = static_cast<int64_t>(Nz);
#ifdef _OPENMP
#if _OPENMP >= 201511
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t n = start; n < end; n++) {
		const int64_t z = n / Nxy;
		const int64_t y = (n - z * Nxy) / static_cast<int64_t>(Nx);
		const int64_t x = n - z * Nxy - y * static_cast<int64_t>(Nx);
		const T uijk = u[n];
		if (type == 4) { // START JPTV || SATV
			float2<T> ux;
			float2<T> uy;
			float2<T> uz;
			// Volume edge is clamped (Neumann boundary), identical to the GPU kernels
			ux.x = tvRead(u, x + 1, y, z, Nx64, Ny64, Nz64);
			ux.y = tvRead(u, x - 1, y, z, Nx64, Ny64, Nz64);
			uy.x = tvRead(u, x, y + 1, z, Nx64, Ny64, Nz64);
			uy.y = tvRead(u, x, y - 1, z, Nx64, Ny64, Nz64);
			uz.x = tvRead(u, x, y, z + 1, Nx64, Ny64, Nz64);
			uz.y = tvRead(u, x, y, z - 1, Nx64, Ny64, Nz64);
			ux = uijk - ux;
			uy = uijk - uy;
			uz = uijk - uz;
			float2<T> uabsx = ux / (fabs2(ux) + epps);
			float2<T> uabsy = uy / (fabs2(uy) + epps);
			float2<T> uabsz = uz / (fabs2(uz) + epps);
			float2<T> output = uabsx - uabsx / (fabs2(ux) / sigma + 1.f) + uabsy - uabsy / (fabs2(uy) / sigma + 1.f) + uabsz - uabsz / (fabs2(uz) / sigma + 1.f);
			grad[n] += beta * (output.x + output.y);
		}
		else {
			float3<T> uijkP;
			float3<T> uijkM;
			float2<T> ui;
			float2<T> uj;
			float2<T> uk;
			// Volume edge is clamped (Neumann boundary), identical to the GPU kernels
			uijkP.x = tvRead(u, x + 1, y, z, Nx64, Ny64, Nz64);
			uijkP.y = tvRead(u, x, y + 1, z, Nx64, Ny64, Nz64);
			uijkP.z = tvRead(u, x, y, z + 1, Nx64, Ny64, Nz64);
			uijkM.x = tvRead(u, x - 1, y, z, Nx64, Ny64, Nz64);
			uijkM.y = tvRead(u, x, y - 1, z, Nx64, Ny64, Nz64);
			uijkM.z = tvRead(u, x, y, z - 1, Nx64, Ny64, Nz64);

			ui.x = tvRead(u, x - 1, y + 1, z, Nx64, Ny64, Nz64);
			ui.y = tvRead(u, x - 1, y, z + 1, Nx64, Ny64, Nz64);
			uj.x = tvRead(u, x + 1, y - 1, z, Nx64, Ny64, Nz64);
			uj.y = tvRead(u, x, y - 1, z + 1, Nx64, Ny64, Nz64);
			uk.x = tvRead(u, x + 1, y, z - 1, Nx64, Ny64, Nz64);
			uk.y = tvRead(u, x, y + 1, z - 1, Nx64, Ny64, Nz64);
			float3<T> u1;
			float3<T> u2;
			float3<T> u3;
			MFLOAT3(uijk - uijkM.x, ui.x - uijkM.x, ui.y - uijkM.x, u1);
			MFLOAT3(uj.x - uijkM.y, uijk - uijkM.y, uj.y - uijkM.y, u2);
			MFLOAT3(uk.x - uijkM.z, uk.y - uijkM.z, uijk - uijkM.z, u3);
			if (type == 6) { // START TVW1
				// const float3 u4 = uijk - uijkM;
				float3<T> u4 = uijkP - uijk;
				float3<T> w4 = (u4) / sigma;
				w4 = exp3(-w4 * w4);
				const T pvalijk = sqrtVal(u4, epps, type, w4);
				float3<T> w1 = (u1) / sigma;
				w1 = exp3(-w1 * w1);
				float3<T> w2 = (u2) / sigma;
				w2 = exp3(-w2 * w2);
				float3<T> w3 = (u3) / sigma;
				w3 = exp3(-w3 * w3);
				grad[n] += beta * (-(w4.x * u4.x + w4.y * u4.y + w4.z * u4.z) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, type, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, type, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, type, w3));
			}
			else if (type == 1 && anatomical) {
				// phi_n = sqrt(Q_n(g_n) + epps) with the weights of voxel n, gradient = -(q_m.x + q_m.y + q_m.z) + q_{m-ex}.x + q_{m-ey}.y + q_{m-ez}.z
				// The image is clamped at the volume edge, so the differences that cross the edge are zero and do not contribute
				const float3<T> q0 = anatTVQ(S, x, y, z, Nx64, Ny64, Nz64, uijkP - uijk, epps);
				const float3<T> q1 = anatTVQ(S, x - 1, y, z, Nx64, Ny64, Nz64, u1, epps);
				const float3<T> q2 = anatTVQ(S, x, y - 1, z, Nx64, Ny64, Nz64, u2, epps);
				const float3<T> q3 = anatTVQ(S, x, y, z - 1, Nx64, Ny64, Nz64, u3, epps);
				T res = (T)0.;
				if (x < Nx64 - 1)
					res -= q0.x;
				if (y < Ny64 - 1)
					res -= q0.y;
				if (z < Nz64 - 1)
					res -= q0.z;
				if (x > 0)
					res += q1.x;
				if (y > 0)
					res += q2.y;
				if (z > 0)
					res += q3.z;
				grad[n] += beta * (res + (T)1e-7);
			}
			else if (type == 2 && anatomical) {
				// The reference image is clamped at the volume edge (same as the image)
				float3<T> uijkR;
				uijkR.x = tvRead(S, x + 1, y, z, Nx64, Ny64, Nz64);
				uijkR.y = tvRead(S, x, y + 1, z, Nx64, Ny64, Nz64);
				uijkR.z = tvRead(S, x, y, z + 1, Nx64, Ny64, Nz64);
				float3<T> apuS = (uijkR - S[n]);
				float3<T> apu = uijkP - uijk;
				const T pvalijk = std::sqrt(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
				// Reference image differences at the three backward neighbors, needed for their phi = sqrt(|g|^2 + C |r|^2 + epps)
				const T vC = tvRead(S, x, y, z, Nx64, Ny64, Nz64);
				const T vX = tvRead(S, x - 1, y, z, Nx64, Ny64, Nz64);
				const T vY = tvRead(S, x, y - 1, z, Nx64, Ny64, Nz64);
				const T vZ = tvRead(S, x, y, z - 1, Nx64, Ny64, Nz64);
				float3<T> r1, r2, r3;
				MFLOAT3(vC - vX, tvRead(S, x - 1, y + 1, z, Nx64, Ny64, Nz64) - vX, tvRead(S, x - 1, y, z + 1, Nx64, Ny64, Nz64) - vX, r1);
				MFLOAT3(tvRead(S, x + 1, y - 1, z, Nx64, Ny64, Nz64) - vY, vC - vY, tvRead(S, x, y - 1, z + 1, Nx64, Ny64, Nz64) - vY, r2);
				MFLOAT3(tvRead(S, x + 1, y, z - 1, Nx64, Ny64, Nz64) - vZ, tvRead(S, x, y + 1, z - 1, Nx64, Ny64, Nz64) - vZ, vC - vZ, r3);
				const T p1 = std::sqrt(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + C * (r1.x * r1.x + r1.y * r1.y + r1.z * r1.z) + epps);
				const T p2 = std::sqrt(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + C * (r2.x * r2.x + r2.y * r2.y + r2.z * r2.z) + epps);
				const T p3 = std::sqrt(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + C * (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) + epps);
				grad[n] += beta * (((T)3. * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / p1 + (uijk - uijkM.y) / p2 + (uijk - uijkM.z) / p3 + 1e-7f);
			}
			else if (type == 5 && anatomical) {
				// Asymmetric parallel level sets (Ehrhardt et al. 2016): P(u) = sum sqrt(|g|^2 - (xi.g)^2 + epps), xi = r / sqrt(|r|^2 + eta^2)
				// The gradient is the sum of the vectors q of the voxel itself and of its three backward neighbors
				// The reference image is clamped at the volume edge, so its differences vanish across the boundary
				const T vC = tvRead(S, x, y, z, Nx64, Ny64, Nz64);
				const T vX = tvRead(S, x - 1, y, z, Nx64, Ny64, Nz64);
				const T vY = tvRead(S, x, y - 1, z, Nx64, Ny64, Nz64);
				const T vZ = tvRead(S, x, y, z - 1, Nx64, Ny64, Nz64);
				float3<T> r0, r1, r2, r3;
				MFLOAT3(tvRead(S, x + 1, y, z, Nx64, Ny64, Nz64) - vC, tvRead(S, x, y + 1, z, Nx64, Ny64, Nz64) - vC, tvRead(S, x, y, z + 1, Nx64, Ny64, Nz64) - vC, r0);
				MFLOAT3(vC - vX, tvRead(S, x - 1, y + 1, z, Nx64, Ny64, Nz64) - vX, tvRead(S, x - 1, y, z + 1, Nx64, Ny64, Nz64) - vX, r1);
				MFLOAT3(tvRead(S, x + 1, y - 1, z, Nx64, Ny64, Nz64) - vY, vC - vY, tvRead(S, x, y - 1, z + 1, Nx64, Ny64, Nz64) - vY, r2);
				MFLOAT3(tvRead(S, x + 1, y, z - 1, Nx64, Ny64, Nz64) - vZ, tvRead(S, x, y + 1, z - 1, Nx64, Ny64, Nz64) - vZ, vC - vZ, r3);
				const float3<T> q0 = aplsQ(uijkP - uijk, r0, C, epps);
				const float3<T> q1 = aplsQ(u1, r1, C, epps);
				const float3<T> q2 = aplsQ(u2, r2, C, epps);
				const float3<T> q3 = aplsQ(u3, r3, C, epps);
				grad[n] += beta * (-(q0.x + q0.y + q0.z) + q1.x + q2.y + q3.z + 1e-7f);
			}
			else {
				const T pvalijk = sqrtVal(uijkP - uijk, epps, type);
				grad[n] += beta * (((T)3. * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps, type) + (uijk - uijkM.y) / sqrtVal(u2, epps, type) + (uijk - uijkM.z) / sqrtVal(u3, epps, type) + 1e-7f);
			}
		}
	}
}