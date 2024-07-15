#pragma once
#include "projector_functions.h"

template <typename T>
T sqrtVal(const float3<T> input, const T epps, const int type, const float3<T> w = float3<T>{}) {
if (type == 6)
	return std::sqrt(w.x * input.x * input.x + w.y * input.y * input.y + w.z * input.z * input.z + epps);
else
	return std::sqrt(input.x * input.x + input.y * input.y + input.z * input.z + epps);
}

template <typename T>
void TVKernel(T* grad, const T* u, const int Nx, const int Ny, const int Nz, const int NxOrig, const int NyOrig, const int NzOrig, const T sigma, const T epps, const T beta, 
	const int type, const bool anatomical, const T C = 0.f, const T* S = nullptr, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr) {

	setThreads();
	int64_t start = 0;
	int64_t end = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny) * static_cast<int64_t>(Nz);
	int64_t Nxy = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny);
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
		//if (fovIndices[z] == 0)
		//	return;
		//const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(x, y)).w;
		//if (maskVal == 0)
		//	return;
		const T uijk = u[n];
		if (type == 4) { // START JPTV || SATV
			float2<T> ux;
			float2<T> uy;
			float2<T> uz;
			if (x < Nx - 1)
				ux.x = u[(x + 1) + (y)*Nx + (z)*Nxy];
			if (x > 0)
				ux.y = u[(x - 1) + (y)*Nx + (z)*Nxy];
			if (y < Ny - 1)
				uy.x = u[(x)+(y + 1) * Nx + (z)*Nxy];
			if (y > 0)
				uy.y = u[(x)+(y - 1) * Nx + (z)*Nxy];
			if (z < Nz - 1)
				uz.x = u[(x)+(y)*Nx + (z + 1) * Nxy];
			if (z > 0)
				uz.y = u[(x)+(y)*Nx + (z - 1) * Nxy];
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
			if (x < Nx - 1)
				uijkP.x = u[(x + 1) + (y)*Nx + (z)*Nxy];
			if (y < Ny - 1)
				uijkP.y = u[(x)+(y + 1) * Nx + (z)*Nxy];
			if (z < Nz - 1)
				uijkP.z = u[(x)+(y)*Nx + (z + 1) * Nxy];
			if (x > 0)
				uijkM.x = u[(x - 1) + (y)*Nx + (z)*Nxy];
			if (y > 0)
				uijkM.y = u[(x)+(y - 1) * Nx + (z)*Nxy];
			if (z > 0)
				uijkM.z = u[(x)+(y)*Nx + (z - 1) * Nxy];

			if (x > 0 && y < Ny - 1)
				ui.x = u[(x - 1) + (y + 1) * Nx + (z)*Nxy];
			if (x > 0 && z < Nz - 1)
				ui.y = u[(x - 1) + (y)*Nx + (z + 1) * Nxy];
			if (y > 0 && x < Nx - 1)
				uj.x = u[(x + 1) + (y - 1) * Nx + (z)*Nxy];
			if (y > 0 && z < Nz - 1)
				uj.y = u[(x)+(y - 1) * Nx + (z + 1) * Nxy];
			if (z > 0 && x < Nx - 1)
				uk.x = u[(x + 1) + (y)*Nx + (z - 1) * Nxy];
			if (z > 0 && y < Ny - 1)
				uk.y = u[(x)+(y + 1) * Nx + (z - 1) * Nxy];
			//const T pvalijk = distance(uijkP, uijkM) + eppsSqrt;
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
			if (type == 1 && anatomical) {
				const int64_t NN = Nxy * Nz;
				T s[9];
				for (int kk = 0; kk < 9; kk++)
					s[kk] = S[n + NN * kk];
				// s1*(f - fi)^2 + s5*(f - fj)^2 + s9*(f - fk)^2 + s2*(f - fi)*(f - fj) + s4*(f - fi)*(f - fj) + s3*(f - fi)*(f - fk) + s7*(f - fi)*(f - fk) + s6*(f - fj)*(f - fk) + s8*(f - fj)*(f - fk)
				float3<T> val = uijkP - uijk;
				const T pvalijk = sqrt(val.x * val.x * s[0] + val.y * val.y * s[4] + val.z * val.z * s[8] + s[1] * (val.x) * (val.y) + s[3] * (val.x) * (val.y) + s[2] * (val.x) * (val.z) + s[6] * (val.x) * (val.z) +
					s[5] * (val.y) * (val.z) + s[7] * (val.y) * (val.z) + epps);
				const T pvalijkX = sqrt(u1.x * u1.x * s[0] + u1.y * u1.y * s[4] + u1.z * u1.z * s[8] + s[1] * (u1.x) * (u1.y) + s[3] * (u1.x) * (u1.y) + s[2] * (u1.x) * (u1.z) + s[6] * (u1.x) * (u1.z) +
					s[5] * (u1.y) * (u1.z) + s[7] * (u1.y) * (u1.z) + epps);
				const T pvalijkY = sqrt(u2.x * u2.x * s[0] + u2.y * u2.y * s[4] + u2.z * u2.z * s[8] + s[1] * (u2.x) * (u2.y) + s[3] * (u2.x) * (u2.y) + s[2] * (u2.x) * (u2.z) + s[6] * (u2.x) * (u2.z) +
					s[5] * (u2.y) * (u2.z) + s[7] * (u2.y) * (u2.z) + epps);
				const T pvalijkZ = sqrt(u3.x * u3.x * s[0] + u3.y * u3.y * s[4] + u3.z * u3.z * s[8] + s[1] * (u3.x) * (u3.y) + s[3] * (u3.x) * (u3.y) + s[2] * (u3.x) * (u3.z) + s[6] * (u3.x) * (u3.z) +
					s[5] * (u3.y) * (u3.z) + s[7] * (u3.y) * (u3.z) + epps);
				const T dx = s[0] * (2.f * (uijk - uijkM.x)) + s[3] * u1.y + s[2] * u1.z + s[6] * u1.z + s[1] * u1.y;
				const T dy = s[4] * (2.f * (uijk - uijkM.y)) + s[5] * u2.z + s[3] * u2.x + s[1] * u2.x + s[7] * u2.z;
				const T dz = s[8] * (2.f * (uijk - uijkM.z)) + s[6] * u3.x + s[5] * u3.y + s[7] * u3.y + s[2] * u3.x;
				const T d = s[1] * val.x + s[2] * val.x + s[3] * val.x + s[6] * val.x + s[1] * val.y + s[3] * val.y + s[5] * val.y + s[7] * val.y + s[2] * val.z + s[5] * val.z + s[6] * val.z + s[7] * val.z + s[0] * 2.f * val.x + s[4] * 2.f * val.y + s[8] * 2.f * val.z;
				// d = s2*(f - fi) + s3*(f - fi) + s4*(f - fi) + s7*(f - fi) + s2*(f - fj) + s4*(f - fj) + s6*(f - fj) + s8*(f - fj) + s3*(f - fk) + s6*(f - fk) + s7*(f - fk) + s8*(f - fk) + s1*(2*f - 2*fi) + s5*(2*f - 2*fj) + s9*(2*f - 2*fk)
				// dx = s2*(fij - fji) + s4*(fij - fji) + s3*(fik - fki) + s7*(fik - fki) + s1*(2*f - 2*fi2)
				// dy = s6*(fjk - fkj) - s4*(fij - fji) - s2*(fij - fji) + s8*(fjk - fkj) + s5*(2*f - 2*fj2)
				// dz = s9*(2*f - 2*fk2) - s7*(fik - fki) - s6*(fjk - fkj) - s8*(fjk - fkj) - s3*(fik - fki)
				grad[n] += beta * .5f * (d / pvalijk + dx / pvalijkX + dy / pvalijkY + dz / pvalijkZ);
			}
			else if (type == 2 && anatomical) {
				float3<T> uijkR;
				if (x < Nx - 1)
					uijkR.x = S[(x + 1) + (y)*Nx + (z)*Nxy];
				if (y < Ny - 1)
					uijkR.y = S[(x)+(y + 1) * Nx + (z)*Nxy];
				if (z < Nz - 1)
					uijkR.z = S[(x)+(y)*Nx + (z + 1) * Nxy];
				float3<T> apuS = (uijkR - S[n]);
				float3<T> apu = uijkP - uijk;
				const T pvalijk = std::sqrt(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
				grad[n] += beta * (((T)3. * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps, type) + (uijk - uijkM.y) / sqrtVal(u2, epps, type) + (uijk - uijkM.z) / sqrtVal(u3, epps, type) + 1e-7f);
			}
			if (type == 5 && anatomical) {
				float3<T> uijkR;
				if (x < Nx - 1)
					uijkR.x = S[(x + 1) + (y)*Nx + (z)*Nxy];
				if (y < Ny - 1)
					uijkR.y = S[(x)+(y + 1) * Nx + (z)*Nxy];
				if (z < Nz - 1)
					uijkR.z = S[(x)+(y)*Nx + (z + 1) * Nxy];
				float3<T> epsilon = (uijkR - S[n]);
				epsilon = epsilon / std::sqrt(epsilon.x * epsilon.x + epsilon.y * epsilon.y + epsilon.z * epsilon.z + C * C);
				float3<T> apu = uijkP - uijk;
				const T apuR = uijkR.x * apu.x + uijkR.y * apu.y + uijkR.z * apu.z;
				// (vx*(f - fi) + vy*(f - fj) + vz*(f - fk) + (f - fi)^2 + (f - fj)^2 + (f - fk)^2
				const T pvalijk = std::sqrt(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z - apuR * apuR + epps);
				T apuRXYZ = uijkR.x * u1.x + uijkR.y * u1.y + uijkR.z * u1.z;
				const T pvalijkX = std::sqrt(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + apuRXYZ * apuRXYZ + epps);
				apuRXYZ = uijkR.x * u2.x + uijkR.y * u2.y + uijkR.z * u2.z;
				const T pvalijkY = std::sqrt(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + apuRXYZ * apuRXYZ + epps);
				apuRXYZ = uijkR.x * u3.x + uijkR.y * u3.y + uijkR.z * u3.z;
				const T pvalijkZ = std::sqrt(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + apuRXYZ * apuRXYZ + epps);
				// d = -(2*fi - 6*f + 2*fj + 2*fk + 2*(vx*(f - fi) + vy*(f - fj) + vz*(f - fk))*(vx + vy + vz))
				// dx = 2*f - 2*fi2 - 2*vx*(vx*(f - fi2) + vy*(fij - fji) + vz*(fik - fki))
				grad[n] += beta * .5f * ((6.f * uijk - 2.f * uijkP.x - 2.f * uijkP.y - 2.f * uijkP.z + 2.f * (epsilon.x * (uijk - uijkP.x) + epsilon.y * (uijk - uijkP.y) + epsilon.z * (uijk - uijkP.z)) * (epsilon.x + epsilon.y + epsilon.z)) / pvalijk +
					2.f * (u1.x - epsilon.x * (epsilon.x * u1.x + epsilon.y * u1.y + epsilon.z * u1.z)) / pvalijkX + 2.f * (u2.y - epsilon.y * (epsilon.x * u2.x + epsilon.y * u2.y + epsilon.z * u2.z)) / pvalijkY +
					2.f * (u3.z - epsilon.z * (epsilon.x * u3.x + epsilon.y * u3.y + epsilon.z * u3.z)) / pvalijkZ + 1e-7f);
			}
			else {
				const T pvalijk = sqrtVal(uijkP - uijk, epps, type);
				grad[n] += beta * (((T)3. * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps, type) + (uijk - uijkM.y) / sqrtVal(u2, epps, type) + (uijk - uijkM.z) / sqrtVal(u3, epps, type) + 1e-7f);
			}
		}
	}
}