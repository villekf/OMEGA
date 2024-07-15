#pragma once
#include "projector_functions.h"

template <typename T>
void RDPKernel(T* grad, const T* u,	const int Nx, const int Ny, const int Nz, const int NxOrig, const int NyOrig, const int NzOrig, const T gamma, const T epps, const T beta, 
	const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr) {

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
		// Current voxel
		float2<T> ux, uy, uz;
		const T uj = u[x + y * Nx + z * Nx * Ny];
		MFLOAT2(0.f, 0.f, ux);
		MFLOAT2(0.f, 0.f, uy);
		MFLOAT2(0.f, 0.f, uz);
		// Left-right
		if (x < Nx - 1)
			ux.x = u[(x + 1) + y * Nx + z * Nx * Ny];
		if (x > 0)
			ux.y = u[(x - 1) + y * Nx + z * Nx * Ny];
		// Top-bottom
		if (y < Ny - 1)
			uy.x = u[(x)+(y + 1) * Nx + (z)*Nx * Ny];
		if (y > 0)
			uy.y = u[(x)+(y - 1) * Nx + (z)*Nx * Ny];
		// Front-back
		if (z < Nz - 1)
			uz.x = u[(x)+(y)*Nx + (z + 1) * Nx * Ny];
		if (z > 0)
			uz.y = u[(x)+(y)*Nx + (z - 1) * Nx * Ny];
		float2<T> uxy;
		float2<T> uyx;
		float2<T> uxz;
		float2<T> uzx;
		float2<T> uyz;
		float2<T> uzy;
		//if (x < Nx - 1 && y < Ny - 1)
		//	uxy.x = u[(x + 1) + (y + 1) * Nx + z * Nx * Ny];
		//if (x > 0 && y < Ny - 1)
		//	uxy.y = u[(x - 1) + (y + 1) * Nx + z * Nx * Ny];
		//if (y > 0 && x < Nx - 1)
		//	uyx.x = u[(x + 1) + (y - 1) * Nx + (z) * Nx * Ny];
		//if (y > 0 && x > 0)
		//	uyx.y = u[(x - 1) + (y - 1) * Nx + (z) * Nx * Ny];
		//if (z < N.z - 1 && x < Nx - 1)
		//	uxz.x = u[(x + 1) + (y) * Nx + (z + 1) * Nx * Ny];
		//if (z < N.z - 1 && x > 0)
		//	uxz.y = u[(x - 1) + (y) * Nx + (z + 1) * Nx * Ny];
		//if (x < Nx - 1 && z > 0)
		//	uzx.x = u[(x + 1) + (y) * Nx + (z - 1) * Nx * Ny];
		//if (x > 0 && z > 0)
		//	uzx.y = u[(x - 1) + (y) * Nx + (z - 1) * Nx * Ny];
		//if (y < Ny - 1 && z < N.z - 1)
		//	uyz.x = u[(x) + (y + 1) * Nx + (z + 1) * Nx * Ny];
		//if (y > 0 && z < N.z - 1)
		//	uyz.y = u[(x) + (y - 1) * Nx + (z + 1) * Nx * Ny];
		//if (z > 0 && y < Ny - 1)
		//	uzy.x = u[(x) + (y + 1) * Nx + (z - 1) * Nx * Ny];
		//if (y > 0 && z > 0)
		//	uzy.y = u[(x) + (y - 1) * Nx + (z - 1) * Nx * Ny];
		float2<T> uuxy = (uj - uxy);
		float2<T> uuyx = (uj - uyx);
		float2<T> uuxz = (uj - uxz);
		float2<T> uuzx = (uj - uzx);
		float2<T> uuyz = (uj - uyz);
		float2<T> uuzy = (uj - uzy);
		float2<T> divPow2XY = (uj + uuxy + gamma * fabs2(uuxy));
		float2<T> divPow2YX = (uj + uuyx + gamma * fabs2(uuyx));
		float2<T> divPow2XZ = (uj + uuxz + gamma * fabs2(uuxz));
		float2<T> divPow2ZX = (uj + uuzx + gamma * fabs2(uuzx));
		float2<T> divPow2YZ = (uj + uuyz + gamma * fabs2(uuyz));
		float2<T> divPow2ZY = (uj + uuzy + gamma * fabs2(uuzy));
		float2<T> uj_ux = uj - ux;
		float2<T> uj_uy = uj - uy;
		float2<T> uj_uz = uj - uz;
		float2<T> divPow2X = (uj + ux + gamma * fabs2(uj_ux));
		float2<T> divPow2Y = (uj + uy + gamma * fabs2(uj_uy));
		float2<T> divPow2Z = (uj + uz + gamma * fabs2(uj_uz));
		float2<T> output = uj_ux * (gamma * fabs2(uj_ux) + uj + 3.f * ux + epps * epps) / (divPow2X * divPow2X + epps)
			+ uj_uy * (gamma * fabs2(uj_uy) + uj + 3.f * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
			+ uj_uz * (gamma * fabs2(uj_uz) + uj + 3.f * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
		//output += uuxy * (gamma * fabs2(uuxy) + uj + 3.f * uxy + epps * epps) / (divPow2XY * divPow2XY + epps) * M_SQRT1_2_F + uuyx * (gamma * fabs2(uuyx) + uj + 3.f * uyx + epps * epps) / (divPow2YX * divPow2YX + epps) * M_SQRT1_2_F +
		//	uuxz * (gamma * fabs2(uuxz) + uj + 3.f * uxz + epps * epps) / (divPow2XZ * divPow2XZ + epps) * M_SQRT1_2_F + uuzx * (gamma * fabs2(uuzx) + uj + 3.f * uzx + epps * epps) / (divPow2ZX * divPow2ZX + epps) * M_SQRT1_2_F +
		//	uuyz * (gamma * fabs2(uuyz) + uj + 3.f * uyz + epps * epps) / (divPow2YZ * divPow2YZ + epps) * M_SQRT1_2_F + uuzy * (gamma * fabs2(uuzy) + uj + 3.f * uzy + epps * epps) / (divPow2ZY * divPow2ZY + epps) * M_SQRT1_2_F;
		if (std::isnan(output.x))
			output.x = (T)0.;
		if (std::isnan(output.y))
			output.y = (T)0.;
		grad[n] += beta * (output.x + output.y);
	}
}