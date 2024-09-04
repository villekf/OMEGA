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
		if (std::isnan(output.x))
			output.x = (T)0.;
		if (std::isnan(output.y))
			output.y = (T)0.;
		grad[n] += beta * (output.x + output.y);
	}
}