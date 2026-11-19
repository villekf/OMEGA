#pragma once
#include "projector_functions.h"

// Neighbor read, each coordinate is clamped to the volume (clamp-to-edge, i.e. Neumann boundary), same as the GPU kernels
template <typename T>
inline T rdpRead(const T* buf, const int64_t x, const int64_t y, const int64_t z, const int64_t Nx, const int64_t Ny, const int64_t Nz) {
	const int64_t xc = std::max<int64_t>(0, std::min<int64_t>(x, Nx - 1));
	const int64_t yc = std::max<int64_t>(0, std::min<int64_t>(y, Ny - 1));
	const int64_t zc = std::max<int64_t>(0, std::min<int64_t>(z, Nz - 1));
	return buf[xc + yc * Nx + zc * Nx * Ny];
}

template <typename T>
void RDPKernel(T* grad, const T* u,	const int Nx, const int Ny, const int Nz, const int NxOrig, const int NyOrig, const int NzOrig, const T gamma, const T epps, const T beta, 
	const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr, const uint32_t maskBPZ = 1U) {

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
		// Voxels outside the extended FOV or the prior mask do not contribute, same as the GPU kernel
		if (fovIndices != nullptr && fovIndices[z] == 0)
			continue;
		if (maskBP != nullptr) {
			const int64_t maskInd = maskBPZ > 1U ? n : x + y * static_cast<int64_t>(Nx);
			if (maskBP[maskInd] == 0)
				continue;
		}
		// Current voxel
		float2<T> ux, uy, uz;
		const int64_t Nx64 = static_cast<int64_t>(Nx);
		const int64_t Ny64 = static_cast<int64_t>(Ny);
		const int64_t Nz64 = static_cast<int64_t>(Nz);
		const T uj = u[n];
		// Neighbors outside the volume are clamped to the edge (Neumann boundary), identical to the GPU kernel
		// Left-right
		ux.x = rdpRead(u, x + 1, y, z, Nx64, Ny64, Nz64);
		ux.y = rdpRead(u, x - 1, y, z, Nx64, Ny64, Nz64);
		// Top-bottom
		uy.x = rdpRead(u, x, y + 1, z, Nx64, Ny64, Nz64);
		uy.y = rdpRead(u, x, y - 1, z, Nx64, Ny64, Nz64);
		// Front-back
		uz.x = rdpRead(u, x, y, z + 1, Nx64, Ny64, Nz64);
		uz.y = rdpRead(u, x, y, z - 1, Nx64, Ny64, Nz64);
		// TODO: Potentially implement RDP with corners later
		//float2<T> uxy;
		//float2<T> uyx;
		//float2<T> uxz;
		//float2<T> uzx;
		//float2<T> uyz;
		//float2<T> uzy;
		//float2<T> uuxy = (uj - uxy);
		//float2<T> uuyx = (uj - uyx);
		//float2<T> uuxz = (uj - uxz);
		//float2<T> uuzx = (uj - uzx);
		//float2<T> uuyz = (uj - uyz);
		//float2<T> uuzy = (uj - uzy);
		//float2<T> divPow2XY = (uj + uuxy + gamma * fabs2(uuxy));
		//float2<T> divPow2YX = (uj + uuyx + gamma * fabs2(uuyx));
		//float2<T> divPow2XZ = (uj + uuxz + gamma * fabs2(uuxz));
		//float2<T> divPow2ZX = (uj + uuzx + gamma * fabs2(uuzx));
		//float2<T> divPow2YZ = (uj + uuyz + gamma * fabs2(uuyz));
		//float2<T> divPow2ZY = (uj + uuzy + gamma * fabs2(uuzy));
		float2<T> uj_ux = uj - ux;
		float2<T> uj_uy = uj - uy;
		float2<T> uj_uz = uj - uz;
		float2<T> divPow2X = (uj + ux + gamma * fabs2(uj_ux));
		float2<T> divPow2Y = (uj + uy + gamma * fabs2(uj_uy));
		float2<T> divPow2Z = (uj + uz + gamma * fabs2(uj_uz));
		float2<T> output = uj_ux * (gamma * fabs2(uj_ux) + uj + (T)3. * ux + epps * epps) / (divPow2X * divPow2X + epps)
			+ uj_uy * (gamma * fabs2(uj_uy) + uj + (T)3. * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
			+ uj_uz * (gamma * fabs2(uj_uz) + uj + (T)3. * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
		if (std::isnan(output.x))
			output.x = (T)0.;
		if (std::isnan(output.y))
			output.y = (T)0.;
		grad[n] += beta * (output.x + output.y);
	}
}