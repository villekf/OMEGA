/**************************************************************************
* Header for Metal error compatibility constants
*
* Copyright(C) 2026 Niilo Saarlemo
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
#include <iostream>
#include <simd/simd.h>
#include "kernelParams.hpp"

#ifndef CL_SUCCESS
#define CL_SUCCESS 0
#endif
#ifndef CL_MEM_READ_ONLY
#define CL_MEM_READ_ONLY 0
#endif
#ifndef CL_MEM_WRITE_ONLY
#define CL_MEM_WRITE_ONLY 0
#endif
#ifndef CL_MEM_READ_WRITE
#define CL_MEM_READ_WRITE 0
#endif

// For 64-bit integer atomics
#define TH 100000000000.f
// For 32-bit integer atomics
#define TH32 100000.f
// For projector type 4 backprojection
// How many voxels are computed per thread
#define NVOXELS 1
#define NVOXELSHELICAL 1
// How many voxels are computed per thread for BDD backprojection
#define NVOXELS5 1
// How many slices are computed in the same thread for BDD forward projection
#define NVOXELSFP 8

#define getErrorString(ans) { metalAssert((ans), __FILE__, __LINE__); }

inline void metalAssert(const int code, const char* file, const int line)
{
	if (code != CL_SUCCESS)
		std::cerr << "Metal assert: " << code << ", " << file << ", line " << line << std::endl;
}
