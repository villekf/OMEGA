/**************************************************************************
* Header for ArrayFire functions. Used by implementation 2.
*
* Copyright(C) 2019-2024 Ville - Veikko Wettenhovi
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
#define NOMINMAX
#ifdef CUDA
#include "ProjectorClassCUDA.h"
#define transferAF(varA) varA.device<CUdeviceptr>()
#elif defined(CPU)
#include "ProjectorClassCPU.h"
#define transferAF(varA) varA.device<float>()
#else
#include "ProjectorClass.h"
#define transferAF(varA) cl::Buffer(*varA.device<cl_mem>(), true)
#endif
#define VAL 0.00001f
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#pragma pack(1) 
#pragma warning(disable : 4996)

// This function sets the variables needed for the special large dimensional reconstruction method
inline void largeDimCreate(scalarStruct& inputScalars, const RecMethods& MethodList, const Weighting& w_vec) {
	inputScalars.lDimStruct.Nz.resize(inputScalars.subsets);
	inputScalars.lDimStruct.NzPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.startPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.endPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.imDim.resize(inputScalars.subsets);
	inputScalars.lDimStruct.imDimPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.cumDim.resize(inputScalars.subsets + 1);
	inputScalars.lDimStruct.cumDimPr.resize(inputScalars.subsets + 1);
	inputScalars.lDimStruct.bz.resize(inputScalars.subsets);
	inputScalars.lDimStruct.bmaxZ.resize(inputScalars.subsets);
	inputScalars.lDimStruct.d_Scale4Z.resize(inputScalars.subsets);
	inputScalars.lDimStruct.cumDim[0] = 0LL;
	inputScalars.lDimStruct.cumDimPr[0] = 0LL;
	const uint32_t intZ = inputScalars.Nz[0] / inputScalars.subsets;
	const uint32_t remainder = inputScalars.Nz[0] % inputScalars.subsets;
	for (int kk = 0; kk < inputScalars.subsets; kk++) {
		if (kk == 0) {
			inputScalars.lDimStruct.Nz[kk] = intZ + remainder;
			inputScalars.lDimStruct.bz[kk] = inputScalars.bz[0];
		}
		else {
			inputScalars.lDimStruct.Nz[kk] = intZ;
			inputScalars.lDimStruct.bz[kk] = inputScalars.lDimStruct.bmaxZ[kk - 1];
		}
		inputScalars.lDimStruct.bmaxZ[kk] = inputScalars.lDimStruct.bz[kk] + inputScalars.lDimStruct.Nz[kk] * inputScalars.dz[0];
		if (MethodList.NLM || MethodList.RDP || MethodList.hyperbolic || MethodList.GGMRF || MethodList.TV) {
			inputScalars.lDimStruct.imDim[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk]);
			inputScalars.lDimStruct.cumDim[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk];
			if (kk < inputScalars.subsetsUsed - 1 && kk > 0) {
				if (MethodList.NLM) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U + w_vec.Nlz * 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U + w_vec.Nlz * 2U;
				}
				else if (MethodList.TV && (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + 2U;
				}
				else {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U;
				}
			}
			else {
				if (MethodList.NLM) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz + w_vec.Nlz));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz + w_vec.Nlz;
				}
				else if (MethodList.TV && (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + 1));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + 1;
				}
				else {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz;
				}
			}
			if (kk == 0) {
				if (MethodList.NLM)
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (w_vec.Ndz + w_vec.Nlz);
				else if (MethodList.TV && (MethodList.RDP && !w_vec.RDPLargeNeighbor))
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * w_vec.Ndz;
				inputScalars.lDimStruct.startPr[kk] = 0;
			}
			else
				if (MethodList.NLM) {
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (w_vec.Ndz + w_vec.Nlz);
					inputScalars.lDimStruct.startPr[kk] = (w_vec.Ndz + w_vec.Nlz) * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				}
				else if (MethodList.TV && (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.startPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				}
				else {
					inputScalars.lDimStruct.startPr[kk] = w_vec.Ndz * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * w_vec.Ndz;
				}
			if (kk == inputScalars.subsetsUsed - 1)
				inputScalars.lDimStruct.endPr[kk] = 0U;
			else
				if (MethodList.NLM)
					inputScalars.lDimStruct.endPr[kk] = (w_vec.Ndz + w_vec.Nlz) * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else if (MethodList.TV && (MethodList.RDP && !w_vec.RDPLargeNeighbor))
					inputScalars.lDimStruct.endPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else
					inputScalars.lDimStruct.endPr[kk] = w_vec.Ndz * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
		}
		else {
			inputScalars.lDimStruct.imDim[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk]);
			inputScalars.lDimStruct.cumDim[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk];
			inputScalars.lDimStruct.imDimPr[kk] = inputScalars.lDimStruct.imDim[kk];
			inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk + 1];
			inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk];
		}
		inputScalars.lDimStruct.d_Scale4Z[kk] = 1.f / (static_cast<float>(inputScalars.lDimStruct.Nz[kk]) * inputScalars.dz[0]);
		if (DEBUG) {
			mexPrintBase("inputScalars.lDimStruct.NzPr[kk] = %u\n", inputScalars.lDimStruct.NzPr[kk]);
			mexPrintBase("inputScalars.lDimStruct.Nz[kk] = %u\n", inputScalars.lDimStruct.Nz[kk]);
			mexPrintBase("inputScalars.lDimStruct.imDim[kk] = %u\n", inputScalars.lDimStruct.imDim[kk]);
			mexPrintBase("inputScalars.lDimStruct.imDimPr[kk] = %u\n", inputScalars.lDimStruct.imDimPr[kk]);
			mexEval();
		}
	}
}

inline void largeDimFirst(scalarStruct& inputScalars, ProjectorClass& proj, const int iter) {
	if (iter == 0) {
		inputScalars.lDimStruct.NzOrig = inputScalars.Nz[0];
		inputScalars.lDimStruct.imDimOrig = inputScalars.im_dim[0];
#if defined(OPENCL)
		inputScalars.lDimStruct.bzOrig = proj.b[0].s[2];
		inputScalars.lDimStruct.bmaxZOrig = proj.bmax[0].s[2];
		inputScalars.lDimStruct.d_Scale4ZOrig = inputScalars.d_Scale4[0].s[2];
#elif defined(CUDA)
		inputScalars.lDimStruct.bzOrig = proj.b[0].z;
		inputScalars.lDimStruct.bmaxZOrig = proj.bmax[0].z;
		inputScalars.lDimStruct.d_Scale4ZOrig = inputScalars.d_Scale4[0].z;
#endif
	}
	inputScalars.Nz[0] = inputScalars.lDimStruct.Nz[iter];
	inputScalars.im_dim[0] = inputScalars.lDimStruct.imDim[iter];
#if defined(OPENCL)
	proj.d_N[0].s[2] = inputScalars.lDimStruct.Nz[iter];
	proj.b[0].s[2] = inputScalars.lDimStruct.bz[iter];
	proj.bmax[0].s[2] = inputScalars.lDimStruct.bmaxZ[iter];
	inputScalars.d_Scale4[0].s[2] = inputScalars.lDimStruct.d_Scale4Z[iter];
#elif defined(CUDA)
	proj.d_N[0].z = inputScalars.lDimStruct.Nz[iter];
	proj.b[0].z = inputScalars.lDimStruct.bz[iter];
	proj.bmax[0].z = inputScalars.lDimStruct.bmaxZ[iter];
	inputScalars.d_Scale4[0].z = inputScalars.lDimStruct.d_Scale4Z[iter];
#endif
}

inline void largeDimLast(scalarStruct& inputScalars, ProjectorClass& proj) {
	inputScalars.Nz[0] = inputScalars.lDimStruct.NzOrig;
	inputScalars.im_dim[0] = inputScalars.lDimStruct.imDimOrig;
#if defined(OPENCL)
	proj.d_N[0].s[2] = inputScalars.lDimStruct.NzOrig;
	proj.b[0].s[2] = inputScalars.lDimStruct.bzOrig;
	proj.bmax[0].s[2] = inputScalars.lDimStruct.bmaxZOrig;
	inputScalars.d_Scale4[0].s[2] = inputScalars.lDimStruct.d_Scale4ZOrig;
#elif defined(CUDA)
	proj.d_N[0].z = inputScalars.lDimStruct.NzOrig;
	proj.b[0].z = inputScalars.lDimStruct.bzOrig;
	proj.bmax[0].z = inputScalars.lDimStruct.bmaxZOrig;
	inputScalars.d_Scale4[0].z = inputScalars.lDimStruct.d_Scale4ZOrig;
#endif
}

// Fill the backprojection (aka right-hand side) with zeros
inline void initializeRHS(AF_im_vectors& vec, const scalarStruct& inputScalars, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Initialize the backprojection output");
#ifdef OPENCL
	if (inputScalars.atomic_64bit)
		vec.rhs_os[ii] = (af::constant(0LL, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1, s64));
	else if (inputScalars.atomic_32bit)
		vec.rhs_os[ii] = (af::constant(0, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1, s32));
	else
#endif
		vec.rhs_os[ii] = (af::constant(0.f, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1));
	vec.rhs_os[ii].eval();
}

// Zero and symmetric padding of an array
inline af::array padding(const af::array& im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz,
	const bool zero_pad = false, const uint32_t Nw = 1)
{
	af::array out = im;
	if (zero_pad == 1) {
		af::dtype type = out.type();
		if (Nz == 1) {
			if (out.dims(1) == 1)
				out = moddims(out, Nx, Ny, Nz, Nw);
			out = af::constant(0, out.dims(0) + 2 * Ndx, out.dims(1) + 2 * Ndy, 1, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))), af::span, af::span) = out;
		}
		else {
			if (out.dims(2) == 1)
				out = moddims(out, Nx, Ny, Nz, Nw);
			af::array out = af::constant(0, out.dims(0) + 2 * Ndx, out.dims(1) + 2 * Ndy, out.dims(2) + 2 * Ndz, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))),
				static_cast<double>(Ndz) + af::seq(static_cast<double>(out.dims(2))), af::span) = out;
		}
	}
	else {
		if (out.dims(1) == 1)
			out = moddims(out, Nx, Ny, Nz, Nw);
		if (Ndx > 0)
			out = af::join(0, af::flip(out(af::seq(static_cast<double>(Ndx)), af::span, af::span, af::span), 0), out,
				af::flip(out(af::seq(static_cast<double>(out.dims(0) - Ndx), static_cast<double>(out.dims(0) - 1)), af::span, af::span, af::span), 0));
		if (Ndy > 0)
			out = af::join(1, af::flip(out(af::span, af::seq(static_cast<double>(Ndy)), af::span, af::span), 1), out,
				af::flip(out(af::span, af::seq(static_cast<double>(out.dims(1) - Ndy), static_cast<double>(out.dims(1) - 1)), af::span, af::span), 1));
		if (Nz == 1 || Ndz == 0) {
		}
		else {
			out = af::join(2, af::flip(out(af::span, af::span, af::seq(static_cast<double>(Ndz)), af::span), 2), out,
				af::flip(out(af::span, af::span, af::seq(static_cast<double>(out.dims(2) - Ndz), static_cast<double>(out.dims(2) - 1)), af::span), 2));
		}
	}
	return out;
}

// Compute the convolution of an array
inline af::array computeConvolution(const af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec,
	const uint32_t nRekos = 1, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Starting PSF blurring");
	af::array apu = af::moddims(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], nRekos);
	padding(apu, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1,
		false, inputScalars.nRekos2);
	apu = af::convolve3(apu, g);
	if (inputScalars.verbose >= 3)
		mexPrint("PSF blurring complete");
	return af::flat(apu);
}

/// <summary>
/// Transfer the sensitivity image from an ArrayFire array to an OpenCL buffer or CUDA device pointer
/// </summary>
/// <param name="apuSum the sensitivity image"></param>
/// <param name="proj the projector class object"></param>
/// <returns></returns>
inline void transferSensitivityImage(af::array& apuSum, ProjectorClass& proj) {
	apuSum.eval();
	af::sync();
	if (proj.d_Summ.size() < 1)
		proj.d_Summ.emplace_back(transferAF(apuSum));
	else
		proj.d_Summ[0] = transferAF(apuSum);
}

/// <summary>
/// Transfer the backprojection from an ArrayFire array to an OpenCL buffer or CUDA device pointer
/// </summary>
/// <param name="rhs_os the backprojection"></param>
/// <param name="proj the projector class object"></param>
/// <returns></returns>
inline int transferRHS(af::array& rhs_os, ProjectorClass& proj) {
	af::sync();
	if (DEBUG) {
		mexPrintBase("proj.vec_opencl.d_rhs_os.size() = %u\n", proj.vec_opencl.d_rhs_os.size());
		mexEval();
	}
	if (proj.vec_opencl.d_rhs_os.size() < 1)
		proj.vec_opencl.d_rhs_os.emplace_back(transferAF(rhs_os));
	else
		proj.vec_opencl.d_rhs_os[0] = transferAF(rhs_os);
	if (DEBUG) {
		mexPrintBase("proj.vec_opencl.d_rhs_os.size() = %u\n", proj.vec_opencl.d_rhs_os.size());
		mexEval();
	}
	return 0;
}

/// <summary>
/// Copy the current estimates from an ArrayFire array to an OpenCL image/CUDA texture. For branchless distance-driven, compute the integral images before copy
/// </summary>
/// <param name="vec image estimates and backprojection"></param>
/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
/// <param name="proj the projector class object"></param>
/// <param name="ii optional multi-resolution volume number, default is 0 (main volume)"></param>
/// <returns></returns>
inline int updateInputs(AF_im_vectors& vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii = 0) {

#ifdef CUDA
	CUresult status = CUDA_SUCCESS;
#elif defined(OPENCL)
	int status = 0;
	cl::detail::size_t_array region = { inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii] };
#endif
	if (inputScalars.FPType == 5) {
		af::array im;
		af::sync();
		af::deviceGC();
		mexEval();
		af::array intIm = af::constant(0.f, inputScalars.Ny[ii] + 1, inputScalars.Nz[ii] + 1, inputScalars.Nx[ii]);
		if (inputScalars.meanFP) {
			im = af::reorder(af::moddims(vec.im_os[ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 1, 2, 0);
			vec.meanFP = af::constant(0.f, inputScalars.Nx[ii] + inputScalars.Ny[ii]);
			vec.meanFP(af::seq(0, inputScalars.Nx[ii] - 1)) = af::flat(af::mean(af::mean(im, 0), 1));
			im -= af::tile(vec.meanFP(af::seq(0, inputScalars.Nx[ii] - 1)), im.dims(0), im.dims(1), 1);
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
			im.eval();
		}
		else
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os[ii], inputScalars.Nx[ii],
				inputScalars.Ny[ii], inputScalars.Nz[ii] * inputScalars.nRekos), 1, 2, 0));
		intIm.eval();
		dim_t dim0 = intIm.dims(0);
		dim_t dim1 = intIm.dims(1);
		dim_t dim2 = intIm.dims(2);
		if (DEBUG) {
			mexPrintBase("dim0 = %u\n", dim0);
			mexPrintBase("dim1 = %u\n", dim1);
			mexPrintBase("dim2 = %u\n", dim2);
			mexPrintBase("af::sum<float>(intIm) = %f\n", af::sum<float>(intIm));
			mexPrintBase("af::sum<float>(vec.im_os[ii]) = %f\n", af::sum<float>(vec.im_os[ii]));
			mexEval();
		}
		intIm = af::flat(intIm);
		af::sync();
#ifdef CUDA
		CUDA_TEXTURE_DESC texDesc;
		CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
		CUDA_RESOURCE_DESC resDesc;
		CUDA_RESOURCE_VIEW_DESC viewDesc;
		std::memset(&texDesc, 0, sizeof(texDesc));
		std::memset(&resDesc, 0, sizeof(resDesc));
		std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
		std::memset(&viewDesc, 0, sizeof(viewDesc));
		arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
		arr3DDesc.NumChannels = 1;
		arr3DDesc.Height = dim1;
		arr3DDesc.Width = dim0;
		arr3DDesc.Depth = dim2;
		status = cuArray3DCreate(&proj.integArrayXY, &arr3DDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to create integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray creation completed\n");
		CUDA_MEMCPY3D cpy3d;
		std::memset(&cpy3d, 0, sizeof(cpy3d));
		cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
		cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(intIm.device<CUdeviceptr>());
		cpy3d.srcPitch = dim0 * sizeof(float);
		cpy3d.srcHeight = dim1;
		cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		cpy3d.dstArray = proj.integArrayXY;
		cpy3d.WidthInBytes = dim0 * sizeof(float);
		cpy3d.Height = dim1;
		cpy3d.Depth = dim2;
		status = cuMemcpy3D(&cpy3d);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to copy integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray copy completed\n");
		resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
		resDesc.res.array.hArray = proj.integArrayXY;
		texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
		texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
		viewDesc.height = dim1;
		viewDesc.width = dim0;
		viewDesc.depth = dim2;
		viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
		status = cuTexObjectCreate(&proj.vec_opencl.d_image_os_int, &resDesc, &texDesc, &viewDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Integral image xz copy failed\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Texture creation completed\n");
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Synchronization failed\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Synchronization completed\n");
#elif defined(OPENCL)
		region = { static_cast<cl_ulong>(dim0), static_cast<cl_ulong>(dim1), static_cast<cl_ulong>(dim2) };
		proj.vec_opencl.d_image_os_int = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
		af::sync();
		cl_int status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*intIm.device<cl_mem>(), true), proj.vec_opencl.d_image_os_int, 0, proj.origin, region);
		if (status != 0) {
			getErrorString(status);
			mexPrint("Integral image xz copy failed\n");
			return -1;
		}
#endif
		af::sync();
		intIm.unlock();
		intIm = af::constant(0.f, inputScalars.Nx[ii] + 1, inputScalars.Nz[ii] + 1, inputScalars.Ny[ii]);
		if (inputScalars.meanFP) {
			im = af::reorder(af::moddims(vec.im_os[ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 0, 2, 1, 3);
			vec.meanFP(af::seq(inputScalars.Nx[ii], inputScalars.Nx[ii] + inputScalars.Ny[ii])) = af::flat(af::mean(af::mean(im, 0), 1));
			im -= af::tile(vec.meanFP(af::seq(inputScalars.Nx[ii], inputScalars.Nx[ii] + inputScalars.Ny[ii])), im.dims(0), im.dims(1), 1);
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
		}
		else
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os[ii], inputScalars.Nx[ii],
				inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 0, 2, 1, 3));
		dim0 = intIm.dims(0);
		dim1 = intIm.dims(1);
		dim2 = intIm.dims(2);
		intIm = af::flat(intIm);
		af::sync();
#ifdef CUDA
		arr3DDesc.Height = dim1;
		arr3DDesc.Width = dim0;
		arr3DDesc.Depth = dim2;
		status = cuArray3DCreate(&proj.FPArray, &arr3DDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to create integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray creation completed\n");
		cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(intIm.device<CUdeviceptr>());
		cpy3d.srcPitch = dim0 * sizeof(float);
		cpy3d.srcHeight = dim1;
		cpy3d.dstArray = proj.FPArray;
		cpy3d.WidthInBytes = dim0 * sizeof(float);
		cpy3d.Height = dim1;
		cpy3d.Depth = dim2;
		status = cuMemcpy3D(&cpy3d);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to copy integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray copy completed\n");
		resDesc.res.array.hArray = proj.FPArray;
		viewDesc.height = dim1;
		viewDesc.width = dim0;
		viewDesc.depth = dim2;
		status = cuTexObjectCreate(&proj.vec_opencl.d_image_os, &resDesc, &texDesc, &viewDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Integral image yz copy failed\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Texture creation completed\n");
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Synchronization failed\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Synchronization completed\n");
#elif defined(OPENCL)
		region = { static_cast<cl_ulong>(dim0), static_cast<cl_ulong>(dim1), static_cast<cl_ulong>(dim2) };
		proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
		af::sync();
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*intIm.device<cl_mem>(), true), proj.vec_opencl.d_image_os, 0, proj.origin, region);
		if (status != 0) {
			getErrorString(status);
			mexPrint("Integral image yz copy failed\n");
			return -1;
		}
		status = proj.CLCommandQueue[0].finish();
#endif
		intIm.unlock();
		af::deviceGC();
	}
	else {
		af::sync();
#ifdef CUDA
		if (inputScalars.useBuffers) {
			if (inputScalars.use_psf)
				proj.vec_opencl.d_im = transferAF(vec.im_os_blurred[ii]);
			else
				proj.vec_opencl.d_im = transferAF(vec.im_os[ii]);
		}
		else {
			CUdeviceptr* im;
			if (inputScalars.use_psf)
				im = vec.im_os_blurred[ii].device<CUdeviceptr>();
			else
				im = vec.im_os[ii].device<CUdeviceptr>();
			CUDA_TEXTURE_DESC texDesc;
			CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
			CUDA_RESOURCE_DESC resDesc;
			CUDA_RESOURCE_VIEW_DESC viewDesc;
			std::memset(&texDesc, 0, sizeof(texDesc));
			std::memset(&resDesc, 0, sizeof(resDesc));
			std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
			std::memset(&viewDesc, 0, sizeof(viewDesc));
			arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
			arr3DDesc.NumChannels = 1;
			arr3DDesc.Height = inputScalars.Ny[ii];
			arr3DDesc.Width = inputScalars.Nx[ii];
			arr3DDesc.Depth = inputScalars.Nz[ii];
			status = cuArray3DCreate(&proj.FPArray, &arr3DDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create image array\n");
				return -1;
			}
			status = cuCtxSynchronize();
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (DEBUG) {
				mexPrintBase("vec.im_os[ii].elements() = %u\n", vec.im_os[ii].elements());
				mexPrintBase("inputScalars.Nx[ii] * sizeof(float) = %u\n", inputScalars.Nx[ii] * sizeof(float));
				mexPrintBase("inputScalars.Ny[ii] = %u\n", inputScalars.Ny[ii]);
				mexPrintBase("inputScalars.Nz[ii] = %u\n", inputScalars.Nz[ii]);
				mexPrintBase("elements = %u\n", inputScalars.Nz[ii] * inputScalars.Nx[ii] * inputScalars.Ny[ii]);
				mexEval();
			}
			CUDA_MEMCPY3D cpy3d;
			std::memset(&cpy3d, 0, sizeof(cpy3d));
			cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
			cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(im);
			cpy3d.srcPitch = inputScalars.Nx[ii] * sizeof(float);
			cpy3d.srcHeight = inputScalars.Ny[ii];
			cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
			cpy3d.dstArray = proj.FPArray;
			cpy3d.WidthInBytes = inputScalars.Nx[ii] * sizeof(float);
			cpy3d.Height = inputScalars.Ny[ii];
			cpy3d.Depth = inputScalars.Nz[ii];
			status = cuMemcpy3D(&cpy3d);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to copy image array\n");
				status = cuArrayDestroy(proj.FPArray);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
				return -1;
			}
			status = cuCtxSynchronize();
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
			resDesc.res.array.hArray = proj.FPArray;
			if (inputScalars.FPType == 4) {
				texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
				texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
			}
			else {
				texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
			}
			viewDesc.height = inputScalars.Ny[ii];
			viewDesc.width = inputScalars.Nx[ii];
			viewDesc.depth = inputScalars.Nz[ii];
			viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
			status = cuTexObjectCreate(&proj.vec_opencl.d_image_os, &resDesc, &texDesc, &viewDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Integral image xz copy failed\n");
				status = cuArrayDestroy(proj.FPArray);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
				return -1;
			}
			status = cuCtxSynchronize();
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Synchronization failed\n");
				status = cuArrayDestroy(proj.FPArray);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
				return -1;
			}
			else if (DEBUG)
				mexPrint("Synchronization completed\n");
			if (inputScalars.use_psf)
				vec.im_os_blurred[ii].unlock();
			else
				vec.im_os[ii].unlock();
		}
#elif defined(OPENCL)
		if (inputScalars.useBuffers) {
			if (inputScalars.use_psf)
				proj.vec_opencl.d_im = cl::Buffer(*vec.im_os_blurred[ii].device<cl_mem>(), true);
			else
				proj.vec_opencl.d_im = cl::Buffer(*vec.im_os[ii].device<cl_mem>(), true);
			proj.CLCommandQueue[0].finish();
		}
		else {
			proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create input images\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Input image created\n");
			cl_mem* im;
			if (inputScalars.use_psf)
				im = vec.im_os_blurred[ii].device<cl_mem>();
			else
				im = vec.im_os[ii].device<cl_mem>();
			cl::Buffer d_im_os = cl::Buffer(*im, true);
			proj.CLCommandQueue[0].finish();
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(d_im_os, proj.vec_opencl.d_image_os, 0, proj.origin, region);
			proj.CLCommandQueue[0].finish();
			if (inputScalars.use_psf)
				vec.im_os_blurred[ii].unlock();
			else
				vec.im_os[ii].unlock();
			delete im;
			if (status != 0) {
				getErrorString(status);
				mexPrint("Image copy failed\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Input copy succeeded\n");
		}
#else
		if (inputScalars.use_psf)
			proj.vec_opencl.d_im_os = vec.im_os_blurred[ii].device<float>();
		else
			proj.vec_opencl.d_im_os = vec.im_os[ii].device<float>();
#endif
	}
	af::sync();
	return 0;
}

// Helper function to transfer AF arrays into device vectors, call the actual forward projection function and move memory management back to AF
inline int forwardProjectionAFOpenCL(AF_im_vectors& vec, scalarStruct& inputScalars, Weighting& w_vec, af::array& outputFP, uint32_t osa_iter,
	const std::vector<int64_t>& length, const af::array& g, uint64_t m_size, ProjectorClass& proj, const int ii = 0, const int64_t* pituus = nullptr) {
	int status = 0;
	if (inputScalars.use_psf)
		vec.im_os_blurred[ii] = computeConvolution(vec.im_os[ii], g, inputScalars, w_vec, inputScalars.nRekos2, ii);
	if (DEBUG) {
		mexPrintBase("outputFP.dims(0) = %d\n", outputFP.dims(0));
		mexPrintBase("outputFP.dims(1) = %d\n", outputFP.dims(1));
		mexEval();
	}
	proj.d_output = transferAF(outputFP);
	status = updateInputs(vec, inputScalars, proj, ii);
	if (status != 0) {
		return -1;
	}
	proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
#ifndef CPU
	if (inputScalars.meanFP && inputScalars.FPType == 5)
		proj.d_meanFP = transferAF(vec.meanFP);
	status = proj.forwardProjection(inputScalars, w_vec, osa_iter, length, m_size, ii);
	if (inputScalars.useBuffers)
		if (inputScalars.use_psf)
			vec.im_os_blurred[ii].unlock();
		else
			vec.im_os[ii].unlock();
#else
	status = proj.forwardProjection(inputScalars, w_vec, osa_iter, length, pituus, ii);
	if (inputScalars.use_psf)
		vec.im_os_blurred[ii].unlock();
	else
		vec.im_os[ii].unlock();
#endif
	outputFP.unlock();
	if (inputScalars.meanFP && inputScalars.FPType == 5)
		vec.meanFP.unlock();
	proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
	return status;
}

// Same as above, but for backprojection
inline int backwardProjectionAFOpenCL(AF_im_vectors& vec, scalarStruct& inputScalars, Weighting& w_vec, af::array& outputFP, uint32_t osa_iter,
	std::vector<int64_t>& length, uint64_t m_size, af::array& meanBP, const af::array& g, ProjectorClass& proj, const bool compSens = false, const int ii = 0, 
	const int64_t* pituus = nullptr, const bool FDK = false) {
	int status = 0;
	outputFP.eval();
	if (!FDK)
		initializeRHS(vec, inputScalars, ii);
	proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
	if (DEBUG) {
		mexPrintBase("ii = %u\n", ii);
		mexPrintBase("vec.rhs_os[ii].dims(0) = %u\n", vec.rhs_os[ii].dims(0));
		mexPrintBase("inputScalars.nRekos2 = %u\n", inputScalars.nRekos2);
		mexPrintBase("inputScalars.nRekos = %u\n", inputScalars.nRekos);
		mexPrintBase("outputFP = %f\n", af::sum<float>(outputFP));
		mexPrintBase("min(outputFP) = %f\n", af::min<float>(outputFP));
		mexPrintBase("max(outputFP) = %f\n", af::max<float>(outputFP));
		mexEval();
	}
	proj.d_output = transferAF(outputFP);
	if (DEBUG) {
		mexPrint("Transferring backprojection output\n");
	}
	status = transferRHS(vec.rhs_os[ii], proj);
	if (status != 0) {
		return -1;
	}
	 if (DEBUG)
		mexPrint("Backprojection output transfered\n");
#ifndef CPU
	if (inputScalars.meanBP && inputScalars.BPType == 5)
		proj.d_meanBP = transferAF(meanBP);
	status = proj.backwardProjection(inputScalars, w_vec, osa_iter, length, m_size, compSens, ii);
#else

	 status = proj.backwardProjection(inputScalars, w_vec, osa_iter, length, pituus, compSens, ii);
#endif
	vec.rhs_os[ii].unlock();
	outputFP.unlock();
	if (inputScalars.meanBP && inputScalars.BPType == 5)
		meanBP.unlock();
	if (status != 0) {
		return -1;
	}
#ifdef OPENCL
	if (inputScalars.atomic_64bit)
		vec.rhs_os[ii] = vec.rhs_os[ii].as(f32) / TH;
	else if (inputScalars.atomic_32bit)
		vec.rhs_os[ii] = vec.rhs_os[ii].as(f32) / TH32;
#endif
	if (inputScalars.use_psf)
		vec.rhs_os[ii] = computeConvolution(vec.rhs_os[ii], g, inputScalars, w_vec, inputScalars.nRekos2, ii);
	vec.rhs_os[ii].eval();
	outputFP.eval();
	return status;
}

// Computes custom median root prior, OpenCL, CUDA or CPU
inline int MRPAF(af::array& padd, af::array& grad, const scalarStruct& inputScalars, ProjectorClass& proj, const uint32_t medx, const uint32_t medy, const uint32_t medz) {

	int status = 0;
	if (DEBUG) {
		mexPrintBase("padd = %f\n", af::sum<float>(padd));
	}
	proj.d_W = transferAF(grad);
	proj.d_inputB = transferAF(padd);
#ifndef CPU
	uint64_t global_size[3] = { static_cast<uint64_t>(padd.dims(0)), static_cast<uint64_t>(padd.dims(1)), static_cast<uint64_t>(padd.dims(2)) };
	status = proj.computeMRP(inputScalars, global_size);
	grad.unlock();
	padd.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.computeMRP(inputScalars, medx, medy, medz);
	grad.unlock();
	padd.unlock();
#endif
	return status;
}

// NLM
inline int NLMAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, Weighting& w_vec, ProjectorClass& proj, const float beta, const int kk = 0) {

	int status = 0;
	proj.d_W = transferAF(grad);
	proj.d_gaussianNLM = transferAF(w_vec.gaussianNLM);
#ifndef CPU
	if (inputScalars.useImages) {
#ifdef CUDA
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
			if (DEBUG) {
				mexPrintBase("inputScalars.lDimStruct.NzPr[kk] = %u\n", inputScalars.lDimStruct.NzPr[kk]);
				mexPrintBase("inputScalars.lDimStruct.Nz[kk] = %u\n", inputScalars.lDimStruct.Nz[kk]);
				mexPrintBase("proj.region[2] = %u\n", proj.region[2]);
				mexPrintBase("im.elements() = %u\n", im.elements());
				mexPrintBase("grad.elements() = %u\n", grad.elements());
				mexEval();
			}
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			w_vec.gaussianNLM.unlock();
			mexPrint("Failed to copy NLM image\n");
			return -1;
		}
		if (inputScalars.largeDim)
			proj.region[2] = inputScalars.Nz[0];
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
	}
	status = proj.computeNLM(inputScalars, w_vec, beta, kk);
	grad.unlock();
	im.unlock();
	w_vec.gaussianNLM.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.computeNLM(inputScalars, w_vec, beta);
	grad.unlock();
	im.unlock();
	w_vec.gaussianNLM.unlock();
	grad *= beta;
	grad.eval();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// RDP
inline int RDPAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float gamma, ProjectorClass& proj, const float beta, const af::array& RDPref, const Weighting& w_vec,
	const bool RDPLargeNeighbor = false, const bool useRDPRef = false, const int kk = 0) {
	im.eval();
	int status = 0;
	proj.d_W = transferAF(grad);
#ifndef CPU
	if (inputScalars.useImages) {
#ifdef CUDA
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, Nz);
		if (RDPLargeNeighbor && useRDPRef) {
			CUdeviceptr* inputRef = RDPref.device<CUdeviceptr>();
			status = proj.transferTex(inputScalars, inputRef, true);
		}
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy RDP image\n");
			return -1;
		}
		if (RDPLargeNeighbor && useRDPRef) {
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*RDPref.device<cl_mem>(), true), proj.d_RDPrefI, 0, proj.origin, proj.region);
			if (status != 0) {
				getErrorString(status);
				im.unlock();
				grad.unlock();
				mexPrint("Failed to copy RDP image\n");
				return -1;
			}
		}
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
		if (RDPLargeNeighbor && useRDPRef)
			proj.d_RDPref = transferAF(RDPref);
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.computeRDP(inputScalars, gamma, w_vec, beta, kk, RDPLargeNeighbor, useRDPRef);
	grad.unlock();
	im.unlock();
	if (RDPLargeNeighbor && useRDPRef)
		RDPref.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.computeRDP(inputScalars, gamma, beta, RDPLargeNeighbor, useRDPRef);
	grad.unlock();
	im.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// GGMRF
inline int GGMRFAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, ProjectorClass& proj, 
	const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	proj.d_W = transferAF(grad);
#ifndef CPU
	if (inputScalars.useImages) {
#ifdef CUDA
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy GGMRF image\n");
			return -1;
		}
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.computeGGMRF(inputScalars, p, q, c, pqc, w_vec, beta, kk);
	grad.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.computeGGMRF(inputScalars, p, q, c, pqc, beta);
	grad.unlock();
	im.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// TV
inline int TVAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float sigma, const TVdata& data, ProjectorClass& proj, const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	int type = 0;
	proj.d_W = transferAF(grad);
	if (data.TV_use_anatomical)
		proj.d_refIm = transferAF(data.refIm);
	float C = 0.f;
	if (data.TVtype == 5) {
		C = data.eta;
		type = 3;
	}
	else if (data.TVtype == 2) {
		C = data.C;
		type = 2;
	}
	else if (data.TVtype == 1 && data.TV_use_anatomical)
		type = 1;
#ifndef CPU
	if (inputScalars.useImages) {
#ifdef CUDA
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy TV image\n");
			return -1;
		}
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	const float smooth = data.TVsmoothing;
	status = proj.TVGradient(inputScalars, sigma, smooth, w_vec, beta, kk, C, type);
	grad.unlock();
	im.unlock();
	if (data.TV_use_anatomical)
		data.refIm.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.TVGradient(inputScalars, data, sigma, data.TVsmoothing, beta, C, data.TVtype);
	grad.unlock();
	im.unlock();
	if (data.TV_use_anatomical)
		data.refIm.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// Hyperbolic prior
// No CPU support
inline int hyperAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float sigma, ProjectorClass& proj, const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	int type = 0;
#ifndef CPU
	proj.d_W = transferAF(grad);
	if (inputScalars.useImages) {
#ifdef CUDA
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy hyperbolic prior image\n");
			return -1;
		}
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.hyperGradient(inputScalars, sigma, w_vec, beta, kk);
	grad.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
#else

#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}


#ifndef CPU
// Proximal TV
inline int proxQAF(std::vector<af::array>& q, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = q[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", q[0].elements());
		mexEval();
	}
	proj.d_qX = transferAF(q[0]);
	status = proj.ProxHelperQ(alpha, globalQ);
	q[0].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV
inline int proxTVQAF(std::vector<af::array>& q, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = q[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", q[0].elements());
		mexEval();
	}
	proj.d_qX = transferAF(q[0]);
	proj.d_qY = transferAF(q[1]);
	proj.d_qZ = transferAF(q[2]);
	status = proj.ProxTVHelperQ(alpha, globalQ);
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV
inline int proxTGVQAF(std::vector<af::array>& r, const scalarStruct& inputScalars, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = r[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", r[0].elements());
		mexEval();
	}
	proj.d_rX = transferAF(r[0]);
	proj.d_rY = transferAF(r[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(r[2]);
		proj.d_rXY = transferAF(r[3]);
		proj.d_rXZ = transferAF(r[4]);
		proj.d_rYZ = transferAF(r[5]);
	}
	else {
		proj.d_rXY = transferAF(r[2]);
	}
	status = proj.ProxTGVHelperQ(inputScalars, alpha, globalQ);
	r[0].unlock();
	r[1].unlock();
	r[2].unlock();
	if (!inputScalars.TGV2D) {
		r[3].unlock();
		r[4].unlock();
		r[5].unlock();
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV divergence
inline int proxTVDivAF(const std::vector<af::array>& grad, af::array& input, const scalarStruct& inputScalars, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("input.dims(0) = %u\n", input.dims(0));
		mexPrintBase("grad[0].dims(0) = %u\n", grad[0].dims(0));
		mexPrintBase("grad[0].dims(1) = %u\n", grad[0].dims(1));
		mexPrintBase("grad[0].dims(2) = %u\n", grad[0].dims(2));
		mexEval();
	}
	proj.d_qX = transferAF(grad[0]);
	proj.d_qY = transferAF(grad[1]);
	proj.d_qZ = transferAF(grad[2]);
	proj.vec_opencl.d_rhs_os[0] = transferAF(input);
	status = proj.ProxTVDiv(inputScalars);
	grad[0].unlock();
	grad[1].unlock();
	grad[2].unlock();
	input.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV gradient
inline int proxTVGradAF(const af::array& im, std::vector<af::array>& grad, const scalarStruct& inputScalars, const float sigma2, const std::vector<af::array>& v, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("output.dims(0) = %u\n", grad[0].dims(0));
		mexPrintBase("im.dims(0) = %u\n", im.dims(0));
		mexPrintBase("im.dims(1) = %u\n", im.dims(1));
		mexPrintBase("im.dims(2) = %u\n", im.dims(2));
		mexEval();
	}
	const size_t vSize = v.size();
	proj.d_qX = transferAF(grad[0]);
	proj.d_qY = transferAF(grad[1]);
	proj.d_qZ = transferAF(grad[2]);
	proj.d_inputB = transferAF(im);
	if (v.size() > 0) {
		if (DEBUG) {
			mexPrintBase("v0.dims(0) = %u\n", v[0].dims(0));
			mexEval();
		}
		if (DEBUG) {
			mexPrintBase("v1.dims(0) = %u\n", v[1].dims(0));
			mexEval();
		}
		proj.d_vX = transferAF(v[0]);
		proj.d_vY = transferAF(v[1]);
		if (!inputScalars.TGV2D)
			proj.d_vZ = transferAF(v[2]);
	}
	status = proj.ProxTVGrad(inputScalars, sigma2, vSize);
	grad[0].unlock();
	grad[1].unlock();
	grad[2].unlock();
	im.unlock();
	if (v.size() > 0) {
		v[0].unlock();
		v[1].unlock();
		if (!inputScalars.TGV2D) {
			v[2].unlock();
		}
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV symmetric derivative
inline int proxTGVSymmDerivAF(const std::vector<af::array>& v, std::vector<af::array>& q, const scalarStruct& inputScalars, const float sigma2, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("input.dims(0) = %u\n", v[0].dims(0));
		if (!inputScalars.TGV2D)
			mexPrintBase("input2.dims(0) = %u\n", v[2].dims(0));
		mexPrintBase("im.dims(0) = %u\n", q[0].dims(0));
		mexPrintBase("im.dims(1) = %u\n", q[0].dims(1));
		mexPrintBase("im.dims(2) = %u\n", q[0].dims(2));
		mexPrintBase("q1.dims(0) = %u\n", q[1].dims(0));
		mexPrintBase("q2.dims(0) = %u\n", q[2].dims(0));
		if (!inputScalars.TGV2D) {
			mexPrintBase("q3.dims(0) = %u\n", q[3].dims(0));
			mexPrintBase("q5.dims(0) = %u\n", q[5].dims(0));
		}
		mexPrintBase("v.size() = %u\n", v.size());
		mexPrintBase("q.size() = %u\n", q.size());
		mexEval();
	}
	proj.d_rX = transferAF(q[0]);
	proj.d_rY = transferAF(q[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(q[2]);
		proj.d_rXY = transferAF(q[3]);
		proj.d_rXZ = transferAF(q[4]);
		proj.d_rYZ = transferAF(q[5]);
	}
	else
		proj.d_rXY = transferAF(q[2]);
	proj.d_vX = transferAF(v[0]);
	proj.d_vY = transferAF(v[1]);
	if (!inputScalars.TGV2D)
		proj.d_vZ = transferAF(v[2]);
	status = proj.ProxTGVSymmDeriv(inputScalars, sigma2);
	v[0].unlock();
	v[1].unlock();
	if (!inputScalars.TGV2D)
		v[2].unlock();
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (!inputScalars.TGV2D) {
		q[3].unlock();
		q[4].unlock();
		q[5].unlock();
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV divergence
inline int proxTGVDivAF(const std::vector<af::array>& q, std::vector<af::array>& v, std::vector<af::array>& p, const scalarStruct& inputScalars,
	const float theta, const float tau, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("v.dims(0) = %u\n", v[0].dims(0));
		mexPrintBase("q2.dims(0) = %u\n", q[0].dims(0));
		mexPrintBase("q2.dims(1) = %u\n", q[0].dims(1));
		mexPrintBase("q2.dims(2) = %u\n", q[0].dims(2));
		mexEval();
	}
	proj.d_rX = transferAF(q[0]);
	proj.d_rY = transferAF(q[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(q[2]);
		proj.d_rXY = transferAF(q[3]);
		proj.d_rXZ = transferAF(q[4]);
		proj.d_rYZ = transferAF(q[5]);
	}
	else
		proj.d_rXY = transferAF(q[2]);
	proj.d_vX = transferAF(v[0]);
	proj.d_vY = transferAF(v[1]);
	if (!inputScalars.TGV2D)
		proj.d_vZ = transferAF(v[2]);
	proj.d_qX = transferAF(p[0]);
	proj.d_qY = transferAF(p[1]);
	proj.d_qZ = transferAF(p[2]);
	status = proj.ProxTGVDiv(inputScalars, theta, tau);
	v[0].unlock();
	v[1].unlock();
	if (!inputScalars.TGV2D)
		v[2].unlock();
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (!inputScalars.TGV2D) {
		q[3].unlock();
		q[4].unlock();
		q[5].unlock();
	}
	p[0].unlock();
	p[1].unlock();
	p[2].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Computes elementwise multiplication or division, depending on mult-variable
inline int elementWiseAF(const af::array& vector, af::array& input, const bool mult, ProjectorClass& proj, const bool D2 = false) {
	int status = 0;
	uint64_t gSize[3] = { static_cast<uint64_t>(input.dims(0)), static_cast<uint64_t>(input.dims(1)), static_cast<uint64_t>(input.dims(2)) };
	if (DEBUG) {
		mexPrintBase("vector.dims[0] = %u\n", vector.dims(0));
		mexPrintBase("vector.dims[1] = %u\n", vector.dims(1));
		mexPrintBase("vector.dims[2] = %u\n", vector.dims(2));
		mexPrintBase("input.dims[0] = %u\n", input.dims(0));
		mexPrintBase("input.dims[1] = %u\n", input.dims(1));
		mexPrintBase("input.dims[2] = %u\n", input.dims(2));
		mexEval();
	}
	proj.d_vector = transferAF(vector);
	proj.d_input = transferAF(input);
	status = proj.elementWiseComp(mult, gSize, D2);
	vector.unlock();
	input.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Poisson-estimate update kernel call for PKMA, MBSREM and BSREM
inline int poissonUpdateAF(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, const float lambda, const float epps, const float alpha, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	proj.d_im = transferAF(im);
	proj.d_rhs = transferAF(rhs);
	status = proj.PoissonUpdate(inputScalars, lambda, epps, alpha, ii);
	rhs.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Same as above, but for PDHG
inline int PDHGUpdateAF(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, AF_im_vectors& vec, const float epps, const float theta, const float tau, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	proj.d_im = transferAF(im);
	proj.d_rhs = transferAF(rhs);
	proj.d_U = transferAF(vec.uCP[ii]);
	status = proj.PDHGUpdate(inputScalars, epps, theta, tau, ii);
	rhs.unlock();
	im.unlock();
	vec.uCP[ii].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

inline int rotateCustomAF(af::array& imrot, const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angle, const int ii = 0) {
	int status = 0;
	if (!inputScalars.useBuffers) {
#if defined(CUDA)
		CUdeviceptr* input = im.device<CUdeviceptr>();
		status = proj.transferTex(inputScalars, input, false, inputScalars.Nz[0]);
#elif defined(OPENCL)
        status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
        if (status != 0) {
			getErrorString(status);
			im.unlock();
			mexPrint("Failed to copy rotation image\n");
			return -1;
		}
#endif
    } else {
        proj.d_im = transferAF(im);
    }
	proj.d_rhs = transferAF(imrot);
	const float cosa = std::cos(-angle);
	const float sina = std::sin(-angle);
	status = proj.rotateCustom(inputScalars, cosa, sina);
	imrot.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}
#endif

// Various batch functions
inline af::array batchMinus(const af::array& lhs, const af::array& rhs) {
	return lhs - rhs;
}

inline af::array batchPlus(const af::array& lhs, const af::array& rhs) {
	return lhs + rhs;
}

inline af::array batchMul(const af::array& lhs, const af::array& rhs) {
	return lhs * rhs;
}

inline af::array batchDiv(const af::array& lhs, const af::array& rhs) {
	return lhs / rhs;
}

inline af::array batchNotEqual(const af::array& lhs, const af::array& rhs) {
	return lhs != rhs;
}

// Compute the epsilon value for the MBSREM/MRAMLA
inline float MBSREM_epsilon(const af::array& Sino, const af::array& D, const float epps = 1e-5f, const uint32_t randoms_correction = false, const af::array& rand = af::constant(0.1, 1, 1),
	const bool TOF = false, const int64_t nBins = 1, const bool CT = false) {
	float eps;
	if (CT) {
		af::array hk_summa = -af::exp(-Sino) / Sino - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (randoms_correction == 1u) {
			Iind = (Sino > 0.f & rand == 0.f);
			if (af::sum<float>(Iind) == 0.f)
				return 1e8f;
			P_Sino = Sino(Iind);
			apu = D + rand;
			apu = af::sum(-af::exp(-apu) / Sino - apu);
		}
		else {
			Iind = (Sino > 0.f);
			P_Sino = Sino(Iind);
			apu = af::sum(-af::exp(-D) / Sino - D);
		}
		if (randoms_correction == 1u)
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		else
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		af::array epsilon = (af::min)(P_Sino, af::log(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
	else {
		af::array hk_summa = Sino * af::log(Sino) - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (TOF && randoms_correction) {
			af::array rInd = rand == 0.f;
			P_Sino = Sino(Sino > 0.f & af::tile(rInd, nBins));
			apu = D + af::tile(rand, nBins);
			apu = af::sum(Sino * af::log(apu) - apu);
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Sino > 0.f & af::tile(rInd, nBins)), batchMinus);
		}
		else {
			if (randoms_correction == 1u) {
				Iind = (Sino > 0.f & rand == 0.f);
				if (af::sum<float>(Iind) == 0.f)
					return 1e8f;
				P_Sino = Sino(Iind);
				apu = D + rand;
				apu = af::sum(Sino * af::log(apu) - apu);
			}
			else {
				Iind = (Sino > 0.f);
				P_Sino = Sino(Iind);
				apu = af::sum(Sino * af::log(D) - D);
			}
			if (randoms_correction == 1u)
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
			else
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		}
		af::array epsilon = (af::min)(P_Sino, af::exp(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
	eps = eps <= 0.f ? epps : eps;
	return eps;
}

// SPECT forward projection (projector type 6)
inline void forwardProjectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const int ii = 0, const float* atten = nullptr) {
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT forward projection");
	int64_t u1 = uu;
	const af::array apuArr = af::moddims(vec.im_os[ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	if (DEBUG)
		mexPrint("step 1");
	for (int kk = 0; kk < length; kk++) {
		af::array attenuationImage;
		af::array kuvaRot;
#ifndef CPU
		kuvaRot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
        /* apuFP = imtranslate(apuFP, [-P0(1); P0(2); 0]', 'bilinear', 'FillValues', 0); % Translate FP */
		rotateCustomAF(kuvaRot, apuArr, inputScalars, proj, (180-w_vec.angles[u1])*M_PI/180., ii);
#else
        /* apuFP = imtranslate(apuFP, [-P0(1); P0(2); 0]', 'bilinear', 'FillValues', 0); % Translate FP */
		kuvaRot = af::rotate(apuArr, (180-w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
		kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
		if (DEBUG)
			mexPrint("step 3");
		if (inputScalars.attenuation_correction && (atten != nullptr)) {
			attenuationImage = af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], atten);
			if (DEBUG)
				mexPrint("step 4");
#ifndef CPU
			af::array attn = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
            /* attenuationImage = imtranslate(attenuationImage,  [-P0(1); P0(2); 0]', 'bilinear', 'FillValues', 0); % Translate attenuation image */
			rotateCustomAF(attn, attenuationImage, inputScalars, proj, (180-w_vec.angles[u1])*M_PI/180., ii);
			attenuationImage = attn.copy();
#else
            /* attenuationImage = imtranslate(attenuationImage,  [-P0(1); P0(2); 0]', 'bilinear', 'FillValues', 0); % Translate attenuation image */
			attenuationImage = af::rotate(attenuationImage, (180-w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
			attenuationImage = af::accum(attenuationImage, 0);
			attenuationImage = af::exp(-w_vec.dPitchX * attenuationImage);
			if (DEBUG)
				mexPrint("step 5");
			//kuvaRot = kuvaRot * attenuationImage;
			attenuationImage = af::reorder(attenuationImage, 2, 1, 0);
			if (DEBUG) {
				mexPrintBase("af::sum(attenuationImage) = %f\n", af::sum<float>(attenuationImage));
				mexPrintBase("attenuationImageFP.dims(0) = %d\n", attenuationImage.dims(0));
				mexEval();
			}
		}
		kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
		kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
		if (inputScalars.attenuation_correction && (atten != nullptr)) {
			attenuationImage = af::convolve2(attenuationImage, w_vec.gFilter(af::span, af::span, af::span, u1));
			attenuationImage = af::reorder(attenuationImage, 2, 1, 0);
			kuvaRot *= attenuationImage;
		}
		kuvaRot = kuvaRot(af::seq(w_vec.distInt[u1], af::end), af::span, af::span);
		kuvaRot = af::sum(kuvaRot, 0);
		kuvaRot = af::reorder(kuvaRot, 1, 2, 0);
		fProj(af::span, af::span, kk) += kuvaRot.copy();
		u1++;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("SPECT forward projection completed in %f seconds\n", tDiff);
		mexEval();
	}
}

// SPECT backprojection (projector type 6)
inline void backprojectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, const uint32_t osa_iter = 0, const uint32_t iter = 0,
	const uint8_t compute_norm_matrix = 0, const uint32_t iter0 = 0, const int ii = 0, const float* atten = nullptr) {
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT backprojection");
	fProj = af::moddims(fProj, inputScalars.nRowsD, inputScalars.nColsD, length);
	af::array apuBP2 = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], length);
	int64_t u1 = uu;
	if (DEBUG) {
		mexPrintBase("ii = %d\n", ii);
		mexPrintBase("length = %d\n", length);
		mexEval();
	}
	for (int kk = 0; kk < length; kk++) {
		af::array apuBP = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
		af::array kuvaRot = fProj(af::span, af::span, kk).copy();
		kuvaRot = af::reorder(kuvaRot, 1, 0, 2);
		kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
		kuvaRot = kuvaRot(af::span, af::span, af::seq(w_vec.distInt[u1], af::end));
		kuvaRot = reorder(kuvaRot, 2, 1, 0);
		af::eval(kuvaRot);
		apuBP(af::seq(w_vec.distInt[u1], af::end), af::span, af::span) = kuvaRot.copy();
#ifndef CPU
		kuvaRot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
		rotateCustomAF(kuvaRot, apuBP, inputScalars, proj, (180+w_vec.angles[u1])*M_PI/180., ii);
		apuBP = kuvaRot.copy();
		//apuBP = af::rotate(apuBP, w_vec.angles[u1]);
#else
		apuBP = af::rotate(apuBP, (180+w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
		if (DEBUG) {
			mexPrintBase("w_vec.angles[u1] = %f\n", w_vec.angles[u1]);
			mexEval();
		}
		if (inputScalars.attenuation_correction && (atten != nullptr)) {
			af::array attenuationImage = af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], atten);
#ifndef CPU
			kuvaRot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
			rotateCustomAF(kuvaRot, attenuationImage, inputScalars, proj, (180+w_vec.angles[u1])*M_PI/180., ii);
			attenuationImage = kuvaRot.copy();
			//attenuationImage = af::rotate(attenuationImage, w_vec.angles[u1]);
#else
			attenuationImage = af::rotate(attenuationImage, (180+w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
			attenuationImage = af::accum(attenuationImage, 0);
			attenuationImage = af::exp(-w_vec.dPitchX * attenuationImage);
			apuBP *= attenuationImage;
			af::eval(apuBP);
			if (DEBUG) {
				mexPrintBase("af::sum(attenuationImage) = %f\n", af::sum<float>(attenuationImage));
				mexPrintBase("attenuationImage.dims(0) = %d\n", attenuationImage.dims(0));
				mexPrintBase("w_vec.dPitchX = %f\n", w_vec.dPitchX);
				mexEval();
			}
		}
		af::eval(apuBP);
		apuBP2(af::span, kk) = af::flat(apuBP).copy();
		u1++;
	}
	af::sync();
	if (DEBUG) {
		mexPrintBase("u1 = %d\n", u1);
		mexPrintBase("af::sum(apuBP2) = %f\n", af::sum<float>(apuBP2));
		mexPrintBase("vec.rhs_os[ii].dims(0) = %d\n", vec.rhs_os[ii].dims(0));
		mexEval();
	}
	vec.rhs_os[ii] = af::sum(apuBP2, 1);
	vec.rhs_os[ii](vec.rhs_os[ii] < inputScalars.epps && vec.rhs_os[ii] >= 0.f) = inputScalars.epps;
	if ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Computing sensitivity image");
		apuBP2 = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], length);
		u1 = uu;
		for (int kk = 0; kk < length; kk++) {
			af::array apuSumm = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
			af::array kuvaRot = af::constant(1.f, inputScalars.nColsD, inputScalars.nRowsD);
			kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
			kuvaRot = kuvaRot(af::span, af::span, af::seq(w_vec.distInt[u1], af::end));
			kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
			apuSumm(af::seq(w_vec.distInt[u1], af::end), af::span, af::span) = kuvaRot.copy();
#ifndef CPU
			kuvaRot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
			rotateCustomAF(kuvaRot, apuSumm, inputScalars, proj, (180+w_vec.angles[u1])*M_PI/180., ii);
			apuSumm = kuvaRot.copy();
			//apuSumm = af::rotate(apuSumm, w_vec.angles[u1]);
#else
			apuSumm = af::rotate(apuSumm, (180+w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
			if (inputScalars.attenuation_correction && (atten != nullptr)) {
				af::array attenuationImage = af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], atten);
#ifndef CPU
				kuvaRot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
				rotateCustomAF(kuvaRot, attenuationImage, inputScalars, proj, (180+w_vec.angles[u1])*M_PI/180., ii);
				attenuationImage = kuvaRot.copy();
				//attenuationImage = af::rotate(attenuationImage, w_vec.angles[u1]);
#else
				attenuationImage = af::rotate(attenuationImage, (180+w_vec.angles[u1])*M_PI/180., true, AF_INTERP_BILINEAR);
#endif
				attenuationImage = af::accum(attenuationImage, 0);
				attenuationImage = af::exp(-w_vec.dPitchX * attenuationImage);
				apuSumm = apuSumm * attenuationImage;
				af::eval(apuSumm);
			}
			apuBP2(af::span, kk) = af::flat(apuSumm);
			u1++;
		}
		if (DEBUG) {
			mexPrintBase("af::sum(apuBP2, 1) = %f\n", af::sum(apuBP2, 1));
			mexEval();
		}
		if (compute_norm_matrix == 2) {
			vec.Summ[ii][osa_iter] = af::sum(apuBP2, 1);
			vec.Summ[ii][osa_iter](vec.Summ[ii][osa_iter] < inputScalars.epps) = 1.f;
		}
		else {
			vec.Summ[ii][0] = af::sum(apuBP2, 1);
			vec.Summ[ii][0](vec.Summ[ii][0] < inputScalars.epps) = 1.f;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Sensitivity image computed");
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("SPECT backprojection complete");
}


// 1D filtering
inline int filtering(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	af::array temp = af::fft(input, dimmi);
	temp.eval();
	if (DEBUG) {
		mexPrintBase("temp = %f\n", af::sum<float>(temp));
		mexEval();
	}
#ifndef CPU
	status = elementWiseAF(filter, temp, true, proj);
	if (status != 0)
		return -1;
#else
	temp *= af::tile(filter, 1, temp.dims(1), temp.dims(2));
#endif
	af::sync();
	if (DEBUG) {
		mexPrintBase("temp = %f\n", af::sum<float>(temp));
		mexEval();
	}
	af::ifftInPlace(temp);
	temp.eval();
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::span, af::span)));
	return 0;
}

// 2D filtering
inline int filtering2D(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	const dim_t inDim0 = input.dims(0);
	if (DEBUG) {
		mexPrintBase("filter.dims(0) = %d\n", filter.dims(0));
		mexPrintBase("filter.dims(1) = %d\n", filter.dims(1));
		mexPrintBase("input.dims(0) = %d\n", input.dims(0));
		mexPrintBase("input.dims(1) = %d\n", input.dims(1));
		mexPrintBase("input.dims(2) = %d\n", input.dims(2));
		mexPrintBase("dimmi = %d\n", dimmi);
		mexEval();
	}
	af::array temp = af::fft2(input, dimmi, dimmi);
#ifndef CPU
	status = elementWiseAF(filter, temp, true, proj, true);
	if (status != 0)
		return -1;
#else
	temp *= filter;
#endif
	ifft2InPlace(temp);
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::seq(0, input.dims(0) - 1), af::span)));
	input.eval();
	af::deviceGC();
	return 0;
}

// 1D filtering with division
inline int filteringInverse(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	af::array temp = af::fft(input, dimmi);
	temp.eval();
#ifndef CPU
	status = elementWiseAF(filter, temp, false, proj);
	if (status != 0)
		return -1;
#else
	temp /= af::tile(filter, 1, temp.dims(1), temp.dims(2));
	temp.eval();
#endif
	af::sync();
	af::ifftInPlace(temp);
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::span, af::span)));
	return 0;
}

// Image-based preconditioner 2
inline af::array precondIm2(const af::array& im, const af::array& D) {
	return im / D;
}

// Image-based preconditioner 3
inline af::array precondIm3(const af::array& im, const af::array& D, const af::array& ref) {
	return (af::max)(im, (af::max)(VAL, ref)) / D;
}

// Compute the gradient of an array
inline void computeGradient(const af::array& im, const scalarStruct& inputScalars, af::array& f, af::array& g, af::array& h, const int type = 0) {
	if (DEBUG) {
		mexPrintBase("im.dims(0) = %d\n", im.dims(0));
		mexPrintBase("im.dims(1) = %d\n", im.dims(1));
		mexPrintBase("im.dims(2) = %d\n", im.dims(2));
		mexEval();
	}
	// 1st order forward differences
	if (type == 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting forward difference gradient");
		f(af::seq(0, af::end - 1LL), af::span, af::span) = af::diff1(im);
		f(af::end, af::span, af::span) = -1.f * im(af::end, af::span, af::span);
		g(af::span, af::seq(0, af::end - 1LL), af::span) = af::diff1(im, 1);
		g(af::span, af::end, af::span) = -1.f * im(af::span, af::end, af::span);
		h(af::span, af::span, af::seq(0, af::end - 1LL)) = af::diff1(im, 2);
		h(af::span, af::span, af::end) = -1.f * im(af::span, af::span, af::end);
	}
	// 1st order backward differences
	else if (type == 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting backward difference gradient");
		f(af::seq(1, af::end), af::span, af::span) = -af::diff1(im);
		f(0, af::span, af::span) = -1.f * im(0, af::span, af::span);
		g(af::span, af::seq(1, af::end), af::span) = -af::diff1(im, 1);
		g(af::span, 0, af::span) = -1.f * im(af::span, 0, af::span);
		h(af::span, af::span, af::seq(1, af::end)) = -af::diff1(im, 2);
		h(af::span, af::span, 0) = -1.f * im(af::span, af::span, 0);
	}
	// 1st order central differences
	else {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting central difference gradient");
		f = (af::shift(im, -1) - af::shift(im, 1)) * .5f;
		f(0, af::span, af::span) = im(1, af::span, af::span) - im(0, af::span, af::span);
		f(af::end, af::span, af::span) = im(af::end, af::span, af::span) - im(af::end - 1, af::span, af::span);
		g = (af::shift(im, 0, -1) - af::shift(im, 0, 1)) * .5f;
		g(af::span, 0, af::span) = im(af::span, 1, af::span) - im(af::span, 0, af::span);
		g(af::span, af::end, af::span) = im(af::span, af::end, af::span) - im(af::span, af::end - 1, af::span);
		h = (af::shift(im, 0, 0, -1) - af::shift(im, 0, 0, 1)) * .5f;
		h(af::span, af::span, 0) = im(af::span, af::span, 1) - im(af::span, af::span, 0);
		h(af::span, af::span, af::end) = im(af::span, af::span, af::end) - im(af::span, af::span, af::end - 1);
	}
	f = af::flat(f);
	g = af::flat(g);
	h = af::flat(h);
	f.eval();
	g.eval();
	h.eval();
	if (inputScalars.verbose >= 3)
		mexPrint("Gradient computed");
}

// Second order gradient
inline void computeSecondOrderGradient(af::array& f, af::array& g, af::array& h, const int type = 0) {
	// Second order derivatives
	if (type == 0) {
		f = af::shift(f, 1);
		g = af::shift(g, 0, 1);
		h = af::shift(h, 0, 0, 1);
	}
	else if (type == 1) {
		f = af::shift(f, -1);
		g = af::shift(g, 0, -1);
		h = af::shift(h, 0, 0, -1);
	}
	else {
		f = (af::shift(f, -1) - af::shift(f, 1)) * .5f;
		g = (af::shift(g, 0, -1) - af::shift(g, 0, 1)) * .5f;
		h = (af::shift(h, 0, 0, -1) - af::shift(h, 0, 0, 1)) * .5f;
	}
}

// Image-based preconditioner 4
inline void gradientPreconditioner(const scalarStruct& inputScalars, Weighting& w_vec, const af::array& input, const int ii = 0) {
	af::array f = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	af::array g = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	af::array h = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	computeGradient(input, inputScalars, f, g, h, w_vec.derivType);
	if (DEBUG) {
		mexPrintBase("g.dims(0) = %d\n", g.dims(0));
		mexEval();
	}
	f = (af::max)(VAL, af::sqrt(f * f + g * g + h * h) / af::mean<float>(af::flat(input)));
	f = af::mean<float>(af::flat(f)) / f;
	f.eval();
	if (DEBUG) {
		mexPrintBase("f.dims(0) = %d\n", f.dims(0));
		mexEval();
	}
	w_vec.gradF[ii] = (af::min)(w_vec.gradV2, (af::max)(w_vec.gradV1, f));
}

// Compute the image-based preconditioning
inline int applyImagePreconditioning(Weighting& w_vec, const scalarStruct& inputScalars, af::array& input, const af::array& im, ProjectorClass& proj, const int kk = 0, const int ii = 0) {
	int status = 0;
	if (w_vec.precondTypeIm[4] && kk >= w_vec.gradInitIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying gradient-based preconditioner, type 4");
		if (kk <= w_vec.gradFinalIter)
			gradientPreconditioner(inputScalars, w_vec, af::moddims(im, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]), ii);
		input *= w_vec.gradF[ii];
	}
	if (w_vec.precondTypeIm[3]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying momentum-like preconditioner, type 3");
		input *= w_vec.alphaPrecond[kk];
	}
	if (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2]) {
		if (w_vec.precondTypeIm[0]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying diagonal normalization preconditioner , type 0");
			input /= w_vec.D[ii];
		}
		else if (w_vec.precondTypeIm[1]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying EM preconditioner, type 1");
			input *= precondIm2(im, w_vec.D[ii]);
		}
		else if (w_vec.precondTypeIm[2]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying IEM preconditioner, type 2");
			input *= precondIm3(im, w_vec.D[ii], w_vec.preRef[ii]);
		}
	}
	if (w_vec.precondTypeIm[6]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying curvature preconditioner , type 6");
		input *= w_vec.dP[ii];
	}
	if (w_vec.precondTypeIm[5] && kk <= w_vec.filterIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying filtering-based preconditioner, type 5");
		af::deviceGC();
		input = af::moddims(input, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
		status = filtering2D(w_vec.filterIm, input, proj, inputScalars.Nf);
		if (status != 0)
			return -1;
	}
	input.eval();
	af::deviceGC();
	if (inputScalars.verbose >= 3 && (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2] || w_vec.precondTypeIm[3] ||
		(w_vec.precondTypeIm[4] && kk >= w_vec.gradInitIter) || w_vec.precondTypeIm[5] || w_vec.precondTypeIm[6]))
		mexPrint("Image-based preconditioning applied");
	return 0;
}

// Compute measurement-based preconditioning
inline int applyMeasPreconditioning(const Weighting& w_vec, const scalarStruct& inputScalars, af::array& input, ProjectorClass& proj, const uint32_t subIter = 0) {
	int status = 0;
	if (w_vec.precondTypeMeas[0] || w_vec.precondTypeMeas[1]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying measurement-based preconditioning");
		if (w_vec.precondTypeMeas[1]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying filtering-based preconditioner, type 1");
			if (DEBUG) {
				mexPrintBase("input.elements() = %d\n", input.elements());
				mexPrintBase("w_vec.filter.elements() = %d\n", w_vec.filter.elements());
				mexPrintBase("w_vec.filter = %f\n", af::sum<float>(w_vec.filter));
				mexPrintBase("inputScalars.nRowsD = %d\n", inputScalars.nRowsD);
				mexPrintBase("inputScalars.nColsD = %d\n", inputScalars.nColsD);
				mexPrintBase("input.elements() / (inputScalars.nRowsD * inputScalars.nColsD) = %d\n", input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
				mexEval();
			}
			if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType == 5 || inputScalars.subsetType == 4)) {
				if (inputScalars.subsetType == 4)
					input = af::moddims(input, inputScalars.nRowsD, input.elements() / inputScalars.nRowsD);
				else
					input = af::moddims(input, inputScalars.nColsD, input.elements() / inputScalars.nColsD);
			}
			else
				input = af::moddims(input, inputScalars.nRowsD, inputScalars.nColsD, input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
			input.eval();
			status = filtering(w_vec.filter, input, proj, inputScalars.Nf);
			if (status != 0)
				return -1;
		}
		if (w_vec.precondTypeMeas[0]) {
			if (DEBUG) {
				mexPrintBase("w_vec.M[subIter].dims(0) = %d\n", w_vec.M[subIter].dims(0));
				mexPrintBase("input.dims(0) = %d\n", input.dims(0));
				mexEval();
			}
			if (inputScalars.verbose >= 3)
				mexPrint("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
			input /= w_vec.M[subIter];
		}
		input.eval();
		af::deviceGC();
		if (inputScalars.verbose >= 3)
			mexPrint("Measurement-based preconditioning applied");
	}
	return 0;
}

inline void computeIntegralImage(const scalarStruct& inputScalars, const Weighting& w_vec, const int64_t length, af::array& outputFP, af::array& meanBP) {
	if (inputScalars.BPType == 5) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Computing integral image for backprojection");
		if (DEBUG) {
			mexPrintBase("outputFP.dims(0) = %d\n", outputFP.dims(0));
			mexPrintBase("inputScalars.nRowsD = %d\n", inputScalars.nRowsD);
			mexPrintBase("inputScalars.nColsD = %d\n", inputScalars.nColsD);
			mexPrintBase("length = %d\n", length);
			mexEval();
		}
		af::sync();
		af::deviceGC();
		outputFP = af::moddims(outputFP, inputScalars.nRowsD, inputScalars.nColsD, length);
		if (inputScalars.meanBP) {
			meanBP = af::mean(af::mean(outputFP, 0), 1);
			outputFP -= af::tile(meanBP, inputScalars.nRowsD, inputScalars.nColsD, 1);
			outputFP.eval();
		}
		outputFP = af::sat(outputFP);
		outputFP = af::join(0, af::constant(0.f, 1, outputFP.dims(1), outputFP.dims(2)), outputFP);
		outputFP = af::flat(af::join(1, af::constant(0.f, outputFP.dims(0), 1, outputFP.dims(2)), outputFP));
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Integral images computed");
		af::sync();
		af::deviceGC();
	}
}

// Compute deconvolution (deblur/sharpening)
inline void deblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Starting deconvolution");
	af::array jelppi = padding(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	af::array apu = padding(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	for (uint32_t kk = 0U; kk < inputScalars.deblur_iterations; kk++) {
		af::array apu2 = af::convolve3(jelppi, g) + inputScalars.epps;
		apu2 = apu2(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::seq(inputScalars.g_dim_z + 1, inputScalars.Nz[ii] + inputScalars.g_dim_z));
		apu2 = padding(apu2, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
		jelppi *= af::convolve3(apu / apu2, g);
		jelppi = jelppi(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::seq(inputScalars.g_dim_z + 1, inputScalars.Nz[ii] + inputScalars.g_dim_z));
		if (kk < inputScalars.deblur_iterations - 1)
			jelppi = padding(jelppi, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	}
	vec = af::flat(jelppi);
	if (inputScalars.verbose >= 3)
		mexPrint("Deconvolution computed");
}

// The initialization steps for LSQR, CGLS, CP and FISTA algorithms
// Apply PSF blurring if applicable
inline int initializationStep(Weighting& w_vec, af::array& mData, AF_im_vectors& vec, ProjectorClass& proj, scalarStruct& inputScalars,
	std::vector<int64_t> length, uint64_t m_size, const RecMethods& MethodList, uint32_t curIter, af::array& meanBP, 
	const af::array& g = af::constant(0.f, 1, 1), const uint32_t subIter = 0, const int ii = 0) {

	if (MethodList.FISTA || MethodList.FISTAL1) {
		if (curIter == 0 && subIter == 0)
			vec.uFISTA.emplace_back(vec.im_os[ii]);
		else {
			if (inputScalars.subsetsUsed == 1 || (subIter == 0 && curIter > 0))
				vec.im_os[ii] = vec.uFISTA[ii].copy();
		}
		vec.uFISTA[ii].eval();
	}
	if (curIter == 0) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Starting initialization step");
		af::sync();
		int status = 0;
		uint64_t yy = 0u;
		af::array mDataApu;
		if (MethodList.LSQR && subIter == 0) {
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing LSQR");
			vec.fLSQR.emplace_back(vec.im_os[ii].copy());
			if (ii == 0) {
				w_vec.betaLSQR = af::norm(mData);
				mData = mData / w_vec.betaLSQR;
			}
			if (DEBUG) {
				mexPrintBase("!!!!!!!!!!!!!!!mData = %f\n", af::sum<float>(mData));
				mexPrintBase("w_vec.betaLSQR = %f\n", w_vec.betaLSQR);
				mexEval();
			}
			if (inputScalars.projector_type == 6)
				backprojectionType6(mData, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, ii);
			else {
				if (inputScalars.BPType == 5) {
					mDataApu = mData;
					computeIntegralImage(inputScalars, w_vec, length[0], mData, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, mData, 0, length, m_size, meanBP, g, proj, false, ii);
				if (status != 0) {
					return -1;
				}
				af::sync();
				if (inputScalars.BPType == 5)
					mData = mDataApu;
			}
			af::sync();
			if (DEBUG) {
				mexPrintBase("!!!!!!!!!!!!!!!!!!!!!!!vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os[ii]));
				mexEval();
			}
			if (ii == inputScalars.nMultiVolumes) {
				af::array temp;
				temp = vec.rhs_os[0];
				for (int ll = 1; ll <= inputScalars.nMultiVolumes; ll++)
					temp = af::join(0, temp, vec.rhs_os[ii]);
				w_vec.alphaLSQR = af::norm(temp);
				for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++) {
					vec.im_os[ll] = vec.rhs_os[ll] / w_vec.alphaLSQR;
					vec.wLSQR.emplace_back(vec.im_os[ii].copy());
				}
				if (DEBUG) {
					mexPrintBase("!!!!!!vec.im_os = %f\n", af::sum<float>(vec.im_os[ii]));
					mexPrintBase("w_vec.alphaLSQR = %f\n", w_vec.alphaLSQR);
					mexEval();
				}
				w_vec.phiLSQR = w_vec.betaLSQR;
				w_vec.rhoLSQR = w_vec.alphaLSQR;
				af::sync();
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("LSQR initialization complete");
			}
		}
		else if (MethodList.BB && subIter == 0){
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing BB");
			if (inputScalars.projector_type == 6)
				backprojectionType6(mData, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, ii);
			else {
				if (inputScalars.BPType == 5) {
					computeIntegralImage(inputScalars, w_vec, length[0], mData, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, mData, 0, length, m_size, meanBP, g, proj, false, ii);
				if (status != 0) {
					return -1;
				}
			}
			af::sync();
			if (vec.gradBB.size() < ii + 1)
				vec.gradBB.emplace_back( -vec.rhs_os[ii].copy());
			else
				vec.gradBB[ii] = -vec.rhs_os[ii].copy();
			if (vec.imBB.size() < ii + 1)
				vec.imBB.emplace_back(vec.im_os[ii].copy());
			else
				vec.imBB[ii] = vec.im_os[ii].copy();
			if (w_vec.alphaBB.size() < ii + 1)
				w_vec.alphaBB.emplace_back(1e-4f);
			vec.im_os[ii] = vec.im_os[ii] - w_vec.alphaBB[ii] * vec.gradBB[ii];
			vec.gradBB[ii].eval();
			vec.im_os[ii].eval();
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("BB initialization complete");
		}
		else if (MethodList.CGLS && subIter == 0 && !inputScalars.largeDim) {
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing CGLS");
			if (ii == 0)
				vec.rCGLS = mData;
			mDataApu = mData.copy();
			vec.fCGLS.emplace_back(vec.im_os[ii].copy());
			if (inputScalars.projector_type == 6)
				backprojectionType6(mDataApu, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, ii);
			else {
				if (inputScalars.BPType == 5) {
					computeIntegralImage(inputScalars, w_vec, length[0], mDataApu, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, mDataApu, 0, length, m_size, meanBP, g, proj, false, ii);
				if (status != 0) {
					return -1;
				}
				af::sync();
			}
			af::sync();
			vec.im_os[ii] = vec.rhs_os[ii].copy();
			if (ii == inputScalars.nMultiVolumes) {
				for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++)
					w_vec.gammaCGLS += af::dot<float>(vec.rhs_os[ll], vec.rhs_os[ll]);
					//w_vec.gammaCGLS += af::sum<float>(vec.rhs_os[ll] * vec.rhs_os[ll]);
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("CGLS initialization complete");
			}
		}
		if (MethodList.SAGA && inputScalars.currentSubset == 0) {
			if (ii == 0)
				vec.stochasticHelper.resize(inputScalars.nMultiVolumes + 1);
			vec.SAGASum.emplace_back(af::constant(0.f, vec.im_os[ii].elements()));
			for (int uu = 0; uu < inputScalars.subsetsUsed; uu++)
				vec.stochasticHelper[ii].emplace_back(af::constant(0.f, vec.im_os[ii].elements()));
		}
		if (MethodList.CPType) {
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing PDHG algorithm");
			if (ii == 0 && !inputScalars.largeDim && inputScalars.currentSubset == 0) {
				vec.pCP.resize(inputScalars.subsetsUsed);
				for (int uu = 0; uu < inputScalars.subsetsUsed; uu++) {
					uint64_t mSize = length[uu];
					if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
						mSize = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[uu];
					if (inputScalars.listmode && inputScalars.TOF)
						vec.pCP[uu] = af::constant(0.f, mSize);
					else
						vec.pCP[uu] = af::constant(0.f, mSize * inputScalars.nBins);
					proj.memSize += (sizeof(float) * mSize * inputScalars.nBins) / 1048576ULL;
				}
			}
			else if (ii == 0 && inputScalars.largeDim)
				vec.pCP.resize(1);
			if (DEBUG) {
				mexPrintBase("subIter = %d\n", subIter);
				mexEval();
			}
			if (inputScalars.currentSubset == 0 && !inputScalars.largeDim) {
				vec.uCP.emplace_back(vec.im_os[ii].copy());
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
			}
			else if (inputScalars.currentSubset == 0 && inputScalars.largeDim)
				vec.uCP.resize(1);
			if (inputScalars.verbose >= 3)
				mexPrint("PDHG initialization complete");
		}
	}
	af::sync();
	af::deviceGC();
	return 0;
}

// Compute the specific weight needed by ACOSEM
inline int computeACOSEMWeight(scalarStruct& inputScalars, std::vector<int64_t>& length, float& uu, uint32_t osa_iter, const af::array& mData,
	uint64_t m_size, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const int64_t subSum, const af::array& g = af::constant(0.f, 1, 1)) {
	if (inputScalars.verbose >= 3)
		mexPrint("Computing ACOSEM weight");
	int status = 0;
	uu = af::sum<float>(mData);
	af::array outputFP;
	if (inputScalars.projector_type != 6) {
		if (inputScalars.listmode == 0)
			outputFP = af::constant(0.f, m_size * inputScalars.nBins);
		else
			outputFP = af::constant(0.f, m_size);
		af::sync();
		status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, g, m_size, proj);
		af::sync();
		if (status != 0)
			return -1;
	}
	else {
		outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[osa_iter]);
		forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[osa_iter], subSum, proj);
	}
	if (inputScalars.CT)
		w_vec.ACOSEM_rhs = af::sum<float>(af::exp(-outputFP));
	else
		w_vec.ACOSEM_rhs = af::sum<float>(outputFP);
	if (inputScalars.verbose >= 3)
		mexPrint("ACOSEM weight computed");
	return 0;
}

// The power method
inline int powerMethod(scalarStruct& inputScalars, Weighting& w_vec, std::vector<int64_t>& length, ProjectorClass& proj,
	AF_im_vectors& vec, const RecMethods& MethodList, const af::array& g = af::constant(0.f, 1, 1), float* F = nullptr, const float* atten = nullptr) {
	int status = 0;
	std::vector<af::array> Summ;
	af::array meanBP;
	uint64_t m_size = length[0];
	if (inputScalars.verbose > 0) {
		mexPrint("Starting power method\n");
	}
	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[0];
	std::vector<float> tauCP(inputScalars.nMultiVolumes + 1, 0.f);
	af::randomEngine r(AF_RANDOM_ENGINE_DEFAULT, 1);
	if (!inputScalars.largeDim) {
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (ii > 0 && ii % 2 == 0)
				vec.im_os[ii] = vec.im_os[ii - 1];
			else
				vec.im_os[ii] = af::abs(af::randn(inputScalars.im_dim[ii], f32, r));
			vec.im_os[ii] = vec.im_os[ii] / af::norm(vec.im_os[ii]);
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
			vec.im_os[ii].eval();
			if (DEBUG) {
				mexPrintBase("ii = %d\n", ii);
				mexEval();
			}
		}
	}
	if (!inputScalars.largeDim) {
		for (int kk = 0; kk < w_vec.powerIterations; kk++) {
			proj.memSize += (sizeof(float) * m_size) / 1048576ULL;
			af::sync();
			af::array outputFP;
			if (inputScalars.projector_type == 6) {
				outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, atten);
				outputFP.eval();
				outputFP = af::flat(outputFP);
			}
			else {
				if (inputScalars.listmode && inputScalars.TOF)
					outputFP = af::constant(0.f, m_size);
				else
					outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, g, m_size, proj, 0);
			}
			af::sync();
			if (status != 0)
				return -1;
			if (DEBUG) {
				mexPrintBase("m_size = %d\n", m_size);
				mexEval();
			}
			status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj);
			if (status != 0)
				return -1;
			computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
			if (inputScalars.projector_type == 6)
				backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, 0);
			else
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, m_size, meanBP, g, proj, false, 0);
			af::sync();
			if (status != 0)
				return -1;
			status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[0], vec.im_os[0], proj, kk, 0);
			if (status != 0)
				return -1;
			tauCP[0] = (af::dot<float>(vec.im_os[0], vec.rhs_os[0]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[0], vec.im_os[0]));
			vec.im_os[0] = vec.rhs_os[0];
			vec.im_os[0] /= af::norm(vec.im_os[0]);
			vec.im_os[0].eval();
			vec.rhs_os[0].eval();
			if ((inputScalars.verbose >= 2 || DEBUG)) {
				mexPrintBase("Largest eigenvalue for the main volume at iteration %d is %f\n", kk, tauCP[0]);
				mexEval();
			}
			proj.memSize -= (sizeof(float) * inputScalars.im_dim[0]) / 1048576ULL;
			proj.memSize -= (sizeof(float) * m_size) / 1048576ULL;
		}
		if (inputScalars.nMultiVolumes > 0) {
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				proj.memSize += (sizeof(float) * m_size) / 1048576ULL;
				af::sync();
				af::array outputFP;
				if (inputScalars.projector_type == 6)
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				else
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6) {
						forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, ii, atten);
						outputFP.eval();
						outputFP = af::flat(outputFP);
					}
					else {
						status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, g, m_size, proj, ii);
					}
					af::sync();
					if (status != 0)
						return -1;
				}
				if (DEBUG) {
					mexPrintBase("m_size = %d\n", m_size);
					mexEval();
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6)
						backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, ii);
					else
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, m_size, meanBP, g, proj, false, ii);
					af::sync();
					if (status != 0)
						return -1;
					if (ii == 0) {
						status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[ii], vec.im_os[ii], proj, kk, ii);
						if (status != 0)
							return -1;
					}
					if (ii > 0)
						tauCP[ii] = (af::dot<float>(vec.im_os[ii], vec.rhs_os[ii]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[ii], vec.im_os[ii]));
					vec.im_os[ii] = vec.rhs_os[ii];
					vec.im_os[ii] /= af::norm(vec.im_os[ii]);
					vec.im_os[ii].eval();
					vec.rhs_os[ii].eval();
					if ((inputScalars.verbose >= 2 || DEBUG) && ii > 0) {
						mexPrintBase("Largest eigenvalue for volume %d at iteration %d is %f\n", ii, kk, tauCP[ii]);
						mexEval();
					}
					proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
				}
				proj.memSize -= (sizeof(float) * m_size) / 1048576ULL;
			}
		}
	}
	else {
		for (int kk = 0; kk < w_vec.powerIterations; kk++) {
			tauCP[0] = 0.f;
			af::array outputFP = af::constant(0.f, m_size * inputScalars.nBins);
			if (inputScalars.listmode && inputScalars.TOF)
				outputFP = af::constant(0.f, m_size);
			if (DEBUG) {
				mexPrint("Starting largeDim\n");
			}
			for (int ii = 0; ii < inputScalars.subsets; ii++) {
				largeDimFirst(inputScalars, proj, ii);
				if (kk == 0) {
					vec.im_os[0] = af::abs(af::randn(inputScalars.lDimStruct.imDim[ii], f32, r));
					vec.im_os[0] = vec.im_os[0] / (af::norm(vec.im_os[0]) * static_cast<float>(inputScalars.subsets));
					vec.im_os[0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
					af::sync();
				}
				else
					vec.im_os[0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost);
				status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, g, m_size, proj, 0);
				if (status != 0)
					return -1;
				af::sync();
			}
			largeDimLast(inputScalars, proj);
			status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj);
			if (status != 0)
				return -1;
			float upper = 0.f, lower = 0.f;
			for (int ii = 0; ii < inputScalars.subsets; ii++) {
				largeDimFirst(inputScalars, proj, ii);
				vec.im_os[0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost);
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, m_size, meanBP, g, proj, false, 0);
				af::sync();
				if (status != 0)
					return -1;
				upper += af::dot<float>(vec.im_os[0], vec.rhs_os[0]);
				lower += af::dot<float>(vec.im_os[0], vec.im_os[0]);
				vec.im_os[0] = vec.rhs_os[0];
				vec.im_os[0] /= (af::norm(vec.im_os[0]) * static_cast<float>(inputScalars.subsets));
				vec.im_os[0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
				vec.im_os[0].eval();
				vec.rhs_os[0].eval();
			}
			tauCP[0] = upper * static_cast<float>(inputScalars.subsets) / lower;
			if (inputScalars.verbose >= 2 || DEBUG) {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue at iteration %d with filtering is %f\n", kk, tauCP[0]);
				else
					mexPrintBase("Largest eigenvalue at iteration %d is %f\n", kk, tauCP[0]);
				mexEval();
			}
			af::sync();
		}
	}
	std::copy(tauCP.begin(), tauCP.end(), w_vec.tauCP);
	uint32_t subsets = 1;
	if (!inputScalars.stochastic)
		subsets = inputScalars.subsets;
	if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < subsets * inputScalars.Niter) {
		bool apuM = false;
		bool apuI = false;
		if (w_vec.precondTypeMeas[1]) {
			apuM = w_vec.precondTypeMeas[1];
			w_vec.precondTypeMeas[1] = false;
		}
		else if (w_vec.precondTypeIm[5]) {
			apuI = w_vec.precondTypeIm[5];
			w_vec.precondTypeIm[5] = false;
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (ii > 0 && ii % 2 == 0)
				vec.im_os[ii] = vec.im_os[ii - 1];
			else
				vec.im_os[ii] = af::abs(af::randn(inputScalars.im_dim[ii]));
			vec.im_os[ii] = vec.im_os[ii] / af::norm(vec.im_os[ii]);
		}
		for (int kk = 0; kk < w_vec.powerIterations; kk++) {
			af::array outputFP;
			af::sync();
			if (inputScalars.projector_type == 6) {
				outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, atten);
				outputFP.eval();
				outputFP = af::flat(outputFP);
			}
			else {
				if (inputScalars.listmode && inputScalars.TOF)
					outputFP = af::constant(0.f, m_size);
				else
					outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, g, m_size, proj, 0);
			}
			af::sync();
			if (status != 0)
				return -1;
			if (DEBUG) {
				mexPrint("Power forward projection complete\n");
				mexEval();
			}
			af::sync();
			status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj);
			if (status != 0)
				return -1;
			computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
			if (inputScalars.projector_type == 6)
				backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, 0);
			else
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, m_size, meanBP, g, proj, false, 0);
			af::sync();
			if (status != 0)
				return -1;
			status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[0], vec.im_os[0], proj, kk, 0);
			if (status != 0)
				return -1;
			tauCP[0] = (af::dot<float>(vec.im_os[0], vec.rhs_os[0]) * static_cast<float>(inputScalars.subsetsUsed)) / (af::dot<float>(vec.im_os[0], vec.im_os[0]));
			vec.im_os[0] = vec.rhs_os[0];
			vec.im_os[0] /= af::norm(vec.im_os[0]);
			vec.im_os[0].eval();
			vec.rhs_os[0].eval();
			if (inputScalars.verbose >= 2 || DEBUG) {
				mexPrintBase("Largest eigenvalue for the main volume without filtering at iteration %d is %f\n", kk, tauCP[0]);
				mexEval();
			}
		}
		if (inputScalars.nMultiVolumes > 0) {
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				af::array outputFP;
				af::sync();
				if (inputScalars.projector_type == 6)
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				else
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6) {
						forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, ii, atten);
						outputFP.eval();
						outputFP = af::flat(outputFP);
					}
					else {
						status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, g, m_size, proj, ii);
					}
					af::sync();
					if (status != 0)
						return -1;
					if (DEBUG) {
						mexPrint("Power forward projection complete\n");
						mexEval();
					}
					af::sync();
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6)
						backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, 0, 0, 0, ii);
					else
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, length, m_size, meanBP, g, proj, false, ii);
					af::sync();
					if (status != 0)
						return -1;
					if (ii == 0) {
						status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[ii], vec.im_os[ii], proj, kk, ii);
						if (status != 0)
							return -1;
					}
					if (ii > 0)
						tauCP[ii] = (af::dot<float>(vec.im_os[ii], vec.rhs_os[ii]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[ii], vec.im_os[ii]));
					vec.im_os[ii] = vec.rhs_os[ii];
					vec.im_os[ii] /= af::norm(vec.im_os[ii]);
					vec.im_os[ii].eval();
					vec.rhs_os[ii].eval();
					if (inputScalars.verbose >= 2 || DEBUG) {
						mexPrintBase("Largest eigenvalue for volume %d without filtering at iteration %d is %f\n", ii, kk, tauCP[ii]);
						mexEval();
					}
				}
			}
		}
		std::copy(tauCP.begin(), tauCP.end(), w_vec.tauCP2);
		if (apuM) {
			w_vec.precondTypeMeas[1] = apuM;
		}
		else if (apuI) {
			w_vec.precondTypeIm[5] = apuI;
		}
	}
	for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		w_vec.sigmaCP[ii] = 1.f;
		if (ii > 0)
			w_vec.sigma2CP[ii] = 1.f;
		if (inputScalars.verbose > 0) {
			if (ii == 0) {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue for main volume with filtering is %f\n", w_vec.tauCP[ii]);
				else
					mexPrintBase("Largest eigenvalue for main volume is %f\n", w_vec.tauCP[ii]);
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < inputScalars.subsets * inputScalars.Niter)
					mexPrintBase("Largest eigenvalue for main volume without filtering is %f\n", w_vec.tauCP2[ii]);
			}
			else {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue for volume %d with filtering is %f\n", ii, w_vec.tauCP[ii]);
				else
					mexPrintBase("Largest eigenvalue for volume %d is %f\n", ii, w_vec.tauCP[ii]);
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < inputScalars.subsets * inputScalars.Niter)
					mexPrintBase("Largest eigenvalue for volume %d without filtering is %f\n", ii, w_vec.tauCP2[ii]);
			}
			mexEval();
		}
		w_vec.LCP.emplace_back(w_vec.tauCP[ii]);
		w_vec.tauCP[ii] = 1.f / w_vec.tauCP[ii];
		if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5])) {
			w_vec.LCP2.emplace_back(w_vec.tauCP2[ii]);
			w_vec.tauCP2[ii] = 1.f / w_vec.tauCP2[ii];
		}
	}
	return 0;
}

// Allocate vectors for proximal TV, TGV and PDHG
inline void initializeProxPriors(const RecMethods& MethodList, const scalarStruct& inputScalars, AF_im_vectors& vec) {
	if (MethodList.ProxTV || MethodList.ProxTGV) {
		vec.qProxTV.resize(3);
		std::fill(vec.qProxTV.begin(), vec.qProxTV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProxTV.size(); kk++) {
			vec.qProxTV[kk].eval();
		}
	}
	if (MethodList.ProxRDP || MethodList.ProxNLM) {
		vec.qProx.resize(1);
		std::fill(vec.qProx.begin(), vec.qProx.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProx.size(); kk++) {
			vec.qProx[kk].eval();
		}
	}
	if (MethodList.ProxTGV) {
		if (inputScalars.TGV2D) {
			vec.vProxTGV.resize(2);
			vec.qProxTGV.resize(3);
		}
		else {
			vec.vProxTGV.resize(3);
			vec.qProxTGV.resize(6);
		}
		std::fill(vec.qProxTGV.begin(), vec.qProxTGV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		std::fill(vec.vProxTGV.begin(), vec.vProxTGV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProxTGV.size(); kk++) {
			vec.qProxTGV[kk].eval();
		}
		for (int kk = 0; kk < vec.vProxTGV.size(); kk++) {
			vec.vProxTGV[kk].eval();
		}
	}
	if (MethodList.CPType && inputScalars.adaptiveType >= 1)
		vec.rhsCP.resize(inputScalars.nMultiVolumes + 1);
}

// Transfer memory control back to ArrayFire
inline void transferControl(AF_im_vectors& vec, const scalarStruct& inputScalars, const af::array& g, const Weighting& w_vec, const uint8_t compute_norm_matrix = 2, const uint8_t no_norm = 1,
	const uint32_t osa_iter = 0, const int ii = 0) {
	if (compute_norm_matrix == 1u) {
		vec.Summ[ii][0].unlock();
		if (no_norm == 0u) {
#ifdef OPENCL
			if (inputScalars.atomic_64bit)
				vec.Summ[ii][0] = vec.Summ[ii][0].as(f32) / TH;
			else if (inputScalars.atomic_32bit)
				vec.Summ[ii][0] = vec.Summ[ii][0].as(f32) / TH32;
#endif
			if (inputScalars.use_psf) {
				vec.Summ[ii][0] = computeConvolution(vec.Summ[ii][0], g, inputScalars, w_vec, 1, ii);
			}
			// Prevent division by zero
			vec.Summ[ii][0](vec.Summ[ii][0] < inputScalars.epps) = inputScalars.epps;
			vec.Summ[ii][0].eval();
			if (DEBUG) {
				mexPrint("Sens image steps 1 done\n");
			}
		}
	}
	else if (compute_norm_matrix == 2) {
		vec.Summ[ii][osa_iter].unlock();
		if (no_norm == 0u) {
#ifdef OPENCL
			if (inputScalars.atomic_64bit) {
				vec.Summ[ii][osa_iter] = vec.Summ[ii][osa_iter].as(f32) / TH;
			}
			else if (inputScalars.atomic_32bit) {
				vec.Summ[ii][osa_iter] = vec.Summ[ii][osa_iter].as(f32) / TH32;
			}
#endif
			if (inputScalars.use_psf) {
				vec.Summ[ii][osa_iter] = computeConvolution(vec.Summ[ii][osa_iter], g, inputScalars, w_vec, 1, ii);
				af::sync();
			}
			vec.Summ[ii][osa_iter](vec.Summ[ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
			vec.Summ[ii][osa_iter].eval();
			if (DEBUG) {
				mexPrint("Sens image steps 2 done\n");
			}
			if (DEBUG) {
				mexPrintBase("inputScalars.epps = %f\n", inputScalars.epps);
				mexPrintBase("min(Summ) = %f\n", af::min<float>(vec.Summ[ii][osa_iter]));
				mexEval();
			}
		}
	}
	if (DEBUG && inputScalars.atomic_64bit) {
		mexPrintBase("min(rhs_os) = %d\n", af::min<int64_t>(vec.rhs_os[ii]));
		mexPrintBase("inputScalars.atomic_64bit = %d\n", inputScalars.atomic_64bit);
		mexEval();
	}
}