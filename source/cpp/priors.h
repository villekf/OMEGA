#pragma once
#include "functions.hpp"

inline int MRP(const af::array& im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const scalarStruct& inputScalars,
	ProjectorClass& proj, af::array& dU, const float beta, const bool med_no_norm = false) {
	int status = 0;
	// Prevent out of bounds cases for the (rare) 2D cases
	const uint32_t medzEff = inputScalars.Nz[0] == 1 ? 0U : medz;
	af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], medx, medy, medzEff);
	af::array grad = af::constant(0.f, im.elements());
	status = MRPAF(padd, grad, inputScalars, proj, medx, medy, medzEff);
	if (status != 0)
		return -1;
	if (med_no_norm)
		dU += (im - grad) * beta;
	else
		dU += (im - (grad)) / (grad + inputScalars.epps) * beta;
	if (DEBUG) {
		mexPrintBase("min(grad2) = %f\n", af::min<float>(grad));
		mexPrintBase("grad2 = %f\n", af::sum<float>(grad));
		mexPrintBase("min(dU) = %f\n", af::min<float>(dU));
		mexPrintBase("max(dU) = %f\n", af::max<float>(dU));
		mexPrintBase("min(im) = %f\n", af::min<float>(im));
		mexPrintBase("max(im) = %f\n", af::max<float>(im));
	}
	return status;
}

inline af::array Quadratic_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_quad)
{
	const af::array apu_pad = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
	af::array grad;
	af::array weights = weights_quad;
	if (Ndz == 0 || inputScalars.Nz[0] == 1) {
		grad = af::convolve2(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::span);
	}
	else {
		grad = af::convolve3(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz - 1));
	}
	grad = af::flat(grad);
	return grad;
}

// Huber prior gradient: grad_j = sum_k w_k * clip(x_j - x_k, -delta, delta), where k runs over the neighbours of voxel j.
// The weight kernel has the sum of the weights in the centre and -w_k elsewhere, so w_k = -kernel(k). The clipping has to be
// applied to each difference separately (clipping the weighted sum of the differences is not the Huber gradient)
inline af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_huber, const float delta)
{
	const af::array apu_pad = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
	const dim_t Nx = inputScalars.Nx[0], Ny = inputScalars.Ny[0], Nz = inputScalars.Nz[0];
	const af::array im3d = af::moddims(im, Nx, Ny, Nz);
	const bool is2D = (Ndz == 0 || inputScalars.Nz[0] == 1);
	const dim_t kx = 2 * Ndx + 1, ky = 2 * Ndy + 1, kz = 2 * Ndz + 1;
	// The kernel is tiny, so it is copied to the host once
	std::vector<float> kernel(weights_huber.elements());
	weights_huber.as(f32).host(kernel.data());
	af::array grad = af::constant(0.f, Nx, Ny, Nz);
	// Largest absolute difference, accumulated on the device and reduced once at the end (verbose output only)
	af::array maxAbs;
	// A 2D image uses only the central axial slice of the kernel
	const dim_t kBegin = is2D ? static_cast<dim_t>(Ndz) : 0;
	const dim_t kEnd = is2D ? static_cast<dim_t>(Ndz) : kz - 1;
	for (dim_t k = kBegin; k <= kEnd; k++) {
		for (dim_t j = 0; j < ky; j++) {
			for (dim_t i = 0; i < kx; i++) {
				if (i == static_cast<dim_t>(Ndx) && j == static_cast<dim_t>(Ndy) && k == static_cast<dim_t>(Ndz))
					continue;
				const float w = -kernel[i + j * kx + k * kx * ky];
				if (w == 0.f)
					continue;
				af::array shifted;
				if (is2D)
					shifted = apu_pad(af::seq(static_cast<double>(i), static_cast<double>(i + Nx - 1)), af::seq(static_cast<double>(j), static_cast<double>(j + Ny - 1)), af::span);
				else
					shifted = apu_pad(af::seq(static_cast<double>(i), static_cast<double>(i + Nx - 1)), af::seq(static_cast<double>(j), static_cast<double>(j + Ny - 1)),
						af::seq(static_cast<double>(k), static_cast<double>(k + Nz - 1)));
				const af::array diff = im3d - shifted;
				if (inputScalars.verbose > 0)
					maxAbs = maxAbs.isempty() ? af::abs(diff) : (af::max)(maxAbs, af::abs(diff));
				grad += w * af::clamp(diff, -static_cast<double>(delta), static_cast<double>(delta));
			}
		}
	}
	const float maxDiff = (inputScalars.verbose > 0 && !maxAbs.isempty()) ? af::max<float>(maxAbs) : 0.f;
	if (inputScalars.verbose > 0 && maxDiff > 0.f && delta >= maxDiff)
		mexPrint("Delta value of Huber prior larger than all the pixel difference values\n");
	return af::flat(grad);
}

inline af::array FMH(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array& offsets, const af::array& fmh_weights, const uint32_t alku_fmh, const bool med_no_norm = false)
{
	af::array grad;
	af::array indeksi1;
	const af::array padd = af::flat(padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz));
	uint32_t luup;
	const bool is2DFMH = (inputScalars.Nz[0] == 1 || Ndz == 0);
	if (is2DFMH) {
		grad = af::constant(0.f, inputScalars.im_dim[0], 5);
		luup = 4;
	}
	else {
		grad = af::constant(0.f, inputScalars.im_dim[0], 14);
		luup = 13;
	}
	for (uint32_t ii = 0; ii < luup; ii++) {
		// Line direction ii in the unit cube (2D: 4 directions, 3D: 13 directions). The neighborhood is the
		// (2Ndx+1)x(2Ndy+1)x(2Ndz+1) box with x fastest, so the columns of the line are the centre plus k * d, k = N...-N
		const int dx = static_cast<int>(ii % 3U) - 1;
		const int dy = is2DFMH ? static_cast<int>(ii / 3U) - 1 : static_cast<int>((ii / 3U) % 3U) - 1;
		const int dz = is2DFMH ? 0 : static_cast<int>(ii / 9U) - 1;
		std::vector<uint32_t> colIdx;
		colIdx.reserve(2U * Ndx + 1U);
		for (int k = static_cast<int>(Ndx); k >= -static_cast<int>(Ndx); k--)
			colIdx.push_back(static_cast<uint32_t>((static_cast<int>(Ndx) + k * dx) + (static_cast<int>(Ndy) + k * dy) * static_cast<int>(2U * Ndx + 1U)
				+ (static_cast<int>(Ndz) + k * dz) * static_cast<int>((2U * Ndx + 1U) * (2U * Ndy + 1U))));
		const af::array colArr(static_cast<dim_t>(colIdx.size()), colIdx.data());
		indeksi1 = af::flat(offsets(af::span, colArr));
		af::array apu_pad = af::moddims(padd(indeksi1 + 0), inputScalars.im_dim[0], fmh_weights.dims(0));
		grad(af::span, ii) = af::matmul(apu_pad, fmh_weights(af::span, ii));
	}
	indeksi1 = offsets.col(alku_fmh);
	grad(af::span, af::end) = padd(indeksi1 + 0U);
	grad = af::median(grad, 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

inline af::array L_filter(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& offsets, const af::array& a_L, const bool med_no_norm = false)
{
	af::array grad;
	af::array apu_pad = af::flat(padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz));
	apu_pad = apu_pad(af::flat(offsets));
	apu_pad = af::sort(af::moddims(apu_pad, inputScalars.im_dim[0], a_L.dims(0)), 1);
	grad = af::sum(af::batchFunc(apu_pad, af::transpose(a_L), batchMul), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

inline af::array Weighted_mean(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& weighted_weights, const float w_sum, const uint32_t mean_type = 1, const bool med_no_norm = false)
{
	af::array grad = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
	// w_sum is the precomputed sum of weighted_weights (set together with the weights in the prepass);
	// the sum is recomputed only if it is not available
	const float wsum = (std::isfinite(w_sum) && w_sum != 0.f) ? w_sum : af::sum<float>(af::flat(weighted_weights));
	if (mean_type == 1U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = af::convolve2(padd, weighted_weights / wsum);
		else
			grad = af::convolve3(padd, weighted_weights / wsum);
	}
	else if (mean_type == 2U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = 1.f / af::convolve2(1.f / padd, weighted_weights / wsum);
		else
			grad = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
	}
	else if (mean_type == 3U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = af::exp(af::convolve2(af::log(padd), weighted_weights / wsum));
		else
			grad = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
	}
	else if (mean_type == 4U) {
		grad = af::constant(0.f, im.dims(0));
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::convolve3(padd, weighted_weights / wsum);
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		grad = batchFunc(im, mm, batchDiv);
		grad(grad < 1e-6f) = 1e-6f;
		grad = batchFunc(af::log(grad), af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 5U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-6f) = 1e-6f;
		grad = batchFunc(batchFunc(im_t, mm, batchMinus), im_t, batchDiv);
		grad = grad - af::pow(batchFunc((batchFunc(im_t, mm, batchMinus)), std::sqrt(2.f) * (im_t), batchDiv), 2.);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 6U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-6f) = 1e-6f;
		grad = 1.f - batchFunc(mm, im_t, batchDiv);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else {
		mexWarning("Unsupported mean type");
	}
	if (mean_type <= 3U) {
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::span);
		}
		else {
			grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz - 1));
		}
		grad = af::flat(grad);
		if (med_no_norm)
			grad = im - grad;
		else
			grad = (im - grad) / (grad + inputScalars.epps);
	}
	return grad;
}

inline af::array AD(const af::array& im, const scalarStruct& inputScalars, const float TimeStepAD, const float KAD, const uint32_t NiterAD,
	const af_flux_function FluxType = AF_FLUX_EXPONENTIAL, const af_diffusion_eq DiffusionType = AF_DIFFUSION_GRAD, const bool med_no_norm = false)
{
	const af::array padd = af::moddims(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
	af::array grad = af::anisotropicDiffusion(padd, TimeStepAD, KAD, NiterAD, FluxType, DiffusionType);
	grad = af::flat(grad);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}


// Compute the TV prior
inline int TVprior(const scalarStruct& inputScalars, const TVdata& S, const af::array& ima, const Weighting& w_vec, ProjectorClass& proj, af::array& dU, 
	const float beta, const int kk = 0) {
	af::array gradi;
	int status = 0;

	status = TVAF(dU, ima, inputScalars, w_vec.data.SATVPhi, S, proj, beta, w_vec, kk);

	return status;
}

// Compute the hyperbolic prior
inline int hyperbolic(const scalarStruct& inputScalars, const af::array& ima, const Weighting& w_vec, ProjectorClass& proj, af::array& dU, const float beta, const int kk = 0) {
	af::array gradi;
	int status = 0;

	status = hyperAF(dU, ima, inputScalars, w_vec.data.SATVPhi, proj, beta, w_vec, kk);

	return status;
}

inline int proxTV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const float beta, const uint32_t timestep) {
	int status = 0;
#ifndef CPU
	status = proxTVGradAF(im, vec.qProxTV, inputScalars, w_vec.sigma2CP[timestep][0], beta, vec.vProxTGV, proj);
	if (status != 0)
		return -1;
	status = proxTVDivAF(vec.qProxTV, dU, inputScalars, proj, timestep);
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("beta = %f\n", beta);
		mexPrintBase("w_vec.sigma2CP[timestep] = %f\n", w_vec.sigma2CP[timestep][0]);
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexEval();
	}
#else
	mexPrint("Proximal TV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}

inline int proxTGV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const uint32_t osa_iter, const uint32_t timestep) {
	int status = 0;
#ifndef CPU
	status = proxTV(im, inputScalars, vec, proj, w_vec, dU, w_vec.alpha0CPTGV, timestep);
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexPrintBase("vec.qProxTGV = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVSymmDerivAF(vec.vProxTGV, vec.qProxTGV, inputScalars, w_vec.sigma2CP[timestep][0], w_vec.alpha1CPTGV, proj);
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("vec.qCPTGV2 = %f\n", af::sum<float>(vec.qProxTGV[2]));
		mexEval();
	}
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("w_vec.alpha0CPTGV = %f\n", w_vec.alpha0CPTGV);
		mexPrintBase("w_vec.alpha1CPTGV = %f\n", w_vec.alpha1CPTGV);
		mexPrintBase("w_vec.sigma2CP[timestep] = %f\n", w_vec.sigma2CP[timestep][0]);
		mexPrintBase("osa_iter = %d\n", osa_iter);
		mexPrintBase("vec.qProxTGV0 = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVDivAF(vec.qProxTGV, vec.vProxTGV, vec.qProxTV, inputScalars, w_vec.thetaCP[timestep][osa_iter], w_vec.tauCP[timestep][0], proj);
#else
	mexPrint("Proximal TGV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}

inline int RDP(const af::array& im, const scalarStruct& inputScalars, const float gamma, ProjectorClass& proj, af::array& dU, const float beta, const af::array& RDPref, const Weighting& w_vec, const bool RDPLargeNeighbor = false, const bool useRDPRef = false, const int kk = 0) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("im_RDP = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_RDP) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = RDPAF(dU, im, inputScalars, gamma, proj, beta, RDPref, w_vec, RDPLargeNeighbor, useRDPRef, kk);
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("grad = %f\n", af::sum<float>(dU));
		mexPrintBase("min(grad) = %f\n", af::min<float>(dU));
		mexPrintBase("max(grad) = %f\n", af::max<float>(dU));
		mexEval();
	}
	return status;
}

inline int GGMRF(const af::array& im, const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, ProjectorClass& proj, af::array& dU, const Weighting& w_vec, const float beta, const int kk = 0) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("im_GGMRF = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_GGMRF) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = GGMRFAF(dU, im, inputScalars, p, q, c, pqc, proj, beta, w_vec, kk);
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("grad = %f\n", af::sum<float>(dU));
		mexEval();
	}
	return status;
}



inline int NLM(ProjectorClass& proj, const af::array& im, Weighting& w_vec, const scalarStruct& inputScalars, af::array& dU, const float beta, const int kk = 0)
{
	int status = 0;
	status = NLMAF(dU, im, inputScalars, w_vec, proj, beta, kk);
	return status;
}

// Curvature bound of the prior for the separable-surrogate step used by SPS and the image preconditioner type 6.
// For priors whose gradient is linear (or, for Huber, piecewise linear with slope at most one) in the image, grad = K * x with
// a symmetric kernel K, and the Hessian of the prior is bounded by the Gershgorin bound sum(|K|) (the row sum of the absolute values).
// This is added to the data curvature in the step size. Priors without a known bound return zero and set unbounded, in which
// case the step uses only the data curvature. Proximal priors have no gradient step and need no bound.
inline float separablePriorCurvature(const RecMethods& MethodList, const Weighting& w_vec, bool& unbounded)
{
	unbounded = false;
	float kappa = 0.f;
	// The order follows applySpatialPrior, where only the first selected prior is applied
	if (MethodList.MRP)
		unbounded = true;
	else if (MethodList.Quad || MethodList.Huber) {
		const af::array& kernel = MethodList.Quad ? w_vec.weights_quad : w_vec.weights_huber;
		std::vector<float> K(kernel.elements());
		kernel.as(f32).host(K.data());
		for (const float v : K)
			if (std::isfinite(v))
				kappa += std::abs(v);
	}
	else if (MethodList.WeightedMean) {
		// With med_no_norm and the arithmetic mean the gradient is x - conv(x, w / sum(w)), i.e. K = delta - w / sum(w),
		// where w is the full weight kernel including the centre element
		if (w_vec.med_no_norm && w_vec.mean_type == 1U) {
			std::vector<float> K(w_vec.weighted_weights.elements());
			w_vec.weighted_weights.as(f32).host(K.data());
			float wsum = 0.f;
			for (const float v : K)
				if (std::isfinite(v))
					wsum += v;
			if (wsum != 0.f) {
				const size_t centre = K.size() / 2;
				for (size_t i = 0; i < K.size(); i++) {
					if (!std::isfinite(K[i]))
						continue;
					const float val = K[i] / wsum;
					kappa += std::abs((i == centre ? 1.f : 0.f) - val);
				}
			}
			else
				unbounded = true;
		}
		else
			unbounded = true;
	}
	else if (MethodList.L || MethodList.FMH || MethodList.TV || MethodList.hyperbolic || MethodList.AD || MethodList.APLS
		|| MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.GGMRF)
		unbounded = true;
	return kappa;
}


inline int applySpatialPrior(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, ProjectorClass& proj, const float beta, const uint32_t timestep, const uint32_t osa_iter = 0, const uint8_t compute_norm_matrix = 0, 
	const bool iter = false, const int kk = 0) {
	af::array* dU = nullptr;
	int status = 0;
	// The prior gradient is always computed into its own zero-initialised buffer vec.dU, which the caller adds into the gradient
	// (or the one-step-late/BSREM/SART/POCS code reads directly). The buffer is sized like the rhs of the main volume, which equals
	// the image size except for a largeDim chunk, where the prior kernels write the unpadded chunk into the rhs-sized buffer.
	// The iteration-level (BSREM/ROSEM-MAP), SART and POCS calls have no rhs of that size and use the image size
	const bool imageSized = iter || MethodList.POCS || MethodList.SART || vec.rhs_os[timestep].empty() || vec.rhs_os[timestep][0].elements() <= 1;
	vec.dU[timestep] = af::constant(0.f, imageSized ? vec.im_os[timestep][0].elements() : vec.rhs_os[timestep][0].elements());
	dU = &vec.dU[timestep];
	if (MethodList.MRP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing MRP gradient");
		status = MRP(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, proj, *dU, beta, w_vec.med_no_norm);
		if (status != 0)
			return -1;
	}
	else if (MethodList.Quad) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing quadratic prior gradient");
		*dU += beta * Quadratic_prior(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_quad);
	}
	else if (MethodList.Huber) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing Huber prior gradient");
		*dU += beta * Huber_prior(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_huber, w_vec.huber_delta);
	}
	else if (MethodList.L) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing L-filter gradient");
		*dU += beta * L_filter(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets,
			w_vec.a_L, w_vec.med_no_norm);
	}
	else if (MethodList.FMH) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing FMH prior gradient");
		*dU += beta * FMH(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets,
			w_vec.fmh_weights, w_vec.alku_fmh, w_vec.med_no_norm);
	}
	else if (MethodList.WeightedMean) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing weighted mean prior gradient");
		*dU += beta * Weighted_mean(vec.im_os[timestep][0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weighted_weights,
			w_vec.w_sum, w_vec.mean_type, w_vec.med_no_norm);
	}
	else if (MethodList.TV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TV prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[timestep][0], w_vec, proj, *dU, beta, kk);
		if (status != 0)
			return -1;
	}
	else if (MethodList.hyperbolic) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing hyperbolic prior gradient");
		status = hyperbolic(inputScalars, vec.im_os[timestep][0], w_vec, proj, *dU, beta, kk);
		if (status != 0)
			return -1;
	}
	else if (MethodList.AD) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing AD prior gradient");
		if (osa_iter != 0u) {
			*dU += beta * AD(vec.im_os[timestep][0], inputScalars, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
				w_vec.DiffusionType, w_vec.med_no_norm);
		}
	}
	else if (MethodList.APLS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing APLS prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[timestep][0], w_vec, proj, *dU, beta, kk);
		if (status != 0)
			return -1;
	}
	else if (MethodList.ProxTGV || MethodList.TGV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TGV prior");
		if (osa_iter >= 100)
			w_vec.sigma2CP[timestep] = w_vec.sigmaCP[timestep];
		status = proxTGV(vec.im_os[timestep][0], inputScalars, vec, proj, w_vec, *dU, osa_iter, timestep);
		if (status != 0)
			return -1;
	}
	else if (MethodList.ProxTV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing proximal TV prior");
		status = proxTV(vec.im_os[timestep][0], inputScalars, vec, proj, w_vec, *dU, w_vec.betaReg, timestep);
		if (status != 0)
			return -1;
	}
	else if (MethodList.NLM) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing NLM prior gradient");
		status = NLM(proj, vec.im_os[timestep][0], w_vec, inputScalars, *dU, beta, kk);
		if (status != 0)
			return -1;
	}
	else if (MethodList.RDP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing RDP prior gradient");
		// RDPref is empty unless the RDP reference precompute ran
		const af::array emptyRef;
		status = RDP(vec.im_os[timestep][0], inputScalars, w_vec.RDP_gamma, proj, *dU, beta, (w_vec.RDPref.size() > timestep) ? w_vec.RDPref[timestep] : emptyRef, w_vec, w_vec.RDPLargeNeighbor, w_vec.RDP_anatomical, kk);
		if (status != 0)
			return -1;
	}
	else if (MethodList.GGMRF) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing GGMRF prior gradient");
		status = GGMRF(vec.im_os[timestep][0], inputScalars, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c, w_vec.GGMRF_pqc, proj, *dU, w_vec, beta, kk);
		if (status != 0)
			return -1;
	}
	if (inputScalars.verbose >= 3 && (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV 
		|| MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.ProxTGV 
		|| MethodList.ProxTV ||  MethodList.GGMRF))
		mexPrint("Spatial prior computed");
	af::eval(*dU);
	return status;
}

inline int applyTemporalPrior(
    AF_im_vectors& vec,
    Weighting& w_vec,
    const RecMethods& MethodList,
    const scalarStruct& inputScalars,
    ProjectorClass& proj,
    const bool iter = false,
    const std::vector<af::array>* neighbours = nullptr,
    const int onlyTimestep = -1,
    const float betaScale = 1.f
) {
	int status = 0;
    // The temporal prior gradient is always written into its own buffer vec.dUt (positive beta, scaled with betaScale), the caller
    // adds it into the gradient with the sign of the update direction. If onlyTimestep >= 0, only that timestep is computed. If
    // neighbours is given, the frames t - 1 and t + 1 are taken from it instead of vec.im_os (the own frame always comes from vec.im_os)
    const float betaT = w_vec.beta_temporal * betaScale;
    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        if (onlyTimestep >= 0 && static_cast<uint32_t>(onlyTimestep) != timestep)
            continue;
        // The buffer has the size of the image (the same as the rhs of the main volume)
        af::array grad = af::constant(0.f, vec.im_os[timestep][0].elements());
        const af::array xCur = vec.im_os[timestep][0];
        if (MethodList.TemporalSmoothness) { // First derivative smoothness prior
            if (inputScalars.verbose >= 3)
                mexPrint("Computing temporal smoothness prior gradient");
            if (timestep > 0) {
                const af::array xPrev = neighbours ? (*neighbours)[timestep-1] : vec.im_os[timestep-1][0];
                grad += (xCur - xPrev);
            }
            if (timestep < inputScalars.Nt - 1) {
                const af::array xNext = neighbours ? (*neighbours)[timestep+1] : vec.im_os[timestep+1][0];
                grad += (xCur - xNext);
            }
        } else if (MethodList.TemporalTV) { // Total variations prior
            if (inputScalars.verbose >= 3)
                mexPrintBase("Computing temporal total variations prior gradient for timestep %d\n", timestep);
            if (timestep > 0) {
                const af::array xPrev = neighbours ? (*neighbours)[timestep-1] : vec.im_os[timestep-1][0];
                af::array backward_diff = (xCur - xPrev);
                grad += (backward_diff) / (af::abs(backward_diff) + w_vec.TemporalTVsmoothing);
            }
            if (timestep < inputScalars.Nt - 1) {
                const af::array xNext = neighbours ? (*neighbours)[timestep+1] : vec.im_os[timestep+1][0];
                af::array forward_diff = (xCur - xNext);
                grad += (forward_diff) / (af::abs(forward_diff) + w_vec.TemporalTVsmoothing);
            }
        }

        vec.dUt[timestep] = betaT * grad;
        af::eval(vec.dUt[timestep]);
    }

    if (inputScalars.verbose >= 3 && (MethodList.TemporalSmoothness))
        mexPrint("Temporal prior computed");

    return status;
}