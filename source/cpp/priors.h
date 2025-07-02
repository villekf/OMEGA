#pragma once
#include "functions.hpp"

inline int MRP(const af::array& im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const scalarStruct& inputScalars,
	ProjectorClass& proj, af::array& dU, const float beta, const bool med_no_norm = false) {
	int status = 0;
	af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], medx, medy, medz);
	af::array grad = af::constant(0.f, im.elements());
	status = MRPAF(padd, grad, inputScalars, proj, medx, medy, medz);
	if (status != 0)
		return -1;
	if (med_no_norm)
		dU += (im - grad) * beta;
	else
		dU += (im - (grad)) / (grad + inputScalars.epps) * beta;
	af::sync();
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

inline af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_huber, const float delta)
{
	af::array grad = Quadratic_prior(im, Ndx, Ndy, Ndz, inputScalars, weights_huber);
	if (af::sum<dim_t>(delta >= af::abs(af::flat(grad))) == grad.elements() && af::sum<int>(af::flat(grad)) != 0 && inputScalars.verbose > 0)
		mexPrint("Delta value of Huber prior larger than all the pixel difference values\n");
	grad(grad > delta) = delta;
	grad(grad < -delta) = -delta;
	return grad;
}

inline af::array FMH(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array& offsets, const af::array& fmh_weights, const uint32_t alku_fmh, const bool med_no_norm = false)
{
	af::array grad;
	af::array indeksi1;
	const af::array padd = af::flat(padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz));
	uint32_t luup;
	if (inputScalars.Nz[0] == 1 || Ndz == 0) {
		grad = af::constant(0.f, inputScalars.im_dim[0], 5);
		luup = 4;
	}
	else {
		grad = af::constant(0.f, inputScalars.im_dim[0], 14);
		luup = 13;
	}
	for (uint32_t ii = 0; ii < luup; ii++) {
		indeksi1 = af::flat(offsets(af::span, af::seq(Ndx * ii, offsets.dims(1) - Ndx * (ii)-1, alku_fmh / Ndx - ii)));
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
	const float wsum = af::sum<float>(af::flat(weighted_weights));
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

inline int proxTV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const float beta) {
	int status = 0;
#ifndef CPU
	status = proxTVGradAF(im, vec.qProxTV, inputScalars, w_vec.sigma2CP[0], vec.vProxTGV, proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTVQAF(vec.qProxTV, beta, proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTVDivAF(vec.qProxTV, dU, inputScalars, proj);
	if (DEBUG) {
		mexPrintBase("beta = %f\n", beta);
		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP[0]);
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexEval();
	}
#else
	mexPrint("Proximal TV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}

inline int proxTGV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const uint32_t osa_iter = 0) {
	int status = 0;
#ifndef CPU
	status = proxTV(im, inputScalars, vec, proj, w_vec, dU, w_vec.alpha0CPTGV);
	if (DEBUG) {
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexPrintBase("vec.qProxTGV = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVSymmDerivAF(vec.vProxTGV, vec.qProxTGV, inputScalars, w_vec.sigma2CP[0], proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTGVQAF(vec.qProxTGV, inputScalars, w_vec.alpha1CPTGV, proj);
	if (DEBUG) {
		mexPrintBase("vec.qCPTGV2 = %f\n", af::sum<float>(vec.qProxTGV[2]));
		mexEval();
	}
	af::sync();
	if (status != 0)
		return -1;
	if (DEBUG) {
		mexPrintBase("w_vec.alpha0CPTGV = %f\n", w_vec.alpha0CPTGV);
		mexPrintBase("w_vec.alpha1CPTGV = %f\n", w_vec.alpha1CPTGV);
		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP[0]);
		mexPrintBase("osa_iter = %d\n", osa_iter);
		mexPrintBase("vec.qProxTGV0 = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVDivAF(vec.qProxTGV, vec.vProxTGV, vec.qProxTV, inputScalars, w_vec.thetaCP[osa_iter], w_vec.tauCP[0], proj);
#else
	mexPrint("Proximal TGV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}

inline int RDP(const af::array& im, const scalarStruct& inputScalars, const float gamma, ProjectorClass& proj, af::array& dU, const float beta, const af::array& RDPref, const Weighting& w_vec,
	const bool RDPLargeNeighbor = false, const bool useRDPRef = false, const int kk = 0) {

	int status = 0;
	af::sync();
	if (DEBUG) {
		mexPrintBase("im_RDP = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_RDP) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = RDPAF(dU, im, inputScalars, gamma, proj, beta, RDPref, w_vec, RDPLargeNeighbor, useRDPRef, kk);
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
	af::sync();
	if (DEBUG) {
		mexPrintBase("im_GGMRF = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_GGMRF) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = GGMRFAF(dU, im, inputScalars, p, q, c, pqc, proj, beta, w_vec, kk);
	if (DEBUG) {
		mexPrintBase("grad = %f\n", af::sum<float>(dU));
		mexEval();
	}
	return status;
}



inline int NLM(ProjectorClass& proj, const af::array& im, Weighting& w_vec, const scalarStruct& inputScalars, af::array& dU, const float beta, const int kk = 0)
{
	int status = 0;
	af::sync();
	status = NLMAF(dU, im, inputScalars, w_vec, proj, beta, kk);
	return status;
}


inline int applyPrior(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, ProjectorClass& proj, const float beta, const uint32_t osa_iter = 0, const uint8_t compute_norm_matrix = 0, 
	const bool iter = false, const int kk = 0) {
	af::array* dU = nullptr;
	int status = 0;
	if (iter)
		dU = &vec.im_os[0];
	else if (MethodList.RBIOSL || MethodList.OSLOSEM || MethodList.OSLCOSEM || MethodList.POCS || MethodList.SAGA || MethodList.SART) {
		vec.dU = af::constant(0.f, vec.im_os[0].elements());
		dU = &vec.dU;
	}
	else
		dU = &vec.rhs_os[0];
	if (MethodList.MRP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing MRP gradient");
		status = MRP(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, proj, *dU, beta, w_vec.med_no_norm);
	}
	else if (MethodList.Quad) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing quadratic prior gradient");
		*dU += beta * Quadratic_prior(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_quad);
	}
	else if (MethodList.Huber) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing Huber prior gradient");
		*dU += beta * Huber_prior(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_huber, w_vec.huber_delta);
	}
	else if (MethodList.L) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing L-filter gradient");
		*dU += beta * L_filter(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets,
			w_vec.a_L, w_vec.med_no_norm);
	}
	else if (MethodList.FMH) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing FMH prior gradient");
		*dU += beta * FMH(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets,
			w_vec.fmh_weights, w_vec.alku_fmh, w_vec.med_no_norm);
	}
	else if (MethodList.WeightedMean) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing weighted mean prior gradient");
		*dU += beta * Weighted_mean(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weighted_weights,
			w_vec.w_sum, w_vec.mean_type, w_vec.med_no_norm);
	}
	else if (MethodList.TV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TV prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[0], w_vec, proj, *dU, beta, kk);
	}
	else if (MethodList.hyperbolic) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing hyperbolic prior gradient");
		status = hyperbolic(inputScalars, vec.im_os[0], w_vec, proj, *dU, beta, kk);
	}
	else if (MethodList.AD) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing AD prior gradient");
		if (osa_iter == 0u) {
			*dU += af::constant(0.f, inputScalars.im_dim[0], 1);
		}
		else {
			*dU += beta * AD(vec.im_os[0], inputScalars, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
				w_vec.DiffusionType, w_vec.med_no_norm);
		}
	}
	else if (MethodList.APLS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing APLS prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[0], w_vec, proj, *dU, beta, kk);
	}
	else if (MethodList.ProxTGV || MethodList.TGV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TGV prior");
		if (osa_iter >= 100)
			w_vec.sigma2CP = w_vec.sigmaCP;
		status = proxTGV(vec.im_os[0], inputScalars, vec, proj, w_vec, *dU, osa_iter);
	}
	else if (MethodList.ProxTV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing proximal TV prior");
		status = proxTV(vec.im_os[0], inputScalars, vec, proj, w_vec, *dU, w_vec.betaReg);
	}
	else if (MethodList.NLM) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing NLM prior gradient");
		status = NLM(proj, vec.im_os[0], w_vec, inputScalars, *dU, beta, kk);
	}
	else if (MethodList.RDP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing RDP prior gradient");
		status = RDP(vec.im_os[0], inputScalars, w_vec.RDP_gamma, proj, *dU, beta, w_vec.RDPref, w_vec, w_vec.RDPLargeNeighbor, w_vec.RDP_anatomical, kk);
	}
	else if (MethodList.GGMRF) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing GGMRF prior gradient");
		status = GGMRF(vec.im_os[0], inputScalars, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c, w_vec.GGMRF_pqc, proj, *dU, w_vec, beta, kk);
	}
	af::deviceGC();
	if (inputScalars.verbose >= 3 && (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV 
		|| MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.ProxTGV 
		|| MethodList.ProxTV || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF))
		mexPrint("Prior computed");
	af::eval(*dU);
	return status;
}