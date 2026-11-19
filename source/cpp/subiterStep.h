#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"
#include <utility>
#include <vector>

// Swaps the main-volume estimate with the PDHG primal iterate (uCP) of the given timestep, see computeOSEstimates.
// Only swaps when both hold the same number of elements; returns whether the swap was made.
inline bool swapPDHGPrimal(AF_im_vectors& vec, const uint32_t timestep) {
    if (vec.uCP.size() <= timestep || vec.uCP[timestep].empty() || vec.im_os[timestep].empty())
        return false;
    if (vec.uCP[timestep][0].elements() != vec.im_os[timestep][0].elements())
        return false;
    std::swap(vec.im_os[timestep][0], vec.uCP[timestep][0]);
    return true;
}

// Writes the complete-data image of the current subset into the COSEM matrix and keeps the running sum
// (sum over the subsets) of the matrix up to date. The sum is resynchronised exactly at the first subset of each
// iteration to prevent float drift.
inline void updateCOSEMCompleteData(AF_im_vectors& vec, const int ii, const uint32_t osa_iter, const af::array& newColExpr) {
    af::array newCol = newColExpr;
    newCol.eval();
    // The resynchronisation must be computed before column 0 is overwritten
    if (osa_iter == 0)
        vec.C_coSum[ii] = af::sum(vec.C_co[ii], 1);
    vec.C_coSum[ii] += newCol - vec.C_co[ii](af::span, osa_iter);
    // Evaluated before the column is overwritten, so that the old column is read
    vec.C_coSum[ii].eval();
    vec.C_co[ii](af::span, osa_iter) = newCol;
    vec.C_co[ii].eval();
}

inline int computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t iter, uint32_t osa_iter, scalarStruct& inputScalars, std::vector<int64_t>& length, bool& break_iter, const int64_t* pituus, const af::array& g, ProjectorClass& proj, const std::vector<af::array> &mData, uint64_t m_size, const int64_t subSum, const uint8_t compute_norm_matrix, const int kk = 0, const bool largeDim = false, const int vv = 0) {
    // Ignore this function if using fastPDHG
    if (inputScalars.fastPDHG)
        return 0;
    int status = 0;
    af::array OSEMApu, COSEMApu;
    if (DEBUG || inputScalars.verbose >= 3) {
        proj.tStartLocal = std::chrono::steady_clock::now();
    }

    if (DEBUG) {
        mexPrint("Algo start\n");
        mexEval();
    }

    bool MAP = false;
    if (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV || MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.hyperbolic || MethodList.ProxTGV || MethodList.ProxTV || MethodList.GGMRF || MethodList.TemporalSmoothness || MethodList.TemporalTV)
        MAP = true;

    // PDHG without subsets: the forward projection uses the extrapolated image (im_os), but a smooth
    // (gradient-type) prior belongs to the primal step and must be evaluated at the primal iterate u (uCP).
    // The proximal/dual-split TV/TGV priors take their dual step at the extrapolated image and are excluded,
    // as is PDDY, which evaluates the prior at its own point x.
    const bool priorAtPrimalU = (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1) && !MethodList.PDDY
        && inputScalars.subsetsUsed == 1 && !MethodList.ProxTV && !MethodList.ProxTGV && !MethodList.TGV && kk == 0;

    // Subset scaling of the prior strength. These methods combine one subset's data gradient (about 1/S of the full gradient)
    // with the prior on every subiteration, so beta/S per subiteration gives the objective L + beta*R independently of the
    // number of subsets. PDHG (the dual sums over the subsets), OSL-COSEM (full sensitivity), BSREM/ROSEM-MAP (prior once per
    // iteration), SART and ASD-POCS (normalised prior steps) are not scaled
    const float subsetScale = (MethodList.OSLOSEM || MethodList.RBIOSL || MethodList.PKMA || MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS || MethodList.SAGA)
        ? 1.f / static_cast<float>(inputScalars.subsetsUsed) : 1.f;

    // The filter-phase switch is done once, for every timestep and volume, before the loops
    if (w_vec.precondTypeIm[5] && w_vec.filterIter > 0 && osa_iter + inputScalars.subsets * iter  == w_vec.filterIter) {
        if (inputScalars.verbose >= 3)
            mexPrint("Image-based filter iterations complete.");
        if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
            for (uint32_t t = 0; t < inputScalars.Nt; t++) {
                for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
                    if (w_vec.sigmaCP[t][ii] == 1.f)
                        w_vec.tauCP[t][ii] = w_vec.tauCP2[t][ii];
                    else if (w_vec.sigmaCP[t][ii] == w_vec.tauCP[t][ii]) {
                        w_vec.tauCP[t][ii] = w_vec.tauCP2[t][ii];
                        w_vec.sigmaCP[t][ii] = w_vec.tauCP2[t][ii];
                    }
                    else {
                        w_vec.sigmaCP[t][ii] = w_vec.tauCP2[t][ii];
                    }
                }
            }
        }
        if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA) {
            for (uint32_t t = 0; t < inputScalars.Nt; t++)
                w_vec.lambda[t] = w_vec.lambdaFiltered[t];
        }
        w_vec.precondTypeIm[5] = false;
    }

    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        for (int ii = kk; ii <= inputScalars.nMultiVolumes; ii++) {
            // For list-mode data the -1 term of the Poisson gradient is replaced by the sensitivity image D (subset sensitivity):
            // the backprojection contains only BP(y / ybar) over the events
            if (MethodList.PKMA && inputScalars.listmode > 0)
                vec.rhs_os[timestep][ii] = w_vec.D[DIndex(inputScalars, timestep)][ii] - vec.rhs_os[timestep][ii];
            if ((MethodList.RAMLA || MethodList.MRAMLA || MethodList.BSREM || MethodList.MBSREM || MethodList.SPS || MethodList.SAGA || MethodList.DRAMA
                || MethodList.RBI || MethodList.RBIOSL) && inputScalars.listmode > 0)
                vec.rhs_os[timestep][ii] -= w_vec.D[DIndex(inputScalars, timestep)][ii];
            if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.PDDY)
                PDHG1(vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, timestep, osa_iter + inputScalars.subsets * iter, ii);
            // PDDY (Salim, Condat, Mishchenko and Richtarik, JOTA 2022), one subset. p = uCP is the auxiliary primal variable and
            // the image estimate s = im_os was used for the forward projection (dual update). Here rhs_os = A^T y+ (PDHG1 does
            // nothing without subsets). The point x = p - gamma * A^T y+ is where the prior gradient is evaluated, so it is placed
            // into im_os. PDHG2 then forms s+ = prox(x - gamma * (A^T y+ + grad R(x))) and p+ = p + s+ - x
            if (MethodList.PDDY)
                vec.im_os[timestep][ii] = vec.uCP[timestep][ii] - w_vec.tauCP[timestep][ii] * vec.rhs_os[timestep][ii];
            if (status != 0)
                return -1;
            vec.im_os[timestep][ii].eval();
            vec.rhs_os[timestep][ii].eval();
        }

        if (DEBUG) {
            mexPrint("Priori\n");
            mexEval();
        }
        // The spatial prior gradient is recomputed only every regEveryIter-th (sub)iteration and the last computed gradient
        // (vec.dU) is used in between, for all algorithms. The counter below runs from 0 to
        // inputScalars.subsets * inputScalars.Niter - 1; the first and last are always computed.
        const int64_t regCounter = static_cast<int64_t>(osa_iter) + static_cast<int64_t>(inputScalars.subsets) * static_cast<int64_t>(iter);
        const int64_t regCounterMax = static_cast<int64_t>(inputScalars.subsets) * static_cast<int64_t>(inputScalars.Niter) - 1;
        const bool computeReg = inputScalars.regEveryIter <= 1 || regCounter == 0 || regCounter == regCounterMax
            || (regCounter % static_cast<int64_t>(inputScalars.regEveryIter)) == 0;
        if (!MethodList.BSREM && !MethodList.ROSEMMAP && !MethodList.POCS && !MethodList.SART && kk == 0 && computeReg) {
            const bool swapped = priorAtPrimalU && swapPDHGPrimal(vec, timestep);
            // The image-based priors read the forward projection image cache, so give them u (or the PDDY point x) as well
            if (swapped || MethodList.PDDY)
                status = refreshPriorImageCache(vec.im_os[timestep][0], inputScalars, proj, 0);
            // The prior gradient (positive beta, scaled with the subset scale) is always computed into the separate buffer vec.dU.
            // The gradient-type methods add it into their rhs below, with the sign of their update direction
            if (status == 0)
                status = applySpatialPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta * subsetScale, timestep, osa_iter + inputScalars.subsets * iter, compute_norm_matrix, false, vv);
            if (swapped)
                std::swap(vec.im_os[timestep][0], vec.uCP[timestep][0]);
            if (status != 0) return -1;
        }
    }

    if (MethodList.TemporalSmoothness || MethodList.TemporalTV) {
        if (!MethodList.BSREM && !MethodList.ROSEMMAP && !MethodList.POCS && !MethodList.SART && kk == 0) {
            std::vector<uint8_t> swappedT(inputScalars.Nt, 0);
            if (priorAtPrimalU)
                for (uint32_t t = 0; t < inputScalars.Nt; t++)
                    swappedT[t] = swapPDHGPrimal(vec, t) ? 1 : 0;
            status = applyTemporalPrior(vec, w_vec, MethodList, inputScalars, proj, false, nullptr, -1, subsetScale);
            for (uint32_t t = 0; t < inputScalars.Nt; t++)
                if (swappedT[t] != 0)
                    std::swap(vec.im_os[t][0], vec.uCP[t][0]);
            if (status != 0) return -1;
        }
    }

    // Add the prior gradient (spatial plus temporal) into the gradient of the methods that take the prior in their rhs. The
    // gradient is the one computed in this subiteration or, if the spatial prior was skipped because of regEveryIter, the last
    // computed one. MBSREM, MRAMLA and SPS use an ascent-direction rhs (update is im + lam * rhs), so the penalty gradient is
    // subtracted for them. The one-step-late methods (OSL-OSEM, RBI-OSL, OSL-COSEM, SAGA) and BSREM/ROSEM-MAP read the stored
    // gradient themselves. The prior is applied to the main volume only
    if (MAP && kk == 0 && (MethodList.PKMA || MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS || MethodList.PDHG || MethodList.PDHGKL
        || MethodList.PDHGL1 || MethodList.PDDY || MethodList.FISTA || MethodList.FISTAL1)) {
        const float priorSign = (MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS) ? -1.f : 1.f;
        for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
            const af::array priorGrad = totalPriorGradient(vec, MethodList, timestep);
            if (priorGrad.elements() > 1) {
                vec.rhs_os[timestep][0] += priorSign * priorGrad;
                vec.rhs_os[timestep][0].eval();
            }
        }
    }

    // Snapshot of the frames for the SART temporal prior: the neighbour frames must be the ones from the start of this
    // subiteration, not frames that were already updated in the timestep loop below. ArrayFire arrays are copy-on-write and
    // the SART updates assign new arrays, so a shallow copy keeps the old values
    std::vector<af::array> sartFrames;
    if (MethodList.SART && (MethodList.TemporalSmoothness || MethodList.TemporalTV) && inputScalars.Nt > 1) {
        sartFrames.resize(inputScalars.Nt);
        for (uint32_t t = 0; t < inputScalars.Nt; t++)
            sartFrames[t] = vec.im_os[t][0];
    }
    
    // The COSEM family divides the sum of the complete data images by the full sensitivity. w_vec.D is the sensitivity
    // divided by the number of subsets, so it is multiplied back by nSubD in the COSEM calls (Hsiao's COSEM)
    const float nSubD = static_cast<float>(inputScalars.subsetsUsed);
    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        if (largeDim && kk == 0 && vec.im_os[timestep][0].elements() > vec.rhs_os[timestep][0].elements()) {
            if (vv == 0)
                vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(0, vec.rhs_os[timestep][0].elements() - 1));
            else if (vv < inputScalars.subsetsUsed - 1)
                    vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(inputScalars.lDimStruct.startPr[vv], vec.rhs_os[timestep][0].elements() - 1 + inputScalars.lDimStruct.endPr[vv]));
            else
                    vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(inputScalars.lDimStruct.startPr[vv], af::end));
        }

        // Sensitivity image slot for this timestep: 0 for a shared sensitivity image (the common case,
        // and always for Nt == 1), or the real timestep when inputScalars.sensPerTimestep requires a
        // separate slot per timestep (see reconstructionAF.h).
        const uint32_t sIdx = inputScalars.sensPerTimestep ? timestep : 0u;
        // Spatial plus temporal prior gradient for the one-step-late type methods
        const af::array dUtotal = (MethodList.OSLOSEM || MethodList.RBIOSL || MethodList.OSLCOSEM > 0u) ? totalPriorGradient(vec, MethodList, timestep) : af::array();
        for (int ii = kk; ii <= inputScalars.nMultiVolumes; ii++) {
            af::array* Sens = nullptr;
            // Cache of max(Summ / D) for the plain RBI update. Only valid when the per-subset sensitivity images
            // (compute_norm_matrix == 2) are final, i.e. after the first iteration; with a single sensitivity image
            // (compute_norm_matrix == 1) it is recomputed for every subset and the maximum is not constant.
            float* rbiMaxCache = nullptr;
            if ((MethodList.RBI || MethodList.RBIOSL) && iter > 0 && compute_norm_matrix == 2u && !largeDim) {
                const size_t cacheSize = static_cast<size_t>(inputScalars.Nt) * (inputScalars.nMultiVolumes + 1) * inputScalars.subsetsUsed;
                if (w_vec.RBImaxCache.size() != cacheSize)
                    w_vec.RBImaxCache.assign(cacheSize, 0.f);
                // The key includes the timestep, since D can differ per timestep
                rbiMaxCache = &w_vec.RBImaxCache[(static_cast<size_t>(timestep) * (inputScalars.nMultiVolumes + 1) + ii) * inputScalars.subsetsUsed + osa_iter];
            }
            if (compute_norm_matrix == 1u) {
                Sens = &vec.Summ[sIdx][ii][0];
            } else if (compute_norm_matrix == 2u) {
                Sens = &vec.Summ[sIdx][ii][osa_iter];
            }

            // Compute the (matrix free) algorithms
            // Ordered Subsets Expectation Maximization (OSEM)
            if (MethodList.OSEM || MethodList.ECOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing OSEM/ECOSEM");
                if (MethodList.ECOSEM)
                    OSEMApu = EM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii] + inputScalars.epps);
                else
                    if (inputScalars.CT)
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii]);
                    else
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii] + inputScalars.epps);
            }

            // Modfied Row-action Maximum Likelihood (MRAMLA)
            if (MethodList.MRAMLA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing MRAMLA");
                status = MBSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda[timestep], timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
                if (status != 0)
                    return -1;
            }

            // Row-action Maximum Likelihood (RAMLA)
            if (MethodList.RAMLA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RAMLA");
                status = BSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter, inputScalars, proj, ii);
                if (status != 0)
                    return -1;
            }

            // Relaxed OSEM (ROSEM)
            if (MethodList.ROSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ROSEM");
                if (inputScalars.CT)
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter);
                else
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter);
            }

            // Rescaled Block Iterative EM (RBI)
            if (MethodList.RBI) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RBI");
                vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[DIndex(inputScalars, timestep)][ii], 0.f, af::constant(0.f, 1, 1), rbiMaxCache);
            }

            // Dynamic RAMLA
            if (MethodList.DRAMA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing DRAMA");
                vec.im_os[timestep][ii] = DRAMA(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter, osa_iter, inputScalars.Niter);
            }

            // Complete data OSEM
            if (MethodList.COSEM || MethodList.ECOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing COSEM/ECOSEM");
                updateCOSEMCompleteData(vec, ii, osa_iter, vec.rhs_os[timestep][ii] * vec.im_os[timestep][ii]);
                if (MethodList.ECOSEM)
                    COSEMApu = COSEM(vec.im_os[timestep][ii], vec.C_coSum[ii], nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii], w_vec.h_ACOSEM, 2u);
                else
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_coSum[ii], nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii], w_vec.h_ACOSEM, 2u);
            }

            // Enhanced COSEM
            if (MethodList.ECOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ECOSEM");
                vec.im_os[timestep][ii] = ECOSEM(vec.im_os[timestep][ii], w_vec.D[DIndex(inputScalars, timestep)][ii], OSEMApu, COSEMApu, inputScalars.epps);
            }

            // Accelerated COSEM
            if (MethodList.ACOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ACOSEM");
                if (DEBUG) {
                    mexPrintBase("vec.rhs_os[timestep][ii].dims(0) = %u\n", vec.rhs_os[timestep][ii].dims(0));
                    mexPrintBase("vec.im_os[timestep][ii].dims(0) = %u\n", vec.im_os[timestep][ii].dims(0));
                    mexPrintBase("vec.C_co[ii].dims(0) = %u\n", vec.C_co[ii].dims(0));
                    mexPrintBase("vec.C_co[ii].dims(1) = %u\n", vec.C_co[ii].dims(1));
                    mexPrintBase("osa_iter = %u\n", osa_iter);
                    mexEval();
                }
                updateCOSEMCompleteData(vec, ii, osa_iter, vec.rhs_os[timestep][ii] * af::pow(vec.im_os[timestep][ii], w_vec.h_ACOSEM_2));
                if (DEBUG) {
                    mexPrintBase("im_os[timestep] = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
                    mexPrintBase("rhs_os[timestep] = %f\n", af::sum<float>(vec.rhs_os[timestep][ii]));
                    mexPrintBase("D = %f\n", af::sum<float>(w_vec.D[DIndex(inputScalars, timestep)][ii]));
                    mexPrintBase("C_co / D = %f\n", af::sum<float>(af::sum(vec.C_co[ii], 1) / w_vec.D[DIndex(inputScalars, timestep)][ii]));
                    mexPrintBase("C_co = %f\n", af::sum<float>(vec.C_co[ii], 1));
                    mexPrintBase("C_co(:,osa_iter) = %f\n", af::sum<float>(vec.C_co[ii](af::span, osa_iter), 1));
                    mexPrintBase("min(D) = %f\n", af::min<float>(w_vec.D[DIndex(inputScalars, timestep)][ii]));
                    mexPrintBase("h = %f\n", w_vec.h_ACOSEM);
                    mexPrintBase("h2 = %f\n", w_vec.h_ACOSEM_2);
                    mexEval();
                }
                float uu;
                vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_coSum[ii], nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii], w_vec.h_ACOSEM, 1u);
                status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, timestep, mData[timestep], m_size, w_vec, vec, proj, subSum, pituus, g);
                if (status != 0)
                    return -1;
                if (inputScalars.CT)
                    vec.im_os[timestep][ii] *= (w_vec.ACOSEM_rhs / uu);
                else
                    vec.im_os[timestep][ii] *= (uu / w_vec.ACOSEM_rhs);
                vec.im_os[timestep][ii].eval();
                if (DEBUG) {
                    mexPrintBase("w_vec.ACOSEM_rhs = %f\n", w_vec.ACOSEM_rhs);
                    mexPrintBase("uu = %f\n", uu);
                    mexPrintBase("uu / w_vec.ACOSEM_rhs = %f\n", uu / w_vec.ACOSEM_rhs);
                    mexEval();
                }
            }
            if (MethodList.OSLOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing OSL-OSEM");
                // The prior is computed only for the main volume, the extra (multi-resolution) volumes use the OSEM update.
                // The denominator (sensitivity image + prior gradient) is limited to at least 10 % of the sensitivity image (and at least
                // epps where the sensitivity is ~0, i.e. outside the measured region, where a purely relative floor would allow a 10x
                // growth per subset), so that the update cannot change sign or divide by zero with a strong prior and a single update
                // can be at most 10x the OSEM update. Without a prior this is identical to OSEM
                if (inputScalars.CT) {
                    if (ii == 0)
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], (af::max)(*Sens + dUtotal, (af::max)(0.1f * *Sens, static_cast<double>(inputScalars.epps))), inputScalars.flat * vec.rhs_os[timestep][ii]);
                    else
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii]);
                }
                else {
                    if (ii == 0)
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], (af::max)(*Sens + dUtotal, (af::max)(0.1f * *Sens, static_cast<double>(inputScalars.epps))), vec.rhs_os[timestep][ii] + inputScalars.epps);
                    else
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii] + inputScalars.epps);
                }
            }
            else if (MethodList.BSREM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing BSREM");
                status = BSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter, inputScalars, proj, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.MBSREM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing MBSREM");
                status = MBSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda[timestep], timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.ROSEMMAP) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ROSEMMAP");
                if (inputScalars.CT)
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter);
                else
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda[timestep], iter);
            }
            else if (MethodList.RBIOSL) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RBIOSL");
                if (ii == 0)
                    vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[DIndex(inputScalars, timestep)][ii], (w_vec.beta != 0.f || (MethodList.TemporalSmoothness || MethodList.TemporalTV)) ? 1.f : 0.f, dUtotal, rbiMaxCache, inputScalars.epps);
                else
                    vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[DIndex(inputScalars, timestep)][ii], 0.f, af::constant(0.f, 1, 1), rbiMaxCache);
            }
            else if (MethodList.OSLCOSEM > 0u) {
                if (inputScalars.verbose >= 3)
                    mexPrintVar("Computing OSLCOSEM ", MethodList.OSLCOSEM);
                if (MethodList.OSLCOSEM == 1u)
                    updateCOSEMCompleteData(vec, ii, osa_iter, vec.rhs_os[timestep][ii] * pow(vec.im_os[timestep][ii], w_vec.h_ACOSEM_2));
                else
                    updateCOSEMCompleteData(vec, ii, osa_iter, vec.rhs_os[timestep][ii] * vec.im_os[timestep][ii]);
                if (ii == 0)
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_coSum[ii], (af::max)(nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii] + dUtotal, (af::max)(0.1f * nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii], static_cast<double>(inputScalars.epps))), w_vec.h_ACOSEM, MethodList.OSLCOSEM);
                else
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_coSum[ii], nSubD * w_vec.D[DIndex(inputScalars, timestep)][ii], w_vec.h_ACOSEM, MethodList.OSLCOSEM);
                if (MethodList.OSLCOSEM == 1u) {
                    float uu = 0.f;
                    status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, timestep, mData[timestep], m_size, w_vec, vec, proj, subSum, pituus, g);
                    if (status != 0)
                        return -1;
                    if (DEBUG) {
                        mexPrintBase("w_vec.ACOSEM_rhs1 = %f\n", w_vec.ACOSEM_rhs);
                        mexPrintBase("uu = %f\n", uu);
                        mexEval();
                    }
                    if (inputScalars.CT)
                        vec.im_os[timestep][ii] = vec.im_os[timestep][ii] * (w_vec.ACOSEM_rhs / uu);
                    else
                        vec.im_os[timestep][ii] = vec.im_os[timestep][ii] * (uu / w_vec.ACOSEM_rhs);
                    if (DEBUG) {
                        mexPrintBase("w_vec.ACOSEM_rhs3 = %f\n", w_vec.ACOSEM_rhs);
                        mexEval();
                    }
                }
            }
            else if (MethodList.PKMA) {
                if (DEBUG || inputScalars.verbose >= 3)
                    mexPrint("Computing PKMA");
                status = PKMA(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec, inputScalars, timestep, iter, osa_iter, proj, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.SPS) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing SPS");
                status = SPS(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda[timestep], timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.LSQR) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing LSQR");
                LSQR(inputScalars, w_vec, timestep, iter, vec, ii);
            }
            else if (MethodList.CGLS) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing CGLS");
                CGLS(inputScalars, w_vec, timestep, iter, vec, ii, inputScalars.largeDim);
            }
            else if (MethodList.SART || MethodList.POCS) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing SART");
                if (DEBUG) {
                    mexPrintBase("w_vec.lambda[timestep][iter] = %f\n", w_vec.lambda[timestep][iter]);
                    mexEval();
                }
                status = SART(inputScalars, w_vec, MethodList, vec, proj, mData[timestep], g, length, pituus, timestep, osa_iter, iter , *Sens, vec.rhs_os[timestep][ii], w_vec.lambda[timestep][iter], ii, sartFrames.empty() ? nullptr : &sartFrames);
                if (status != 0)
                    return -1;
                if (MethodList.POCS) {
                    status = POCS(inputScalars, w_vec, MethodList, vec, proj, mData[timestep], g, length, pituus, timestep, osa_iter, iter, ii);
                    if (status != 0)
                        return -1;
                }
            }
            else if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.PDDY) {
                if (DEBUG || inputScalars.verbose >= 3)
                    mexPrint("Computing PDHG/PDHGKL/PDHGL1/PDDY");
                status = PDHG2(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, proj, iter, timestep, osa_iter, ii, pituus, g, m_size, length, MethodList.PDDY);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.FISTA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing FISTA");
                status = FISTA(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, proj, timestep, iter, osa_iter, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.FISTAL1) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing FISTAL1");
                status = FISTAL1(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, w_vec.beta, proj, timestep, iter, osa_iter, ii);
                if (status != 0)
                    return -1;
            }
            else if (MethodList.SAGA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing SAGA");
                status = SAGA(vec.im_os[timestep][ii], inputScalars, w_vec, MethodList, vec, proj, timestep, osa_iter, iter, ii);
                if (status != 0)
                    return -1;
            } else if (MethodList.BB) {
                if(inputScalars.verbose>=3)
                    mexPrint("Computing BB");
                status = BB(vec, w_vec, timestep, ii);
                if (status != 0)
                    return -1;
            }
            if (inputScalars.FISTAAcceleration) {
                // Nesterov: x_{k+1} = update(y_k), t_{k+1} = (1 + sqrt(1 + 4 t_k^2)) / 2,
                // y_{k+1} = x_{k+1} + (t_k - 1) / t_{k+1} * (x_{k+1} - x_k). vec.xFISTAAcc holds the previous un-extrapolated iterate x_k
                const float tOld = w_vec.tFISTA[timestep];
                const float tNew = (1.f + std::sqrt(1.f + 4.f * tOld * tOld)) / 2.f;
                const af::array xNew = vec.im_os[timestep][ii].copy();
                if (DEBUG) {
                    mexPrintBase("im_os1 = %f\n", af::sum<float>(xNew));
                    mexEval();
                }
                if (vec.xFISTAAcc[timestep].size() < static_cast<size_t>(inputScalars.nMultiVolumes + 1))
                    vec.xFISTAAcc[timestep].resize(inputScalars.nMultiVolumes + 1);
                if (vec.xFISTAAcc[timestep][ii].isempty()) {
                    // First call for this volume: no extrapolation
                    vec.xFISTAAcc[timestep][ii] = xNew;
                }
                else {
                    vec.im_os[timestep][ii] = xNew + (tOld - 1.f) / tNew * (xNew - vec.xFISTAAcc[timestep][ii]);
                    vec.xFISTAAcc[timestep][ii] = xNew;
                }
                af::eval(vec.im_os[timestep][ii]);
                if (ii == inputScalars.nMultiVolumes)
                    w_vec.tFISTA[timestep] = tNew;
                if (DEBUG) {
                    mexPrintBase("im_os[timestep] = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
                    mexPrintBase("xFISTAAcc = %f\n", af::sum<float>(vec.xFISTAAcc[timestep][ii]));
                    mexPrintBase("t = %f\n", tOld);
                    mexPrintBase("w_vec.tFISTA[timestep] = %f\n", w_vec.tFISTA[timestep]);
                    mexEval();
                }
            }
            vec.im_os[timestep][ii].eval();
        }
        
        if (DEBUG || inputScalars.verbose >= 3) {
            af::sync();
            proj.tEndLocal = std::chrono::steady_clock::now();
            const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
            mexPrintBase("Iterative algorithm computed in %f seconds\n", tDiff.count());
        }
        
    }
    return status;
}