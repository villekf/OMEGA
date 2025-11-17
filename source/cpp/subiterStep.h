#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"

inline int computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t iter, uint32_t osa_iter, scalarStruct& inputScalars, std::vector<int64_t>& length, bool& break_iter, const int64_t* pituus, const af::array& g, ProjectorClass& proj, const std::vector<af::array> &mData, uint64_t m_size, const int64_t subSum, const uint8_t compute_norm_matrix, const int kk = 0, const bool largeDim = false, const int vv = 0) {
    int status = 0;
    af::array OSEMApu, COSEMApu, PDDYApu;
    std::vector<af::array> FISTAApu;
    if (DEBUG || inputScalars.verbose >= 3) {
        proj.tStartLocal = std::chrono::steady_clock::now();
    }

    if (DEBUG) {
        mexPrint("Algo start\n");
        mexEval();
    }

    bool MAP = false;
    if (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV || MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.ProxTGV || MethodList.ProxTV || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF || MethodList.QuadraticSmoothnessTemporal)
        MAP = true;

    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        for (int ii = kk; ii <= inputScalars.nMultiVolumes; ii++) {
            if (inputScalars.FISTAAcceleration)
                FISTAApu.emplace_back(vec.im_os[timestep][ii].copy());

            if (w_vec.precondTypeIm[5] && w_vec.filterIter > 0 && osa_iter + inputScalars.subsets * iter  == w_vec.filterIter) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Image-based filter iterations complete.");
                if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
                    if (w_vec.sigmaCP[ii] == 1.f)
                        w_vec.tauCP[timestep][ii] = w_vec.tauCP2[timestep][ii];
                    else if (w_vec.sigmaCP[ii] == w_vec.tauCP[timestep][ii]) {
                        w_vec.tauCP[timestep][ii] = w_vec.tauCP2[timestep][ii];
                        w_vec.sigmaCP[ii] = w_vec.tauCP2[timestep][ii];
                    }
                    else {
                        w_vec.sigmaCP[ii] = w_vec.tauCP2[timestep][ii];
                    }
                }
                if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA)
                    w_vec.lambda = w_vec.lambdaFiltered;
                w_vec.precondTypeIm[5] = false;
            }
            if (MethodList.PKMA && inputScalars.listmode > 0 && (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2]))
                vec.rhs_os[timestep][ii] = w_vec.D[ii] - vec.rhs_os[timestep][ii];
            if ((MethodList.RAMLA || MethodList.MRAMLA || MethodList.BSREM || MethodList.MBSREM) && inputScalars.listmode > 0 && (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2]))
                vec.rhs_os[timestep][ii] -= w_vec.D[ii];
            if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
                PDHG1(vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, timestep, osa_iter + inputScalars.subsets * iter, ii);
            if (MethodList.PDDY && ii == 0 && MAP) {
                PDDYApu = vec.im_os[timestep][0].copy();
                vec.im_os[timestep][0] -= w_vec.tauCP[timestep][0] * vec.uCP[timestep][0];
                if (inputScalars.verbose == 3) {
                    mexPrint("Computing PDDY step\n");
                    mexEval();
                }
            }
            if (status != 0)
                return -1;
            vec.im_os[timestep][ii].eval();
            vec.rhs_os[timestep][ii].eval();
        }

        if (DEBUG) {
            mexPrint("Priori\n");
            mexEval();
        }
        if (!MethodList.BSREM && !MethodList.ROSEMMAP && !MethodList.POCS && !MethodList.SART && kk == 0) {
            status = applySpatialPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta, osa_iter + inputScalars.subsets * iter, compute_norm_matrix, false, vv);
            if (status != 0) return -1;
        }
    }
    
    if (!MethodList.BSREM && !MethodList.ROSEMMAP && !MethodList.POCS && !MethodList.SART && kk == 0) {
        status = applyTemporalPrior(vec, w_vec, MethodList, inputScalars, proj);
        if (status != 0) return -1;
    }
    
    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        if (MethodList.PDDY && MAP && kk == 0)
            vec.im_os[timestep][0] = PDDYApu.copy();
        if (largeDim && kk == 0 && vec.im_os[timestep][0].elements() > vec.rhs_os[timestep][0].elements()) {
            if (vv == 0)
                vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(0, vec.rhs_os[timestep][0].elements() - 1));
            else if (vv < inputScalars.subsetsUsed - 1)
                    vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(inputScalars.lDimStruct.startPr[vv], vec.rhs_os[timestep][0].elements() - 1 + inputScalars.lDimStruct.endPr[vv]));
            else
                    vec.im_os[timestep][0] = vec.im_os[timestep][0](af::seq(inputScalars.lDimStruct.startPr[vv], af::end));
        }

        for (int ii = kk; ii <= inputScalars.nMultiVolumes; ii++) {
            af::array* Sens = nullptr;
            if (compute_norm_matrix == 1u) {
                Sens = &vec.Summ[ii][0];
            }
            else if (compute_norm_matrix == 2u) {
                Sens = &vec.Summ[ii][osa_iter];
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
                status = MBSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda, timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
            }

            // Row-action Maximum Likelihood (RAMLA)
            if (MethodList.RAMLA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RAMLA");
                status = BSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.lambda, iter, inputScalars, proj, ii);
            }

            // Relaxed OSEM (ROSEM)
            if (MethodList.ROSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ROSEM");
                if (inputScalars.CT)
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii], w_vec.lambda, iter);
                else
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda, iter);
            }

            // Rescaled Block Iterative EM (RBI)
            if (MethodList.RBI) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RBI");
                vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[ii]);
            }

            // Dynamic RAMLA
            if (MethodList.DRAMA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing DRAMA");
                vec.im_os[timestep][ii] = DRAMA(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda, iter, osa_iter, inputScalars.subsets);
            }

            // Complete data OSEM
            if (MethodList.COSEM || MethodList.ECOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing COSEM/ECOSEM");
                vec.C_co(af::span, osa_iter) = vec.rhs_os[timestep][ii] * vec.im_os[timestep][ii];
                if (MethodList.ECOSEM)
                    COSEMApu = COSEM(vec.im_os[timestep][ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 2u);
                else
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 2u);
            }

            // Enhanced COSEM
            if (MethodList.ECOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ECOSEM");
                vec.im_os[timestep][ii] = ECOSEM(vec.im_os[timestep][ii], w_vec.D[ii], OSEMApu, COSEMApu, inputScalars.epps);
            }

            // Accelerated COSEM
            if (MethodList.ACOSEM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ACOSEM");
                if (DEBUG) {
                    mexPrintBase("vec.rhs_os[timestep][ii].dims(0) = %u\n", vec.rhs_os[timestep][ii].dims(0));
                    mexPrintBase("vec.im_os[timestep][ii].dims(0) = %u\n", vec.im_os[timestep][ii].dims(0));
                    mexPrintBase("vec.C_co.dims(0) = %u\n", vec.C_co.dims(0));
                    mexPrintBase("vec.C_co.dims(1) = %u\n", vec.C_co.dims(1));
                    mexPrintBase("osa_iter = %u\n", osa_iter);
                    mexEval();
                }
                vec.C_co(af::span, osa_iter) = vec.rhs_os[timestep][ii] * af::pow(vec.im_os[timestep][ii], w_vec.h_ACOSEM_2);
                if (DEBUG) {
                    mexPrintBase("im_os[timestep] = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
                    mexPrintBase("rhs_os[timestep] = %f\n", af::sum<float>(vec.rhs_os[timestep][ii]));
                    mexPrintBase("D = %f\n", af::sum<float>(w_vec.D[ii]));
                    mexPrintBase("C_co / D = %f\n", af::sum<float>(af::sum(vec.C_co, 1) / w_vec.D[ii]));
                    mexPrintBase("C_co = %f\n", af::sum<float>(vec.C_co, 1));
                    mexPrintBase("C_co(:,osa_iter) = %f\n", af::sum<float>(vec.C_co(af::span, osa_iter), 1));
                    mexPrintBase("min(D) = %f\n", af::min<float>(w_vec.D[ii]));
                    mexPrintBase("h = %f\n", w_vec.h_ACOSEM);
                    mexPrintBase("h2 = %f\n", w_vec.h_ACOSEM_2);
                    mexEval();
                }
                float uu;
                vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 1u);
                status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, timestep, mData[timestep], m_size, w_vec, vec, proj, subSum, pituus, g);
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
                if (ii == 0) {
                    if (inputScalars.CT)
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens + vec.dU, inputScalars.flat * vec.rhs_os[timestep][ii]);
                    else
                        vec.im_os[timestep][ii] = EM(vec.im_os[timestep][ii], *Sens + vec.dU, vec.rhs_os[timestep][ii]);
                }
                else {

                }
            }
            else if (MethodList.BSREM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing BSREM");
                status = BSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.lambda, iter, inputScalars, proj, ii);
            }
            else if (MethodList.MBSREM) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing MBSREM");
                status = MBSREM(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda, timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
            }
            else if (MethodList.ROSEMMAP) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing ROSEMMAP");
                if (inputScalars.CT)
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, inputScalars.flat * vec.rhs_os[timestep][ii], w_vec.lambda, iter);
                else
                    vec.im_os[timestep][ii] = ROSEM(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.lambda, iter);
            }
            else if (MethodList.RBIOSL) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing RBIOSL");
                if (ii == 0)
                    vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[ii], w_vec.beta, vec.dU);
                else
                    vec.im_os[timestep][ii] = RBI(vec.im_os[timestep][ii], *Sens, vec.rhs_os[timestep][ii], w_vec.D[ii], 0.f);
            }
            else if (MethodList.OSLCOSEM > 0u) {
                if (inputScalars.verbose >= 3)
                    mexPrintVar("Computing OSLCOSEM ", MethodList.OSLCOSEM);
                if (MethodList.OSLCOSEM == 1u)
                    vec.C_co(af::span, osa_iter) = vec.rhs_os[timestep][ii] * pow(vec.im_os[timestep][ii], w_vec.h_ACOSEM_2);
                else
                    vec.C_co(af::span, osa_iter) = vec.rhs_os[timestep][ii] * vec.im_os[timestep][ii];
                if (ii == 0)
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_co, w_vec.D[ii] + vec.dU, w_vec.h_ACOSEM, MethodList.OSLCOSEM);
                else
                    vec.im_os[timestep][ii] = COSEM(vec.im_os[timestep][ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, MethodList.OSLCOSEM);
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
            }
            else if (MethodList.SPS) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing SPS");
                status = SPS(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], w_vec.U, w_vec.lambda, timestep, iter, osa_iter, inputScalars, w_vec, proj, ii);
            }
            else if (MethodList.LSQR) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing LSQR");
                LSQR(inputScalars, w_vec, timestep, iter, vec);
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
                    mexPrintBase("w_vec.lambda[iter] = %f\n", w_vec.lambda[iter]);
                    mexEval();
                }
                status = SART(inputScalars, w_vec, MethodList, vec, proj, mData[timestep], g, length, pituus, timestep, osa_iter, iter , *Sens, vec.rhs_os[timestep][ii], w_vec.lambda[iter], ii);
                if (MethodList.POCS)
                    status = POCS(inputScalars, w_vec, MethodList, vec, proj, mData[timestep], g, length, pituus, timestep, osa_iter, iter, ii);
            }
            else if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
                if (DEBUG || inputScalars.verbose >= 3)
                    mexPrint("Computing PDHG/PDHGKL/PDHGL1/PDDY");
                status = PDHG2(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, proj, iter, timestep, osa_iter, ii, pituus, g, m_size, length);
            }
            else if (MethodList.FISTA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing FISTA");
                status = FISTA(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, proj, timestep, iter, osa_iter, ii);
            }
            else if (MethodList.FISTAL1) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing FISTAL1");
                status = FISTAL1(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii], inputScalars, w_vec, vec, w_vec.beta, proj, timestep, iter, osa_iter, ii);
            }
            else if (MethodList.SAGA) {
                if (inputScalars.verbose >= 3)
                    mexPrint("Computing SAGA");
                status = SAGA(vec.im_os[timestep][ii], inputScalars, w_vec, vec, proj, timestep, osa_iter, iter, ii);
            } else if (MethodList.BB) {
                if(inputScalars.verbose>=3)
                    mexPrint("Computing BB");
                status = BB(vec, w_vec, timestep, ii);
            }
            if (inputScalars.FISTAAcceleration) {
                //if ((w_vec.precondTypeIm[5] && w_vec.filterIter > 0 && osa_iter + inputScalars.subsets * iter >= w_vec.filterIter) || !w_vec.precondTypeIm[5]) {
                //if (osa_iter == inputScalars.subsets - 1) {
                    const float t = w_vec.tFISTA;
                    if (DEBUG) {
                        mexPrintBase("im_os1 = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
                        mexEval();
                    }
                    if (osa_iter + iter > 0)
                        vec.im_os[timestep][ii] = vec.im_os[timestep][ii] + (t - 1.f) / w_vec.tFISTA * (vec.im_os[timestep][ii] - FISTAApu[ii]);
                    af::eval(vec.im_os[timestep][ii]);
                    if (ii == inputScalars.nMultiVolumes)
                        w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tFISTA * w_vec.tFISTA)) / 2.f;
                    if (DEBUG) {
                        mexPrintBase("im_os[timestep] = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
                        mexPrintBase("FISTAApu = %f\n", af::sum<float>(FISTAApu[ii]));
                        mexPrintBase("t = %f\n", t);
                        mexPrintBase("w_vec.tFISTA = %f\n", w_vec.tFISTA);
                        mexEval();
                    }
                //}
            }
            vec.im_os[timestep][ii].eval();
        }
        
        if (DEBUG || inputScalars.verbose >= 3) {
            af::sync();
            proj.tEndLocal = std::chrono::steady_clock::now();
            const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
            mexPrintBase("Iterative algorithm computed in %f seconds\n", tDiff);
        }
        
    }
    return status;
}