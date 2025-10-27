// Utilities to build preprocessor macros for different Metal kernels.
// TODO: use this code also for OpenCL/CUDA preprocessor directives through type definitions.

#pragma once
#include "Metal.hpp"
#include "structs.h"

#define TH 100000000000.f
#define TH32 100000.f
#define NVOXELS 8
#define NVOXELS5 1
#define NVOXELSFP 8

static NS::String* S(const std::string& s) {
    return NS::String::string(s.c_str(), NS::UTF8StringEncoding);
}

static NS::String* KS(const char* s) {
    return NS::String::string(s, NS::ASCIIStringEncoding);
}

static inline NS::Number* KTrue() { return NS::Number::number(1); }
static inline NS::Number* KNum(uint64_t v) { return NS::Number::number((unsigned long long)v); }

struct MacroBuilder {
    std::vector<NS::Object*> keys;
    std::vector<NS::Object*> vals;

    void addFlag(const char* name) {
        keys.push_back(KS(name));
        vals.push_back(KTrue());
    }
    void addInt(const char* name, uint64_t v) {
        keys.push_back(KS(name));
        vals.push_back(KNum(v));
    }
    MacroBuilder clone() const {
        MacroBuilder mb;
        mb.keys = keys;
        mb.vals = vals;
        return mb;
    }
    NS::Dictionary* build() const {
        return NS::Dictionary::dictionary(
            const_cast<NS::Object* const*>(vals.data()),
            const_cast<NS::Object* const*>(keys.data()),
            (NS::UInteger)keys.size()
        )->autorelease();
    }
};

struct MacroSets {
    NS::Dictionary* base = nullptr;
    NS::Dictionary* fp = nullptr;
    NS::Dictionary* bp = nullptr;
    NS::Dictionary* sens = nullptr;
    NS::Dictionary* aux = nullptr;
};

static MacroSets BuildMacroDict(
    const scalarStruct& inputScalars,
    const Weighting& w_vec,
    const RecMethods MethodList,
    const int type
) {
    MacroSets out;
    
    // Common preprocessor macros
    MacroBuilder base;
    base.addFlag("METAL");

    if (inputScalars.useHalf) base.addFlag("HALF");
    if (inputScalars.useParallelBeam) base.addFlag("PARALLEL");
    if (inputScalars.raw == 1) base.addFlag("RAW");
    if (inputScalars.useTotLength && !inputScalars.SPECT) base.addFlag("TOTLENGTH");

    if (inputScalars.maskFP) {
        base.addFlag("MASKFP");
        if (inputScalars.maskFPZ > 1) base.addFlag("MASKFP3D");
    }
    if (inputScalars.maskBP) {
        base.addFlag("MASKBP");
        if (inputScalars.maskBPZ > 1) base.addFlag("MASKBP3D");
    }

    if (inputScalars.offset) base.addFlag("OFFSET");
    if (inputScalars.attenuation_correction == 1u && inputScalars.CTAttenuation) base.addFlag("ATN");
    else if (inputScalars.attenuation_correction == 1u && !inputScalars.CTAttenuation) base.addFlag("ATNM");
    if (inputScalars.normalization_correction == 1u) base.addFlag("NORM");
    if (inputScalars.scatter == 1u) base.addFlag("SCATTER");
    if (inputScalars.randoms_correction == 1u) base.addFlag("RANDOMS");

    if (inputScalars.nLayers > 1U) {
        if (inputScalars.listmode > 0 && inputScalars.indexBased)
            base.addInt("NLAYERS", inputScalars.nLayers);
        else
            base.addInt("NLAYERS", inputScalars.nProjections / (inputScalars.nLayers * inputScalars.nLayers));
    }

    if (inputScalars.TOF) base.addFlag("TOF");
    if (inputScalars.CT) base.addFlag("CT");
    else if (inputScalars.PET) base.addFlag("PET");
    else if (inputScalars.SPECT) base.addFlag("SPECT");

    base.addInt("NBINS", inputScalars.nBins);

    if (inputScalars.listmode == 1) base.addFlag("LISTMODE");
    else if (inputScalars.listmode == 2) base.addFlag("LISTMODE2");
    if (inputScalars.listmode > 0 && inputScalars.indexBased) base.addFlag("INDEXBASED");

    const bool siddonVal = (inputScalars.FPType == 1 || inputScalars.BPType == 1 || inputScalars.FPType == 4 || inputScalars.BPType == 4);
    if ((siddonVal && ((inputScalars.n_rays * inputScalars.n_rays3D) > 1)) || inputScalars.SPECT) {
        base.addInt("N_RAYS",  uint64_t(inputScalars.n_rays) * uint64_t(inputScalars.n_rays3D));
        base.addInt("N_RAYS2D", inputScalars.n_rays);
        base.addInt("N_RAYS3D", inputScalars.n_rays3D);
    }

    if (inputScalars.pitch) base.addFlag("PITCH");
    if (((inputScalars.subsets > 1 &&
        (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) &&
        !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
        base.addFlag("SUBSETS");

    if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
        base.addInt("STYPE",    inputScalars.subsetType);
        base.addInt("NSUBSETS", inputScalars.subsets);
    }

    if (inputScalars.FPType == 2 || inputScalars.BPType == 2 ||
        inputScalars.FPType == 3 || inputScalars.BPType == 3)
    {
        if (inputScalars.orthXY) base.addFlag("CRYSTXY");
        if (inputScalars.orthZ) base.addFlag("CRYSTZ");
    }

    // FP
    auto fp = base.clone();
    fp.addFlag("FP");
    if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        fp.addFlag("SIDDON");
        fp.addFlag("ATOMICF");
        if (inputScalars.FPType == 2 || inputScalars.FPType == 3) fp.addFlag("ORTH");
        if (inputScalars.FPType == 3) fp.addFlag("VOL");
    } else if (inputScalars.FPType == 4) { 
        fp.addFlag("PTYPE4");
        if (!inputScalars.largeDim) fp.addInt("NVOXELS", NVOXELS);
    } else if (inputScalars.FPType == 5) { 
        fp.addFlag("PROJ5");
        if (inputScalars.meanFP) fp.addFlag("MEANDISTANCEFP");
        if (inputScalars.pitch) fp.addInt("NVOXELS5", 1);
        else fp.addInt("NVOXELS5", NVOXELS5);
        fp.addInt("NVOXELSFP", NVOXELSFP);
    }
    out.fp = fp.build();

    // BP
    auto bp = base.clone();
    bp.addFlag("BP");
    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
        bp.addFlag("SIDDON");
        bp.addFlag("ATOMICF");
        if (inputScalars.BPType == 2 || inputScalars.BPType == 3) bp.addFlag("ORTH");
        if (inputScalars.BPType == 3) bp.addFlag("VOL");
    } else if (inputScalars.BPType == 4) { 
        bp.addFlag("PTYPE4");
        if (!inputScalars.largeDim) bp.addInt("NVOXELS", NVOXELS);
    } else if (inputScalars.BPType == 5) { 
        bp.addFlag("PROJ5");
        if (inputScalars.meanBP) bp.addFlag("MEANDISTANCEBP");
        if (inputScalars.pitch) bp.addInt("NVOXELS5", 1);
        else bp.addInt("NVOXELS5", NVOXELS5);
        bp.addInt("NVOXELSFP", NVOXELSFP);
    }
    out.bp = bp.build();

    // Sensitivity image
    auto sens = base.clone();
    sens.addFlag("BP");
    sens.addFlag("SENS");
    if (inputScalars.BPType == 3) sens.addFlag("VOL");
    if (inputScalars.BPType == 2 || inputScalars.BPType == 3) sens.addFlag("ORTH");
    if (inputScalars.BPType == 4) {
        sens.addFlag("PTYPE4");
        sens.addInt("NVOXELS", NVOXELS);
    } else {
        sens.addFlag("SIDDON");
    }
    out.sens = sens.build();

    // Auxiliary kernel
    MacroBuilder aux;
    aux.addFlag("METAL");
    if (inputScalars.largeDim) aux.addFlag("LARGEDIM");
    if (inputScalars.useExtendedFOV) aux.addFlag("EFOV");
    if (type == 2) {
        if (inputScalars.use_psf) aux.addFlag("PSF");
    } else if (type == 0) {
        if (inputScalars.CT) aux.addFlag("CT");
        if (inputScalars.randoms_correction) aux.addFlag("RANDOMS");
        if (inputScalars.use_psf) aux.addFlag("PSF");
    } else {
        aux.addFlag("AF");
    }
    if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
        aux.addFlag("MASKPRIOR");
        if (inputScalars.maskBPZ > 1) aux.addFlag("MASKBP3D");
    }
    if (inputScalars.eFOV) aux.addFlag("EFOVZ");
    if (MethodList.MRP) {
        aux.addFlag("MEDIAN");
        aux.addInt("SEARCH_WINDOW_X", w_vec.Ndx);
        aux.addInt("SEARCH_WINDOW_Y", w_vec.Ndy);
        aux.addInt("SEARCH_WINDOW_Z", w_vec.Ndz);
    }
    if (MethodList.NLM) {
        aux.addFlag("NLM_");
        if (w_vec.NLTV) aux.addInt("NLTYPE", 1);
        else if (w_vec.NLM_MRP) aux.addInt("NLTYPE", 2);
        else if (w_vec.NLRD) aux.addInt("NLTYPE", 3);
        else if (w_vec.NLLange) aux.addInt("NLTYPE", 4);
        else if (w_vec.NLLangeFiltered) aux.addInt("NLTYPE", 5);
        else if (w_vec.NLGGMRF) aux.addInt("NLTYPE", 6);
        else aux.addInt("NLTYPE", 0);
        if (w_vec.NLAdaptive) aux.addFlag("NLMADAPTIVE");
        if (w_vec.NLM_anatomical) aux.addFlag("NLMREF");
        aux.addInt("SWINDOWX", w_vec.Ndx);
        aux.addInt("SWINDOWY", w_vec.Ndy);
        aux.addInt("SWINDOWZ", w_vec.Ndz);
        aux.addInt("PWINDOWX", w_vec.Nlx);
        aux.addInt("PWINDOWY", w_vec.Nly);
        aux.addInt("PWINDOWZ", w_vec.Nlz);
    }
    if (MethodList.GGMRF) {
        aux.addFlag("GGMRF");
        aux.addInt("SWINDOWX", w_vec.Ndx);
        aux.addInt("SWINDOWY", w_vec.Ndy);
        aux.addInt("SWINDOWZ", w_vec.Ndz);
    }
    if (MethodList.hyperbolic) {
        aux.addFlag("HYPER");
        aux.addInt("SWINDOWX", w_vec.Ndx);
        aux.addInt("SWINDOWY", w_vec.Ndy);
        aux.addInt("SWINDOWZ", w_vec.Ndz);
    }
    if (MethodList.RDP) {
        aux.addFlag("RDP");
        if (w_vec.RDPLargeNeighbor) {
            aux.addFlag("RDPCORNERS");
            aux.addInt("SWINDOWX", w_vec.Ndx);
            aux.addInt("SWINDOWY", w_vec.Ndy);
            aux.addInt("SWINDOWZ", w_vec.Ndz);
        }
        if (w_vec.RDP_anatomical) aux.addFlag("RDPREF");
    }
    if (MethodList.ProxRDP && w_vec.RDPLargeNeighbor) aux.addFlag("RDPCORNERS");
    if (MethodList.TV && !w_vec.data.TV_use_anatomical) {
        aux.addFlag("TVGRAD");
        if (w_vec.data.TVtype == 6) aux.addFlag("TVW1");
        else if (w_vec.data.TVtype == 4) aux.addFlag("SATV");
        else if (w_vec.data.TVtype == 2) aux.addFlag("JPTV");
        if (w_vec.derivType > 0) aux.addInt("DIFFTYPE", w_vec.derivType);
    } else if ((MethodList.TV && w_vec.data.TV_use_anatomical) || MethodList.APLS) {
        aux.addFlag("TVGRAD");
        if (w_vec.data.TVtype == 1) aux.addFlag("ANATOMICAL1");
        else if (w_vec.data.TVtype == 2) aux.addFlag("ANATOMICAL2");
        else if (w_vec.data.TVtype == 5 || MethodList.APLS) aux.addFlag("ANATOMICAL3");
        if (w_vec.derivType > 0) aux.addInt("DIFFTYPE", w_vec.derivType);
    }
    if (MethodList.ProxTV) {
        aux.addFlag("PROXTV");
        if (w_vec.UseL2Ball) aux.addFlag("L2");
        if (w_vec.derivType > 0) aux.addInt("DIFFTYPE", w_vec.derivType);
    }
    if (MethodList.ProxTGV || MethodList.TGV) {
        aux.addFlag("PROXTV");
        aux.addFlag("PROXTGV");
        if (w_vec.UseL2Ball) aux.addFlag("L2");
        if (w_vec.derivType > 0) aux.addInt("DIFFTYPE", w_vec.derivType);
        if (!inputScalars.TGV2D) aux.addFlag("TGVZ");
    }
    if (MethodList.ProxRDP) aux.addFlag("PROXRDP");
    if (MethodList.ProxNLM) {
        aux.addFlag("PROXNLM");
        if (w_vec.NLTV) aux.addInt("NLTYPE", 1);
        else if (w_vec.NLM_MRP) aux.addInt("NLTYPE", 2);
        else if (w_vec.NLRD) aux.addInt("NLTYPE", 3);
        else if (w_vec.NLLange) aux.addInt("NLTYPE", 4);
        else if (w_vec.NLLangeFiltered) aux.addInt("NLTYPE", 5);
        else aux.addInt("NLTYPE", 0);
        if (w_vec.NLM_anatomical) aux.addFlag("NLMREF");
        aux.addInt("SWINDOWX", w_vec.Ndx);
        aux.addInt("SWINDOWY", w_vec.Ndy);
        aux.addInt("SWINDOWZ", w_vec.Ndz);
        aux.addInt("PWINDOWX", w_vec.Nlx);
        aux.addInt("PWINDOWY", w_vec.Nly);
        aux.addInt("PWINDOWZ", w_vec.Nlz);
    }
    if (MethodList.PKMA) aux.addFlag("PKMA");
    else if (MethodList.MBSREM || MethodList.MRAMLA) aux.addFlag("MBSREM");
    else if (MethodList.BSREM || MethodList.RAMLA) aux.addFlag("BSREM");
    else if (MethodList.CPType) {
        aux.addFlag("PDHG");
        if (inputScalars.subsets > 1) aux.addFlag("SUBSETS");
    }
    if (inputScalars.projector_type == 6) aux.addFlag("ROTATE");
    out.aux = aux.build();

    return out;
}