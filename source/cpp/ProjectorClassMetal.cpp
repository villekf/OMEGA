#include "ProjectorClassMetal.hpp"

int ProjectorClass::createProgram(
    NS::SharedPtr<MTL::Library>& libFP,
    NS::SharedPtr<MTL::Library>& libBP,
    NS::SharedPtr<MTL::Library>& libAux,
    NS::SharedPtr<MTL::Library>& libSens,
    const char* header_directory,
    scalarStruct& inputScalars,
    const RecMethods MethodList,
    const Weighting& w_vec,
    const size_t local_size[],
    const int type
) {
    //NS::AutoreleasePool pool;

    mtlDevice = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
    if (!mtlDevice) {
        mexPrint("No Metal device available");
        return -1;
    }
    std::string kernelFile = header_directory;
    std::string kernel_path, kernel_pathBP;
    std::string contentFP, contentBP;
    std::string contentAux;
    std::string options;
    
    // Load kernel parameter structs
    std::string kernelParamsFile = kernelFile.substr(0, kernelFile.size()-7) + "cpp/kernelParams.hpp";
    std::ifstream sourceHeader(kernelParamsFile);
    std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());

    // Load the header text file
    std::ifstream sourceHeader2(kernelFile + "general_opencl_functions.h");
    std::string contentHeader2((std::istreambuf_iterator<char>(sourceHeader2)), std::istreambuf_iterator<char>());
    contentHeader += contentHeader2;

    // Load orthogonal/volume of intersection headers if applicable
    if (inputScalars.FPType == 2 || inputScalars.BPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 3) {
        std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
        std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
        contentHeader += contentHeader3;
    }

    kernel_path = kernelFile;
    kernel_pathBP = kernelFile;
    if (inputScalars.FPType > 0 && inputScalars.FPType != 6) {
        if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
            kernel_path += "projectorType123.cl";
        }
        else if (inputScalars.FPType == 4)
            kernel_path += "projectorType4.cl";
        else if (inputScalars.FPType == 5)
            kernel_path += "projectorType5.cl";
        std::ifstream sourceFile(kernel_path.c_str());
        std::string contentFFP((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
        contentFP = contentHeader + contentFFP;
    }
    if (inputScalars.BPType > 0 && inputScalars.BPType != 6) {
        if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
            kernel_pathBP += "projectorType123.cl";
        }
        else if (inputScalars.BPType == 4)
            kernel_pathBP += "projectorType4.cl";
        else if (inputScalars.BPType == 5)
            kernel_pathBP += "projectorType5.cl";
        std::ifstream sourceFileBP(kernel_pathBP.c_str());
        std::string contentFBP((std::istreambuf_iterator<char>(sourceFileBP)), std::istreambuf_iterator<char>());
        contentBP = contentHeader + contentFBP;
    }
    
    // Macros TODO: sens, aux
    auto sets = BuildMacroDict(inputScalars, w_vec, MethodList, type, local_size, local_sizePrior);
    
    NS::Error* err;
    NS::SharedPtr<MTL::CompileOptions> optsFP = NS::TransferPtr(MTL::CompileOptions::alloc()->init());
    NS::SharedPtr<MTL::CompileOptions> optsBP = NS::TransferPtr(MTL::CompileOptions::alloc()->init());
    optsFP->setPreprocessorMacros(sets.fp);
    optsBP->setPreprocessorMacros(sets.bp);
    libFP = NS::TransferPtr(mtlDevice->newLibrary(NS::String::string(contentFP.c_str(), NS::UTF8StringEncoding), optsFP.get(), &err));
    if (!libFP) {
        const char* msg = (err && err->localizedDescription())
            ? err->localizedDescription()->utf8String()
            : "unknown Metal compile error";
        mexPrintf("newLibrary failed: %s", msg);
        return -1;
    }
    libBP = NS::TransferPtr(mtlDevice->newLibrary(NS::String::string(contentBP.c_str(), NS::UTF8StringEncoding), optsBP.get(), &err));
    if (!libBP) {
        const char* msg = (err && err->localizedDescription())
            ? err->localizedDescription()->utf8String()
            : "unknown Metal compile error";
        mexPrintf("newLibrary failed: %s", msg);
        return -1;
    }

    if (inputScalars.computeSensImag) {

    }

    // Build prior programs
    if (MethodList.NLM || MethodList.MRP || MethodList.RDP || w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]
        || MethodList.TV || MethodList.APLS || MethodList.hyperbolic || MethodList.ProxTV || MethodList.ProxTGV || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM ||
        MethodList.CPType || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF || inputScalars.projector_type == 6 || type == 0) {
        if (DEBUG) {
            mexPrint("Building aux programs\n");
        }
    }

    err->release();
    return 0;
}

int ProjectorClass::createKernels(
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelFP,
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelBP,
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelNLM,
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelMed,
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelRDP,
    NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelGGMRF,
    const NS::SharedPtr<MTL::Library>& libFP,
    const NS::SharedPtr<MTL::Library>& libBP,
    const NS::SharedPtr<MTL::Library>& libAux,
    const NS::SharedPtr<MTL::Library>& libSens,
    const RecMethods& MethodList,
    const Weighting& w_vec,
    const scalarStruct& inputScalars,
    const int type
) {
    NS::Error* err = nullptr;
    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
            auto fnName = NS::String::string("projectorType123", NS::ASCIIStringEncoding);
            NS::SharedPtr<MTL::Function> fnFP = NS::TransferPtr(libFP->newFunction(fnName));
            NS::SharedPtr<MTL::ComputePipelineState> psoFP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnFP.get(), &err));
            NS::SharedPtr<MTL::CommandQueue> queueFP = NS::TransferPtr(mtlDevice->newCommandQueue());
            commandBufferFP = NS::TransferPtr(queueFP->commandBuffer());
            kernelFP = NS::TransferPtr(commandBufferFP->computeCommandEncoder());
            kernelFP->setComputePipelineState(psoFP.get());
        }
        if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
            auto fnName = NS::String::string("projectorType123", NS::ASCIIStringEncoding);
            NS::SharedPtr<MTL::Function> fnBP = NS::TransferPtr(libBP->newFunction(fnName));
            NS::SharedPtr<MTL::ComputePipelineState> psoBP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnBP.get(), &err));
            NS::SharedPtr<MTL::CommandQueue> queueBP = NS::TransferPtr(mtlDevice->newCommandQueue());
            commandBufferBP = NS::TransferPtr(queueBP->commandBuffer());
            kernelBP = NS::TransferPtr(commandBufferBP->computeCommandEncoder());
            kernelBP->setComputePipelineState(psoBP.get());
        }
    }
    return 0;
}

int ProjectorClass::addProjector(
    scalarStruct& inputScalars,
    Weighting& w_vec,
    const RecMethods& MethodList,
    const char* header_directory,
    const int type
) {
    // Set-up the local group size
    local_size[0] = 64ULL;
    local_size[1] = 1ULL;
    local_size[2] = 1ULL;
    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || (inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT)))
        local_size[0] = 128ULL;
    if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
        if (inputScalars.nColsD > 1 && !(inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT))) {
            local_size[0] = 16ULL;
            local_size[1] = 16ULL;
        }
    }
    if (DEBUG) {
        mexPrintBase("inputScalars.nColsD = %u\n", inputScalars.nColsD);
        mexPrintBase("inputScalars.nRowsD = %u\n", inputScalars.nRowsD);
        mexPrintBase("local_size[0] = %u\n", local_size[0]);
        mexPrintBase("local_size[1] = %u\n", local_size[1]);
        mexEval();
    }
    // Local group for priors
    local_sizePrior[0] = 16ULL;
    local_sizePrior[1] = 16ULL;
    local_sizePrior[2] = 1ULL;
    NS::Integer status = 0;
    proj6 = 0;

    NS::SharedPtr<MTL::Library> libFP, libBP, libAux, libSens;
    status = createProgram(libFP, libBP, libAux, libSens, header_directory, inputScalars, MethodList, w_vec, local_size, type);
    if (status != 0) return -1;
    status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, libFP, libBP, libAux, libSens, MethodList, w_vec, inputScalars, type);
    if (status != 0) return -1;

    if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
        erotus[0] = inputScalars.nRowsD % local_size[0];
        if (inputScalars.FPType == 5)
            erotus[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP) % local_size[1];
        else
            erotus[1] = inputScalars.nColsD % local_size[1];
        if (erotus[1] > 0)
            erotus[1] = (local_size[1] - erotus[1]);
        if (erotus[0] > 0)
            erotus[0] = (local_size[0] - erotus[0]);
    }

    if ((MethodList.ProxTGV || MethodList.ProxTV || MethodList.ProxRDP)) {
        erotusPriorEFOV[0] = inputScalars.NxPrior % local_sizePrior[0];
        erotusPriorEFOV[1] = inputScalars.NyPrior % local_sizePrior[1];
        erotusPriorEFOV[2] = inputScalars.NzPrior % local_sizePrior[2];
        if (erotusPriorEFOV[0] > 0)
            erotusPriorEFOV[0] = (local_sizePrior[0] - erotusPriorEFOV[0]);
        if (erotusPriorEFOV[1] > 0)
            erotusPriorEFOV[1] = (local_sizePrior[1] - erotusPriorEFOV[1]);
        if (erotusPriorEFOV[2] > 0)
            erotusPriorEFOV[2] = (local_sizePrior[1] - erotusPriorEFOV[2]);
        globalPriorEFOV = { (int)inputScalars.NxPrior + (int)erotusPriorEFOV[0], (int)inputScalars.NyPrior + (int)erotusPriorEFOV[1], (int)inputScalars.NzPrior + (int)erotusPriorEFOV[2] };
    }

    erotusBP.resize(2);
    erotusPDHG.resize(2);
    if (MethodList.CPType || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM) {
        for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
            erotusPDHG[0].emplace_back(inputScalars.Nx[ii] % local_sizePrior[0]);
            erotusPDHG[1].emplace_back(inputScalars.Ny[ii] % local_sizePrior[1]);
            if (erotusPDHG[0][ii] > 0)
                erotusPDHG[0][ii] = (local_sizePrior[0] - erotusPDHG[0][ii]);
            if (erotusPDHG[1][ii] > 0)
                erotusPDHG[1][ii] = (local_sizePrior[1] - erotusPDHG[1][ii]);
        }
    }
    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
        erotusBP[0].emplace_back(inputScalars.Nx[ii] % local_size[0]);
        erotusBP[1].emplace_back(inputScalars.Ny[ii] % local_size[1]);
        if (erotusBP[0][ii] > 0)
            erotusBP[0][ii] = (local_size[0] - erotusBP[0][ii]);
        if (erotusBP[1][ii] > 0)
            erotusBP[1][ii] = (local_size[1] - erotusBP[1][ii]);
    }
    local = { (int)local_size[0] , (int)local_size[1], 1 };
    localPrior = { (int)local_sizePrior[0] , (int)local_sizePrior[1], (int)local_sizePrior[2] };
    erotusPrior[0] = inputScalars.Nx[0] % local_sizePrior[0];
    erotusPrior[1] = inputScalars.Ny[0] % local_sizePrior[1];
    erotusPrior[2] = inputScalars.Nz[0] % local_sizePrior[2];
    if (erotusPrior[0] > 0)
        erotusPrior[0] = (local_sizePrior[0] - erotusPrior[0]);
    if (erotusPrior[1] > 0)
        erotusPrior[1] = (local_sizePrior[1] - erotusPrior[1]);
    if (erotusPrior[2] > 0)
        erotusPrior[2] = (local_sizePrior[1] - erotusPrior[2]);
    globalPrior = { (int)inputScalars.Nx[0] + (int)erotusPrior[0], (int)inputScalars.Ny[0] + (int)erotusPrior[1], (int)inputScalars.Nz[0] + (int)erotusPrior[2] };
    //d_NOrig = { static_cast<NS::Integer>(inputScalars.NxOrig), static_cast<NS::Integer>(inputScalars.NyOrig), static_cast<NS::Integer>(inputScalars.NzOrig) };
    //d_NPrior = { static_cast<NS::Integer>(inputScalars.NxPrior), static_cast<NS::Integer>(inputScalars.NyPrior), static_cast<NS::Integer>(inputScalars.NzPrior) };
    //dPitch = { w_vec.dPitchX, w_vec.dPitchY };
    b.resize(inputScalars.nMultiVolumes + 1);
    d.resize(inputScalars.nMultiVolumes + 1);
    d_N.resize(inputScalars.nMultiVolumes + 1);
    bmax.resize(inputScalars.nMultiVolumes + 1);
    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
        b[ii] = { inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii] };
        d[ii] = { inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii] };
        d_N[ii] = { static_cast<uint>(inputScalars.Nx[ii]), static_cast<uint>(inputScalars.Ny[ii]), static_cast<uint>(inputScalars.Nz[ii]) };
        bmax[ii] = { static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
            static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii],
            static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii] };
    }
    if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
        erotusSens[0] = inputScalars.det_per_ring % local_size[0];
        erotusSens[1] = inputScalars.det_per_ring % local_size[1];
        if (erotusSens[1] > 0)
            erotusSens[1] = (local_size[1] - erotusSens[1]);
        if (erotusSens[0] > 0)
            erotusSens[0] = (local_size[0] - erotusSens[0]);
    }
    region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
    return 0;
}

int ProjectorClass::createBuffers(
    scalarStruct& inputScalars,
    Weighting& w_vec,
    const float* x,
    const float* z_det,
    const uint32_t* xy_index,
    const uint16_t* z_index,
    const uint16_t* L,
    const int64_t* pituus,
    const float* atten,
    const float* norm,
    const float* extraCorr, 
    const std::vector<int64_t>& length,
    const RecMethods& MethodList,
    const int type
) {
    if (inputScalars.maskFP)
        d_maskFP.resize(inputScalars.subsetsUsed);
    if (inputScalars.raw)
        d_L.resize(inputScalars.subsetsUsed);
    if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1) {
        d_xyindex.resize(inputScalars.subsetsUsed);
        d_zindex.resize(inputScalars.subsetsUsed);
    }
    if (inputScalars.listmode > 0 && inputScalars.indexBased) {
        d_trIndex.resize(inputScalars.subsetsUsed);
        d_axIndex.resize(inputScalars.subsetsUsed);
    }
    if (inputScalars.listmode > 0 && inputScalars.TOF) {
        d_TOFIndex.resize(inputScalars.subsetsUsed);
    }
    if (inputScalars.normalization_correction)
        d_norm.resize(inputScalars.subsetsUsed);
    if (inputScalars.scatter)
        d_scat.resize(inputScalars.subsetsUsed);
    if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
        d_atten.resize(inputScalars.subsetsUsed);
    if (inputScalars.projector_type != 6) {
        d_x.resize(inputScalars.subsetsUsed);
        d_z.resize(inputScalars.subsetsUsed);
    }

    if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
        d_T.resize(inputScalars.subsetsUsed);

    size_t vecSize = 1;
    if ((inputScalars.PET 
        || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
        vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);

    size_t fpSize = sizeof(float); // Floating point size
    
    const MTL::ResourceOptions sharedOpts = (MTL::ResourceOptions)MTL::ResourceStorageModeShared;

    if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        NS::UInteger bytes = (NS::UInteger)(fpSize * (size_t)inputScalars.size_V);
        d_V = NS::TransferPtr(mtlDevice->newBuffer((const void*)inputScalars.V, bytes, sharedOpts));
    }

    // Detector coordinates
    if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
        //d_x[0] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x, NULL, &status);
    }

    // Attenuation data for image-based attenuation
    if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
        NS::UInteger bytes = (NS::UInteger)(fpSize * (size_t)inputScalars.im_dim[0]);
        d_attenB = NS::TransferPtr(mtlDevice->newBuffer((const void*)atten, bytes, sharedOpts));
    }
    
    // Forward projection mask
    if (inputScalars.maskFP) {
        if (inputScalars.maskFPZ > 1) {
            for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
                NS::UInteger bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD * (size_t)length[kk]);
                const uint8_t* src = &w_vec.maskFP[(size_t)pituus[kk] * (size_t)vecSize];
                d_maskFP[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)src, bytes, sharedOpts));
            }
        } else {
            NS::UInteger bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD);
            d_maskFP[0] = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.maskFP, bytes, sharedOpts));
        }
    }

    // Backprojection mask
    if (inputScalars.maskBP) {
        NS::UInteger bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.Nx[0] * (size_t)inputScalars.Ny[0] * (size_t)inputScalars.maskBPZ);
        d_maskBP = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.maskBP, bytes, sharedOpts));
    }

    // SPECT ray shifts (required for only projector types 1,2 and 3)
    if (inputScalars.SPECT) {
        NS::UInteger bytes = (NS::UInteger)(fpSize * 2ull * (size_t)inputScalars.n_rays * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD * (size_t)inputScalars.nProjections);
        d_rayShiftsDetector = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.rayShiftsDetector, bytes, sharedOpts));
        d_rayShiftsSource = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.rayShiftsSource, bytes, sharedOpts));
    }

    if (inputScalars.eFOV) {
        NS::UInteger bytes = (NS::UInteger)(sizeof(uint8_t) * inputScalars.Nz[0]);
        d_eFOVIndices = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.eFOVIndices, bytes, sharedOpts));
    }

    if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
        NS::UInteger bytes = (NS::UInteger)(sizeof(float) * inputScalars.nProjections);
        d_angle = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.angles, bytes, sharedOpts));
    }

    // TOF bin centers
    if (inputScalars.TOF) {
        NS::UInteger bytes = (NS::UInteger)(sizeof(float) * inputScalars.nBins);
        d_eFOVIndices = NS::TransferPtr(mtlDevice->newBuffer((const void*)inputScalars.TOFCenter, bytes, sharedOpts));
    }

    if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
        NS::UInteger bytesX = (NS::UInteger)(sizeof(float) * inputScalars.size_of_x);
        NS::UInteger bytesZ = (NS::UInteger)(sizeof(float) * inputScalars.size_z);
        d_xFull.emplace_back(NS::TransferPtr(mtlDevice->newBuffer((const void*)x, bytesX, sharedOpts)));
        d_zFull.emplace_back(NS::TransferPtr(mtlDevice->newBuffer((const void*)z_det, bytesZ, sharedOpts)));
    }

    // Per-subset buffers
    for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
        if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
            if (inputScalars.pitch) {
                NS::UInteger bytesZ = (NS::UInteger)(fpSize * (size_t)length[kk] * 6u);
                const float* srcZ = &z_det[(size_t)pituus[kk] * 6u];
                d_z[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
            } else {
                NS::UInteger bytesZ = (NS::UInteger)(fpSize * (size_t)length[kk] * 2u);
                const float* srcZ = &z_det[(size_t)pituus[kk] * 2u];
                d_z[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
            }
        } else {
            if (inputScalars.PET && inputScalars.listmode == 0) {
                if (inputScalars.nLayers > 1) {
                    NS::UInteger bytesZ = (NS::UInteger)(fpSize * (size_t)length[kk] * 3u);
                    const float* srcZ = &z_det[(size_t)pituus[kk] * 3u];
                    d_z[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
                } else {
                    NS::UInteger bytesZ = (NS::UInteger)(fpSize * (size_t)length[kk] * 2u);
                    const float* srcZ = &z_det[(size_t)pituus[kk] * 2u];
                    d_z[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
                }
            } else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
                NS::UInteger bytesZ = (NS::UInteger)(fpSize * (size_t)inputScalars.size_z);
                const float* srcZ = &z_det[(size_t)inputScalars.size_z];
                d_z[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
            }
        }
        
        if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
            //d_T[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
        }

        if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
            NS::UInteger bytesX = (NS::UInteger)(fpSize * (size_t)length[kk] * 6u);
            const float* srcX = &x[(size_t)pituus[kk] * 6u];
            d_x[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcX, bytesX, sharedOpts));
        } else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
            if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF) {
                NS::UInteger bytesX = (NS::UInteger)(fpSize * (size_t)length[kk] * 6u);
                const float* srcX = &w_vec.listCoord[pituus[kk] * 6u];
                d_x[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcX, bytesX, sharedOpts));
            }
        }

        // Normalization
        if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
            NS::UInteger bytesN = (NS::UInteger)(fpSize * (size_t)length[kk] * (size_t)vecSize);
            const float* srcN = &norm[(size_t)pituus[kk] * (size_t)vecSize];
            d_norm[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcN, bytesN, sharedOpts));
        }

        // Scatter
        if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) {
            NS::UInteger bytesS = (NS::UInteger)(fpSize * (size_t)length[kk] * (size_t)vecSize);
            const float* srcS = &extraCorr[(size_t)pituus[kk] * (size_t)vecSize];
            d_scat[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcS, bytesS, sharedOpts));
        }

        if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
            NS::UInteger bytes = sizeof(float) * length[kk] * vecSize;
            const float* src = &atten[pituus[kk] * vecSize];
            d_atten[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)src, bytes, sharedOpts));
        }

        // Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
        if (inputScalars.raw) {
            NS::UInteger bytesL = (NS::UInteger)(sizeof(uint16_t) * (size_t)length[kk] * 2u);
            const uint16_t* srcL = &L[(size_t)pituus[kk] * 2u];
            d_L[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcL, bytesL, sharedOpts));
        } else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
            NS::UInteger bytesXY = (NS::UInteger)(sizeof(uint32_t) * length[kk]);
            NS::UInteger bytesZ = (NS::UInteger)(sizeof(uint16_t) * length[kk]);
            const uint32_t* srcXY = &xy_index[pituus[kk]];
            const uint16_t* srcZ = &z_index[pituus[kk]];
            d_xyindex[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcXY, bytesXY, sharedOpts));
            d_zindex[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
        }

        if (inputScalars.listmode > 0 && (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF))) {
            if (inputScalars.indexBased) {
                NS::UInteger bytes_tr_ax = (NS::UInteger)(sizeof(uint16_t) * length[kk] * 2);
                const uint16_t* srcTR = &w_vec.trIndex[pituus[kk] * 2];
                const uint16_t* srcAX = &w_vec.axIndex[pituus[kk] * 2];
                d_trIndex[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcTR, bytes_tr_ax, sharedOpts));
                d_axIndex[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcAX, bytes_tr_ax, sharedOpts));
            }
            if (inputScalars.TOF) {
                NS::UInteger bytesTOFIdx = (NS::UInteger)(sizeof(uint8_t) * length[kk]);
                const uint8_t* srcTOFIdx = &w_vec.TOFIndices[pituus[kk]];
                d_TOFIndex[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcTOFIdx, bytesTOFIdx, sharedOpts));
            }
        }
    }
    return 0;
}

int ProjectorClass::initializeKernel(
    scalarStruct& inputScalars,
    Weighting& w_vec
) {
    StaticScalarKernelParams params = {};
    params.nRowsD = inputScalars.nRowsD;
    params.nColsD = inputScalars.nColsD;
    params.dPitchX = w_vec.dPitchX;
    params.dPitchY = w_vec.dPitchY;
    params.dL = inputScalars.dL;
    params.global_factor = inputScalars.global_factor;
    params.epps = inputScalars.epps;
    params.det_per_ring = inputScalars.det_per_ring;
    params.sigma_x = inputScalars.sigma_x;
    params.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
    params.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
    params.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
    params.bmin = inputScalars.bmin;
    params.bmax = inputScalars.bmax;
    params.Vmax = inputScalars.Vmax;
    params.rings = inputScalars.rings;

    // Set buffers to kernels
    kernelFP->setBytes((const void*)&params, (NS::UInteger)sizeof(params), 0);
    if (DEBUG) mexPrint("initializeKernel: FP buffer 0 (static params) set");
    kernelBP->setBytes((const void*)&params, (NS::UInteger)sizeof(params), 0);
    if (DEBUG) mexPrint("initializeKernel: BP buffer 0 (static params) set");

    return 0;
}

int ProjectorClass::setDynamicKernelData(
    scalarStruct& inputScalars,
    Weighting& w_vec
) {
    return 0;
}

int ProjectorClass::forwardProjection(
    const scalarStruct& inputScalars,
    Weighting& w_vec,
    const uint32_t osa_iter,
    const std::vector<int64_t>& length,
    const uint64_t m_size,
    const int32_t ii,
    const int uu
) {
    if (DEBUG) mexPrint("forwardProjection: init");
    if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
        global[0] = (inputScalars.nRowsD + erotus[0]);
        global[1] = (inputScalars.nColsD  + erotus[1]);
        global[2] = length[osa_iter];
    } else {
        erotus[0] = length[osa_iter] % local[0];

        if (erotus[0] > 0)
            erotus[0] = (local[0] - erotus[0]);
        global[0] = static_cast<size_t>(length[osa_iter] + erotus[0]);
        global[1] = 1;
        global[2] = 1;
    }

    DynamicScalarKernelParams params = {};
    params.d_N = d_N[ii];
    params.d = d[ii];
    params.b = b[ii];
    params.bmax = bmax[ii];
    params.nProjections = length[osa_iter];
    params.no_norm = no_norm;
    params.m_size = m_size;
    params.currentSubset = osa_iter;
    params.aa = ii;
    if (inputScalars.FPType == 2) params.orthWidth = inputScalars.tube_width;
    else if (inputScalars.FPType == 3) params.orthWidth = inputScalars.cylRadiusProj3;

    kernelFP->setBytes((const void*)&params, (NS::UInteger)sizeof(params), 1);
    if (DEBUG) mexPrint("forwardProjection: buffer 1 (dynamic params) set");

    if (inputScalars.SPECT) {
        kernelBP->setBuffer(d_rayShiftsDetector.get(), (NS::UInteger)0, /*index*/ 2);
        kernelBP->setBuffer(d_rayShiftsSource.get(),   (NS::UInteger)0, /*index*/ 3);
        if (DEBUG) mexPrint("forwardProjection: buffers 2 and 3 set (SPECT shifts)");
    }
    if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        kernelFP->setBuffer(d_V.get(), (NS::UInteger)0, /*index*/ 4);
    }

    if (inputScalars.attenuation_correction && !inputScalars.CT && inputScalars.CTAttenuation) {
        kernelFP->setBuffer(d_attenB.get(), (NS::UInteger)0, /*index*/ 5);
        if (DEBUG) mexPrint("forwardProjection: buffer 5 set");
    }

    if (inputScalars.maskFP) {
        int subset = 0;
        if (inputScalars.maskFPZ > 1) subset = osa_iter;
        kernelFP->setBuffer(d_maskFP[subset].get(), (NS::UInteger)0, /*index*/ 6);
        if (DEBUG) mexPrint("forwardProjection: buffer 6 set");
    }

    if (inputScalars.maskBP) {
        kernelFP->setBuffer(d_maskBP.get(), (NS::UInteger)0, /*index*/ 7);
        if (DEBUG) mexPrint("forwardProjection: buffer 7 set");
    }

    kernelFP->setBuffer(d_x[osa_iter].get(), (NS::UInteger)0, /*index*/ 8);
    kernelFP->setBuffer(d_z[osa_iter].get(), (NS::UInteger)0, /*index*/ 9);
    if (DEBUG) mexPrint("forwardProjection: buffers 8 and 9 set");

    if (inputScalars.normalization_correction) {
        kernelFP->setBuffer(d_norm[osa_iter].get(), (NS::UInteger)0, /*index*/ 10);
        if (DEBUG) mexPrint("forwardProjection: buffer 10 set");
    }
    if (inputScalars.scatter) {
        kernelFP->setBuffer(d_scat[osa_iter].get(), (NS::UInteger)0, /*index*/ 11);
        if (DEBUG) mexPrint("forwardProjection: buffer 11 set");
    }

    kernelFP->setBuffer(d_Summ[uu].get(), (NS::UInteger)0, /*index*/ 12);
    if (DEBUG) mexPrint("forwardProjection: buffer 12 set");

    if ( (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)
        && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0 )
    {
        kernelFP->setBuffer(d_xyindex[osa_iter].get(), (NS::UInteger)0, /*index*/ 13);
        kernelFP->setBuffer(d_zindex[osa_iter].get(),  (NS::UInteger)0, /*index*/ 14);
        if (DEBUG) mexPrint("forwardProjection: buffers 13 and 14 set");
    }

    // if (inputScalars.raw)
    //     kernelFP->setBuffer(d_L[osa_iter].get(), (NS::UInteger)0, /*index*/ 18);

    kernelFP->setBuffer(vec_opencl.d_im.get(), (NS::UInteger)0, /*index*/ 19);
    if (DEBUG) mexPrint("forwardProjection: buffer 19 set");
    kernelFP->setBuffer(d_output.get(), (NS::UInteger)0, /*index*/ 20);
    if (DEBUG) mexPrint("forwardProjection: buffer 20 set");
    if (DEBUG) mexPrint("forwardProjection: all buffers set");

    // Dispatch
    MTL::Size threadsPerThreadgroup = MTL::Size::Make(local[0], local[1], local[2]);
    MTL::Size threadgroupsPerGrid = MTL::Size::Make(global[0] / local[0], global[1] / local[1], global[2] / local[2]);
    
    kernelFP->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
    if (DEBUG) mexPrint("forwardProjection: dispatchThreadGroups ready");
    kernelFP->endEncoding();
    if (DEBUG) mexPrint("forwardProjection: endEncoding ready");

    commandBufferFP->commit();
    if (DEBUG) mexPrint("forwardProjection: commandBufferFP->commit() complete");
    commandBufferFP->waitUntilCompleted();
    if (DEBUG) mexPrint("forwardProjection: commandBufferFP->waitUntilCompleted() complete");
    return 0;
}

int ProjectorClass::backwardProjection(
    const scalarStruct& inputScalars,
    Weighting& w_vec,
    const uint32_t osa_iter,
    const std::vector<int64_t>& length,
    const uint64_t m_size,
    const bool compSens,
    const int32_t ii,
    const int uu,
    int ee
) {
    if (DEBUG) mexPrint("backwardProjection: init");
    if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
        global[0] = (inputScalars.nRowsD + erotus[0]);
        global[1] = (inputScalars.nColsD  + erotus[1]);
        global[2] = length[osa_iter];
    } else {
        erotus[0] = length[osa_iter] % local[0];

        if (erotus[0] > 0)
            erotus[0] = (local[0] - erotus[0]);
        global[0] = static_cast<size_t>(length[osa_iter] + erotus[0]);
        global[1] = 1;
        global[2] = 1;
    }

    DynamicScalarKernelParams params = {};
    params.d_N = d_N[ii];
    params.d = d[ii];
    params.b = b[ii];
    params.bmax = bmax[ii];
    params.nProjections = length[osa_iter];
    params.no_norm = no_norm;
    params.m_size = m_size;
    params.currentSubset = osa_iter;
    params.aa = ii;
    if (inputScalars.BPType == 2) params.orthWidth = inputScalars.tube_width;
    else if (inputScalars.BPType == 3) params.orthWidth = inputScalars.cylRadiusProj3;

    kernelBP->setBytes((const void*)&params, (NS::UInteger)sizeof(params), 1);
    if (DEBUG) mexPrint("backwardProjection: buffer 1 (dynamic params) set");

    if (inputScalars.SPECT) {
        kernelBP->setBuffer(d_rayShiftsDetector.get(), (NS::UInteger)0, /*index*/ 2);
        kernelBP->setBuffer(d_rayShiftsSource.get(),   (NS::UInteger)0, /*index*/ 3);
        if (DEBUG) mexPrint("backwardProjection: buffers 2 and 3 set (SPECT shifts)");
    }
    if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
        kernelBP->setBuffer(d_V.get(), (NS::UInteger)0, /*index*/ 4);
    }

    if (inputScalars.attenuation_correction && !inputScalars.CT && inputScalars.CTAttenuation) {
        kernelBP->setBuffer(d_attenB.get(), (NS::UInteger)0, /*index*/ 5);
        if (DEBUG) mexPrint("backwardProjection: buffer 5 set (attenB)");
    }

    if (inputScalars.maskFP) {
        int subset = 0;
        if (inputScalars.maskFPZ > 1) subset = osa_iter;
        kernelBP->setBuffer(d_maskFP[subset].get(), (NS::UInteger)0, /*index*/ 6);
        if (DEBUG) mexPrint("backwardProjection: buffer 6 set (maskFP)");
    }

    if (inputScalars.maskBP) {
        kernelBP->setBuffer(d_maskBP.get(), (NS::UInteger)0, /*index*/ 7);
        if (DEBUG) mexPrint("backwardProjection: buffer 7 set (maskBP)");
    }

    kernelBP->setBuffer(d_x[osa_iter].get(), (NS::UInteger)0, /*index*/ 8);
    kernelBP->setBuffer(d_z[osa_iter].get(), (NS::UInteger)0, /*index*/ 9);
    if (DEBUG) mexPrint("backwardProjection: buffers 8 and 9 set (x, z)");

    if (inputScalars.normalization_correction) {
        kernelBP->setBuffer(d_norm[osa_iter].get(), (NS::UInteger)0, /*index*/ 10);
        if (DEBUG) mexPrint("backwardProjection: buffer 10 set (norm)");
    }
    if (inputScalars.scatter) {
        kernelBP->setBuffer(d_scat[osa_iter].get(), (NS::UInteger)0, /*index*/ 11);
        if (DEBUG) mexPrint("backwardProjection: buffer 11 set (scat)");
    }

    kernelBP->setBuffer(d_Summ[uu].get(), (NS::UInteger)0, /*index*/ 12);
    if (DEBUG) mexPrint("backwardProjection: buffer 12 set (Summ)");

    if ( (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)
        && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0 )
    {
        kernelBP->setBuffer(d_xyindex[osa_iter].get(), (NS::UInteger)0, /*index*/ 13);
        kernelBP->setBuffer(d_zindex[osa_iter].get(),  (NS::UInteger)0, /*index*/ 14);
        if (DEBUG) mexPrint("backwardProjection: buffers 13 and 14 set (indices)");
    }

    // if (inputScalars.raw)
    //     kernelBP->setBuffer(d_L[osa_iter].get(), (NS::UInteger)0, /*index*/ 18);

    kernelBP->setBuffer(d_output.get(),   (NS::UInteger)0, /*index*/ 19);
    if (DEBUG) mexPrint("backwardProjection: buffer 19 set (output)");
    kernelBP->setBuffer(vec_opencl.d_rhs_os[uu].get(), (NS::UInteger)0, /*index*/ 20);
    if (DEBUG) mexPrint("backwardProjection: buffer 20 set (rhs_os)");

    if (DEBUG) mexPrint("backwardProjection: all buffers set");

    MTL::Size threadsPerThreadgroup = MTL::Size::Make(local[0], local[1], local[2]);
    MTL::Size threadgroupsPerGrid = MTL::Size::Make(global[0] / local[0], global[1] / local[1], global[2] / local[2]);

    kernelBP->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
    if (DEBUG) mexPrint("backwardProjection: dispatchThreadgroups ready");

    kernelBP->endEncoding();
    if (DEBUG) mexPrint("backwardProjection: endEncoding ready");

    commandBufferBP->commit();
    if (DEBUG) mexPrint("backwardProjection: commandBufferBP->commit() complete");
    commandBufferBP->waitUntilCompleted();
    if (DEBUG) mexPrint("backwardProjection: commandBufferBP->waitUntilCompleted() complete");

    return 0;
};