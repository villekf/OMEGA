//#if defined(METAL) // Metal
#include "ProjectorClassMetal.hpp"
//#elif defined(CUDA) // CUDA

//#else // OpenCL

//#endif

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
    if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
        auto fnName = NS::String::string("projectorType123", NS::ASCIIStringEncoding);
        NS::SharedPtr<MTL::Function> fnFP = NS::TransferPtr(libFP->newFunction(fnName));
        NS::SharedPtr<MTL::ComputePipelineState> psoFP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnFP.get(), &err));
        NS::SharedPtr<MTL::CommandQueue> queueFP = NS::TransferPtr(mtlDevice->newCommandQueue());
        commandBufferFP = NS::TransferPtr(queueFP->commandBuffer());
        kernelFP = NS::TransferPtr(commandBufferFP->computeCommandEncoder());
        kernelFP->setComputePipelineState(psoFP.get());
    } else if (inputScalars.FPType == 4) {
        auto fnName = NS::String::string("projectorType4Forward", NS::ASCIIStringEncoding);
        NS::SharedPtr<MTL::Function> fnFP = NS::TransferPtr(libFP->newFunction(fnName));
        NS::SharedPtr<MTL::ComputePipelineState> psoFP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnFP.get(), &err));
        NS::SharedPtr<MTL::CommandQueue> queueFP = NS::TransferPtr(mtlDevice->newCommandQueue());
        commandBufferFP = NS::TransferPtr(queueFP->commandBuffer());
        kernelFP = NS::TransferPtr(commandBufferFP->computeCommandEncoder());
        kernelFP->setComputePipelineState(psoFP.get());
    } else if (inputScalars.FPType == 5) {
        auto fnName = NS::String::string("projectorType5Forward", NS::ASCIIStringEncoding);
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
    } if (inputScalars.BPType == 4) {
        if (inputScalars.FPType == 4 && inputScalars.CT) {
            auto fnName = NS::String::string("projectorType4Backward", NS::ASCIIStringEncoding);
            NS::SharedPtr<MTL::Function> fnBP = NS::TransferPtr(libBP->newFunction(fnName));
            NS::SharedPtr<MTL::ComputePipelineState> psoBP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnBP.get(), &err));
            NS::SharedPtr<MTL::CommandQueue> queueBP = NS::TransferPtr(mtlDevice->newCommandQueue());
            commandBufferBP = NS::TransferPtr(queueBP->commandBuffer());
            kernelBP = NS::TransferPtr(commandBufferBP->computeCommandEncoder());
            kernelBP->setComputePipelineState(psoBP.get());
        } else if (!inputScalars.CT) {
            auto fnName = NS::String::string("projectorType4Forward", NS::ASCIIStringEncoding);
            NS::SharedPtr<MTL::Function> fnBP = NS::TransferPtr(libBP->newFunction(fnName));
            NS::SharedPtr<MTL::ComputePipelineState> psoBP = NS::TransferPtr(mtlDevice->newComputePipelineState(fnBP.get(), &err));
            NS::SharedPtr<MTL::CommandQueue> queueBP = NS::TransferPtr(mtlDevice->newCommandQueue());
            commandBufferBP = NS::TransferPtr(queueBP->commandBuffer());
            kernelBP = NS::TransferPtr(commandBufferBP->computeCommandEncoder());
            kernelBP->setComputePipelineState(psoBP.get());
        } else {
            auto fnName = NS::String::string("projectorType4Backward", NS::ASCIIStringEncoding);
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
        d_maskFPB.resize(inputScalars.subsetsUsed);
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
    if (inputScalars.normalization_correction)
        d_norm.resize(inputScalars.subsetsUsed);
    if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
        d_atten.resize(inputScalars.subsetsUsed);

    d_scat.resize(inputScalars.Nt);
    d_x.resize(inputScalars.Nt);
    d_z.resize(inputScalars.Nt);
    d_trIndex.resize(inputScalars.Nt);
    d_axIndex.resize(inputScalars.Nt);
    d_TOFIndex.resize(inputScalars.Nt);
    for (int tt = 0; tt < inputScalars.Nt; tt++) {
        d_scat[tt].resize(inputScalars.subsetsUsed);
        d_x[tt].resize(inputScalars.subsetsUsed);
        d_z[tt].resize(inputScalars.subsetsUsed);
        d_trIndex[tt].resize(inputScalars.subsetsUsed);
        d_axIndex[tt].resize(inputScalars.subsetsUsed);
        d_TOFIndex[tt].resize(inputScalars.subsetsUsed);
    }

    if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
        d_T.resize(inputScalars.subsetsUsed);

    size_t vecSize = 1;
    if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
        vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
    
    const MTL::ResourceOptions sharedOpts = (MTL::ResourceOptions)MTL::ResourceStorageModeShared;
    NS::UInteger bytes;

    NS::UInteger imX = inputScalars.Nx[0];
    NS::UInteger imY = inputScalars.Ny[0];
    NS::UInteger imZ = inputScalars.Nz[0];

    if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.im_dim[0]);
        //if (inputScalars.useImages)
        //    d_urefIm = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
        //else
            d_uref = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.NLM_ref, bytes, sharedOpts));
    }
    if (MethodList.NLM || MethodList.RDP || MethodList.TV || MethodList.GGMRF || MethodList.APLS || MethodList.hyperbolic || inputScalars.projector_type == 6) {
        if (inputScalars.useImages && !inputScalars.largeDim) {
            //d_inputI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
        }
    }
    if (MethodList.RDP && w_vec.RDPLargeNeighbor && w_vec.RDP_anatomical) {
        if (inputScalars.useImages) {
            //d_RDPrefI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
        }
    }
    if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
        bytes = (NS::UInteger)(sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1);
        d_weights = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.weights, bytes, sharedOpts));
    }
    if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
        if (inputScalars.useBuffers) {
            bytes = (NS::UInteger)(sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ);
            d_maskPriorB = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.maskPrior, bytes, sharedOpts));
        } else {
            imX = inputScalars.Nx[0];
            imY = inputScalars.Ny[0];
            imZ = inputScalars.maskBPZ;
            //if (imZ > 1)
            //    d_maskPrior3 = cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, imZ, 0, 0, NULL, &status);
            //else
            //    d_maskPrior = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
        }
    }

    if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        bytes = (NS::UInteger)(sizeof(float) * (size_t)inputScalars.size_V);
        d_V = NS::TransferPtr(mtlDevice->newBuffer((const void*)inputScalars.V, bytes, sharedOpts));
    }

    // Detector coordinates
    if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.size_of_x);
        d_x[0][0] = NS::TransferPtr(mtlDevice->newBuffer((const void*)x, bytes, sharedOpts));
    }
    
    // Forward projection mask
    if (inputScalars.maskFP) {
        if (inputScalars.useBuffers) {
            if (inputScalars.maskFPZ > 1) {
                for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
                    bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD * (size_t)length[kk]);
                    const uint8_t* src = &w_vec.maskFP[(size_t)pituus[kk] * (size_t)vecSize];
                    d_maskFPB[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)src, bytes, sharedOpts));
                }
            } else {
                bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD);
                d_maskFPB[0] = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.maskFP, bytes, sharedOpts));
            }
        } else {
            imX = inputScalars.nRowsD;
            imY = inputScalars.nColsD;
            imZ = inputScalars.maskFPZ;
            //if (imZ > 1) {
            //    for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
            //        d_maskFP3.emplace_back(cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, length[kk], 0, 0, NULL, &status));
            //}
            //else
            //    d_maskFP = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
        }
    }

    // Backprojection mask
    if (inputScalars.maskBP) {
        if (inputScalars.useBuffers) {
            bytes = (NS::UInteger)((size_t)sizeof(uint8_t) * (size_t)inputScalars.Nx[0] * (size_t)inputScalars.Ny[0] * (size_t)inputScalars.maskBPZ);
            d_maskBPB = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.maskBP, bytes, sharedOpts));
        } else {
            imX = inputScalars.Nx[0];
            imY = inputScalars.Ny[0];
            imZ = inputScalars.maskBPZ;
            //if (imZ > 1)
            //    d_maskBP3 = cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, imZ, 0, 0, NULL, &status);
            //else
            //    d_maskBP = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
        }
    }

    // SPECT ray shifts (required for only projector types 1,2 and 3)
    if (inputScalars.SPECT) {
        bytes = (NS::UInteger)(sizeof(float) * 2ull * (size_t)inputScalars.n_rays * (size_t)inputScalars.nRowsD * (size_t)inputScalars.nColsD * (size_t)inputScalars.nProjections);
        d_rayShiftsDetector = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.rayShiftsDetector, bytes, sharedOpts));
        d_rayShiftsSource = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.rayShiftsSource, bytes, sharedOpts));
    }

    if (inputScalars.eFOV) {
        bytes = (NS::UInteger)(sizeof(uint8_t) * inputScalars.Nz[0]);
        d_eFOVIndices = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.eFOVIndices, bytes, sharedOpts));
    }

    if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.nProjections);
        d_angle = NS::TransferPtr(mtlDevice->newBuffer((const void*)w_vec.angles, bytes, sharedOpts));
    }

    // TOF bin centers
    if (inputScalars.TOF) {
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.nBins);
        d_eFOVIndices = NS::TransferPtr(mtlDevice->newBuffer((const void*)inputScalars.TOFCenter, bytes, sharedOpts));
    }

    if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.size_of_x);
        d_xFull.emplace_back(NS::TransferPtr(mtlDevice->newBuffer((const void*)x, bytes, sharedOpts)));
        bytes = (NS::UInteger)(sizeof(float) * inputScalars.size_z);
        d_zFull.emplace_back(NS::TransferPtr(mtlDevice->newBuffer((const void*)z_det, bytes, sharedOpts)));
    }

    // Per-subset / per-timestep buffers
    for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
        // Attenuation data for image-based attenuation
        if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
            if (inputScalars.useBuffers){
                bytes = (NS::UInteger)(sizeof(float) * (size_t)inputScalars.im_dim[0]);
                d_attenB[timestep] = NS::TransferPtr(mtlDevice->newBuffer((const void*)atten, bytes, sharedOpts));
            } else {
                //d_attenIm[timestep] = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
            }
        }

        for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
            if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
                NS::UInteger bytesX = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * 6u);
                const float* srcX = &x[(size_t)pituus[kk] * 6 + (size_t)pituus[inputScalars.subsets] * 6 * timestep];
                d_x[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcX, bytesX, sharedOpts));
            } else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
                if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF) {
                    NS::UInteger bytesX = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * 6u);
                    const float* srcX = &w_vec.listCoord[pituus[kk] * 6u + inputScalars.kokoNonTOF * 6 * timestep];
                    d_x[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcX, bytesX, sharedOpts));
                }
            }

            NS::UInteger z_coef = 1;
            if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
                if (inputScalars.pitch) {
                    z_coef = 6;
                } else {
                    z_coef = 2;
                }
                bytes = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * z_coef);
                const float* srcZ = &z_det[(size_t)pituus[kk] * z_coef + pituus[inputScalars.subsets] * z_coef * timestep];
                d_z[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytes, sharedOpts));
            } else {
                if (inputScalars.PET && inputScalars.listmode == 0) {
                    if (inputScalars.nLayers > 1) {
                        z_coef = 3;
                    } else {
                        z_coef = 2;
                    }
                    bytes = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * z_coef);
                    const float* srcZ = &z_det[(size_t)pituus[kk] * z_coef + pituus[inputScalars.subsets] * z_coef * timestep];
                    d_z[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytes, sharedOpts));
                } else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
                    NS::UInteger bytesZ = (NS::UInteger)(sizeof(float) * (size_t)inputScalars.size_z);
                    const float* srcZ = &z_det[(size_t)inputScalars.size_z];
                    d_z[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcZ, bytesZ, sharedOpts));
                }
            }

            // Scatter
            if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) {
                NS::UInteger bytesS = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * (size_t)vecSize);
                const float* srcS = &extraCorr[pituus[kk] * vecSize + inputScalars.kokoNonTOF * timestep];
                d_scat[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcS, bytesS, sharedOpts));
            }

            if (inputScalars.listmode > 0 && (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF))) {
                if (inputScalars.indexBased) {
                    NS::UInteger bytes_tr_ax = (NS::UInteger)(sizeof(uint16_t) * length[kk] * 2);
                    const uint16_t* srcTR = &w_vec.trIndex[pituus[kk] * 2 + inputScalars.kokoNonTOF * 2 * timestep];
                    const uint16_t* srcAX = &w_vec.axIndex[pituus[kk] * 2 + inputScalars.kokoNonTOF * 2 * timestep];
                    d_trIndex[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcTR, bytes_tr_ax, sharedOpts));
                    d_axIndex[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcAX, bytes_tr_ax, sharedOpts));
                }
                if (inputScalars.TOF) {
                    NS::UInteger bytesTOFIdx = (NS::UInteger)(sizeof(uint8_t) * length[kk]);
                    const uint8_t* srcTOFIdx = &w_vec.TOFIndices[pituus[kk] + inputScalars.kokoNonTOF * timestep];
                    d_TOFIndex[timestep][kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcTOFIdx, bytesTOFIdx, sharedOpts));
                }
            }
        }
    }

    for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
        if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
            bytes = (NS::UInteger)(sizeof(float) * length[kk]);
            const float* src = &inputScalars.T[pituus[kk]];
            d_T[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)src, bytes, sharedOpts));
        }

        // Normalization
        if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
            NS::UInteger bytesN = (NS::UInteger)(sizeof(float) * (size_t)length[kk] * (size_t)vecSize);
            const float* srcN = &norm[(size_t)pituus[kk] * (size_t)vecSize];
            d_norm[kk] = NS::TransferPtr(mtlDevice->newBuffer((const void*)srcN, bytesN, sharedOpts));
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
    }
    return 0;
}

int ProjectorClass::initializeKernel(
    scalarStruct& inputScalars,
    Weighting& w_vec
) {
    // Set scalar struct
    kParams.nRowsD = inputScalars.nRowsD;
    kParams.nColsD = inputScalars.nColsD;
    kParams.dPitch = {w_vec.dPitchX, w_vec.dPitchY};
    kParams.dL = inputScalars.dL;
    kParams.global_factor = inputScalars.global_factor;
    kParams.epps = inputScalars.epps;
    kParams.det_per_ring = inputScalars.det_per_ring;
    kParams.sigma_x = inputScalars.sigma_x;
    kParams.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
    kParams.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
    kParams.coneOfResponseStdCoeffC = inputScalars.coneOfResponseStdCoeffC;
    kParams.bmin = inputScalars.bmin;
    kParams.bmax = inputScalars.bmax;
    kParams.Vmax = inputScalars.Vmax;
    kParams.rings = inputScalars.rings;
    kParams.helicalRadius = inputScalars.helicalRadius;

    // Set buffers to kernels
    if (inputScalars.FPType >= 1 && inputScalars.FPType <= 3) {
        if (inputScalars.SPECT) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_rayShiftsDetector, 0, 1);
            if (DEBUG) mexPrint("initializeKernel: FP buffer 1 (rayShiftsDetector) set");
            SET_KERNEL_ARG_BUFFER(kernelFP, d_rayShiftsSource, 0, 2);
            if (DEBUG) mexPrint("initializeKernel: FP buffer 2 (rayShiftsSource) set");
        }
        if (inputScalars.TOF) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFCenter, 0, 3);
            if (DEBUG) mexPrint("initializeKernel: FP buffer 3 (d_TOFCenter) set");
        }
        if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_V, 0, 4);
            if (DEBUG) mexPrint("initializeKernel: FP buffer 4 (d_V) set");
        }
    }

    if (inputScalars.BPType >= 1 && inputScalars.BPType <= 3) {
        if (inputScalars.SPECT) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_rayShiftsDetector, 0, 1);
            if (DEBUG) mexPrint("initializeKernel: BP buffer 1 (rayShiftsDetector) set");
            SET_KERNEL_ARG_BUFFER(kernelBP, d_rayShiftsSource, 0, 2);
            if (DEBUG) mexPrint("initializeKernel: BP buffer 2 (rayShiftsSource) set");
        }
        if (inputScalars.TOF) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_TOFCenter, 0, 3);
            if (DEBUG) mexPrint("initializeKernel: BP buffer 3 (d_TOFCenter) set");
        }
        if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_V, 0, 4);
            if (DEBUG) mexPrint("initializeKernel: BP buffer 4 (d_V) set");
        }
        /*if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
            if (inputScalars.SPECT) {
                kernelSensList->setBytes((const void*)&d_rayShiftsDetector, (NS::UInteger)sizeof(d_rayShiftsDetector), 1);
                if (DEBUG) mexPrint("initializeKernel: Sens buffer 1 (d_rayShiftsDetector) set");
                kernelSensList->setBytes((const void*)&d_rayShiftsSource, (NS::UInteger)sizeof(d_rayShiftsSource), 2);
                if (DEBUG) mexPrint("initializeKernel: Sens buffer 2 (d_rayShiftsSource) set");
            }
            if (inputScalars.TOF) {
                kernelSensList->setBytes((const void*)&d_TOFCenter, (NS::UInteger)sizeof(d_TOFCenter), 3);
                if (DEBUG) mexPrint("initializeKernel: Sens buffer 3 (d_TOFCenter) set");
            }
            if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
                kernelSensList->setBytes((const void*)&d_V, (NS::UInteger)sizeof(d_V), 4);
                if (DEBUG) mexPrint("initializeKernel: Sens buffer 4 (d_V) set");
            }
        }*/
    }

    if ((inputScalars.BPType == 4 || inputScalars.FPType == 4) && !inputScalars.CT && inputScalars.TOF) {
        if (inputScalars.FPType == 4) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFCenter, 0, 1);
            if (DEBUG) mexPrint("initializeKernel: FP buffer 1 (d_TOFCenter) set");
        }
        if (inputScalars.BPType == 4) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_TOFCenter, 0, 1);
            if (DEBUG) mexPrint("initializeKernel: BP buffer 1 (d_TOFCenter) set");
            if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
                //kernelSensList.setArg(kernelIndSens++, d_TOFCenter);
            }
        }
    }

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
    const uint32_t timestep,
    const std::vector<int64_t>& length,
    const uint64_t m_size,
    const int32_t ii,
    const int uu
) {
    if (DEBUG) mexPrint("forwardProjection: init");
    if (inputScalars.FPType == 5) {
        global[0] = (inputScalars.nRowsD + erotus[0]);
        global[1] = (inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1];
        global[2] = length[osa_iter];
    } else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
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

    std::chrono::steady_clock::time_point tStart;
    std::chrono::steady_clock::time_point tEnd;
    if (DEBUG || inputScalars.verbose >= 3) {
        tStart = std::chrono::steady_clock::now();
    }

    kParams.d_N = d_N[ii];
    kParams.d = d[ii];
    kParams.b = b[ii];
    kParams.d_bmax = bmax[ii];
    kParams.d_Scale4 = inputScalars.d_Scale4[ii];
    kParams.d_Scale5 = inputScalars.d_Scale[ii];
    kParams.dSize5 = inputScalars.dSize[ii];
    kParams.rings = inputScalars.rings;
    kParams.det_per_ring = inputScalars.det_per_ring;
    kParams.nProjections = length[osa_iter];
    kParams.no_norm = no_norm;
    kParams.m_size = m_size;
    kParams.currentSubset = osa_iter;
    kParams.aa = ii;
    if (inputScalars.FPType == 2) kParams.orthWidth = inputScalars.tube_width; // Set here to allow for projector types 23 and 32
    if (inputScalars.FPType == 3) kParams.orthWidth = inputScalars.cylRadiusProj3; // Set here to allow for projector types 23 and 32

    SET_KERNEL_ARG_BYTES(kernelFP, kParams, sizeof(kParams), 0);
    if (DEBUG) mexPrint("forwardProjection: FP buffer 0 (scalar params) set");

    if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
        if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
            //SET_KERNEL_ARG_BUFFER(kernelFP, d_atten[timestep][osa_iter], 0, 5);
        } else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
            if (inputScalars.useBuffers) {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_attenB[timestep], 0, 5);
            } else {
                //SET_KERNEL_ARG_TEXTURE(kernelFP, d_attenIm[timestep], 0, 5);
            }
        }
        if (DEBUG) mexPrint("forwardProjection: buffer 5 set");
  
        if (inputScalars.maskFP) {
            if (inputScalars.useBuffers) {
                int subset = 0;
                if (inputScalars.maskFPZ > 1) subset = osa_iter;
                SET_KERNEL_ARG_BUFFER(kernelFP, d_maskFPB[subset], 0, 6);
            } else {
                if (inputScalars.maskFPZ > 1) {
                    //kernelFP->setBuffer(d_maskFP3[osa_iter].get(), (NS::UInteger)0, 6);
                } else {
                    //kernelFP->setBuffer(d_maskFP[osa_iter].get(), (NS::UInteger)0, 6);
                }
            }
            if (DEBUG) mexPrint("forwardProjection: buffer 6 set");
        }

        if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[0][0], 0, 8);
        else
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[timestep][osa_iter], 0, 8);
        if (DEBUG) mexPrint("forwardProjection: buffer 8 set");

        if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
            SET_KERNEL_ARG_BUFFER(kernelFP, d_z[timestep][osa_iter], 0, 9);
        else
            SET_KERNEL_ARG_BUFFER(kernelFP, d_z[timestep][inputScalars.osa_iter0], 0, 9);
        if (DEBUG) mexPrint("forwardProjection: buffer 9 set");

        if (inputScalars.normalization_correction) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_norm[osa_iter], 0, 10);
            if (DEBUG) mexPrint("forwardProjection: buffer 10 set");
        }

        if (inputScalars.scatter) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_scat[timestep][osa_iter], 0, 11);
            if (DEBUG) mexPrint("forwardProjection: buffer 11 set");
        }

        SET_KERNEL_ARG_BUFFER(kernelFP, d_Summ[uu], 0, 12);
        if (DEBUG) mexPrint("forwardProjection: buffer 12 set");

        if ( (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0 ) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_xyindex[osa_iter], 0, 13);
            SET_KERNEL_ARG_BUFFER(kernelFP, d_zindex[osa_iter], 0, 14);
            if (DEBUG) mexPrint("forwardProjection: buffers 13 and 14 set");
        }

        if (inputScalars.listmode > 0 && inputScalars.indexBased) {
            if (!inputScalars.loadTOF) {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_trIndex[0][0], 0, 15);
                SET_KERNEL_ARG_BUFFER(kernelFP, d_axIndex[0][0], 0, 16);
            } else {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_trIndex[timestep][osa_iter], 0, 15);
                SET_KERNEL_ARG_BUFFER(kernelFP, d_axIndex[timestep][osa_iter], 0, 16);
            }
        }
        if (inputScalars.listmode > 0 && inputScalars.TOF) {
            if (!inputScalars.loadTOF) {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFIndex[0][0], 0, 17);
            } else {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFIndex[timestep][osa_iter], 0, 17);
            }
        }

        if (inputScalars.raw)
            SET_KERNEL_ARG_BUFFER(kernelFP, d_L[osa_iter], 0, 18);

        if (inputScalars.useBuffers) {
            SET_KERNEL_ARG_BUFFER(kernelFP, vec_opencl.d_im, 0, 19);
        } else {
            SET_KERNEL_ARG_TEXTURE(kernelFP, vec_opencl.d_image_os.get(), 19);
        }
        if (DEBUG) mexPrint("forwardProjection: buffer 19 set");

        SET_KERNEL_ARG_BUFFER(kernelFP, d_output, 0, 20);
        if (DEBUG) mexPrint("forwardProjection: buffer 20 set");
    }
    if (inputScalars.FPType == 4) {
        if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_atten[osa_iter], 0, 2);
        } else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
            if (inputScalars.useBuffers)
                SET_KERNEL_ARG_BUFFER(kernelFP, d_attenB[timestep], 0, 2);
            //else
            //    SET_KERNEL_ARG_TEXTURE(kernelFP, d_attenIm[timestep], 0, 2);
        }

        SET_KERNEL_ARG_TEXTURE(kernelFP, vec_opencl.d_image_os.get(), 3);
        SET_KERNEL_ARG_BUFFER(kernelFP, d_output, 0, 4);

        if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[0][0], 0, 5);
        else
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[timestep][osa_iter], 0, 5);

        if ((inputScalars.CT || inputScalars.PET || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
            SET_KERNEL_ARG_BUFFER(kernelFP, d_z[timestep][osa_iter], 0, 6);
        else
            SET_KERNEL_ARG_BUFFER(kernelFP, d_z[timestep][inputScalars.osa_iter0], 0, 6);

        if (inputScalars.maskFP) {
            if (inputScalars.useBuffers) {
                int subset = 0;
                if (inputScalars.maskFPZ > 1) subset = osa_iter;
                SET_KERNEL_ARG_BUFFER(kernelFP, d_maskFPB[subset], 0, 7);
            } else {
                if (inputScalars.maskFPZ > 1) {
                    //kernelFP->setBuffer(d_maskFP3[osa_iter].get(), (NS::UInteger)0, 7);
                } else {
                    //kernelFP->setBuffer(d_maskFP[osa_iter].get(), (NS::UInteger)0, 7);
                }
            }
        }

        if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_xyindex[osa_iter], 0, 9);
            SET_KERNEL_ARG_BUFFER(kernelFP, d_zindex[osa_iter], 0, 10);
        }
        if (inputScalars.listmode > 0 && inputScalars.indexBased) {
            if (!inputScalars.loadTOF) { // The data is loaded step by step to the first buffer
                SET_KERNEL_ARG_BUFFER(kernelFP, d_trIndex[0][0], 0, 9);
                SET_KERNEL_ARG_BUFFER(kernelFP, d_axIndex[0][0], 0, 10);
            } else { // All buffers are populated with data
                SET_KERNEL_ARG_BUFFER(kernelFP, d_trIndex[timestep][osa_iter], 0, 9);
                SET_KERNEL_ARG_BUFFER(kernelFP, d_axIndex[timestep][osa_iter], 0, 10);
            }
        }
        if (inputScalars.listmode > 0 && inputScalars.TOF) {
            if (!inputScalars.loadTOF) {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFIndex[0][0], 0, 11);
            } else {
                SET_KERNEL_ARG_BUFFER(kernelFP, d_TOFIndex[timestep][osa_iter], 0, 11);
            }
        }
        if (inputScalars.raw) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_L[osa_iter], 0, 12);
        }
        if (inputScalars.normalization_correction)
            SET_KERNEL_ARG_BUFFER(kernelFP, d_norm[osa_iter], 0, 13);
        if (inputScalars.scatter)
            SET_KERNEL_ARG_BUFFER(kernelFP, d_scat[timestep][osa_iter], 0, 14);
    }
    if (inputScalars.FPType == 5) {
        if (!inputScalars.loadTOF && inputScalars.listmode > 0)
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[0][0], 0, 1);
        else
            SET_KERNEL_ARG_BUFFER(kernelFP, d_x[timestep][osa_iter], 0, 1);
        SET_KERNEL_ARG_BUFFER(kernelFP, d_z[timestep][osa_iter], 0, 2);
        SET_KERNEL_ARG_TEXTURE(kernelFP, vec_opencl.d_image_os.get(), 3);
        SET_KERNEL_ARG_TEXTURE(kernelFP, vec_opencl.d_image_os_int.get(), 4);
        SET_KERNEL_ARG_BUFFER(kernelFP, d_output, 0, 5);

        if (inputScalars.meanFP) {
            SET_KERNEL_ARG_BUFFER(kernelFP, d_meanFP, 0, 6);
        }

        if (inputScalars.maskFP) {
            if (inputScalars.useBuffers) {
                int subset = 0;
                if (inputScalars.maskFPZ > 1) subset = osa_iter;
                SET_KERNEL_ARG_BUFFER(kernelFP, d_maskFPB[subset], 0, 7);
            } else {
                if (inputScalars.maskFPZ > 1) {
                    //kernelFP->setBuffer(d_maskFP3[osa_iter].get(), (NS::UInteger)0, 7);
                } else {
                    //kernelFP->setBuffer(d_maskFP[osa_iter].get(), (NS::UInteger)0, 7);
                }
            }
        }

        if (inputScalars.normalization_correction)
            SET_KERNEL_ARG_BUFFER(kernelFP, d_norm[osa_iter], 0, 8);
    }
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
    const uint32_t timestep,
    const std::vector<int64_t>& length,
    const uint64_t m_size,
    const bool compSens,
    const int32_t ii,
    const int uu,
    int ee
) {
    std::chrono::steady_clock::time_point tStart;
    std::chrono::steady_clock::time_point tEnd;
    if (DEBUG || inputScalars.verbose >= 3) {
        tStart = std::chrono::steady_clock::now();
    }

    if (DEBUG) mexPrint("backwardProjection: init");

    kParams.d_N = d_N[ii];
    kParams.d = d[ii];
    kParams.b = b[ii];
    kParams.d_bmax = bmax[ii];
    kParams.d_Scale5 = inputScalars.d_Scale[ii];
    kParams.dSize5 = inputScalars.dSizeBP;
    kParams.nProjections = length[osa_iter];
    kParams.no_norm = no_norm;
    kParams.m_size = m_size;
    kParams.currentSubset = osa_iter;
    kParams.aa = ii;
    if (inputScalars.BPType == 2) kParams.orthWidth = inputScalars.tube_width;
    if (inputScalars.BPType == 3) kParams.orthWidth = inputScalars.cylRadiusProj3;
    if (inputScalars.BPType == 4) kParams.kerroin4 = w_vec.kerroin4[ii];

    SET_KERNEL_ARG_BYTES(kernelBP, kParams, sizeof(kParams), 0);
    if (DEBUG) mexPrint("backwardProjection: BP buffer 0 (static params) set");

    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
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

        if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
            // SET_KERNEL_ARG_BUFFER(kernelBP, d_atten[timestep][osa_iter], 0, 5);
            if (DEBUG) mexPrint("backwardProjection: buffer 5 set (attenB)");
        } else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
            if (inputScalars.useBuffers) {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_attenB[timestep], 0, 5);
            } else {
                //SET_KERNEL_ARG_TEXTURE(kernelBP, d_attenIm[timestep], 0, 5);
            }
            if (DEBUG) mexPrint("backwardProjection: buffer 5 set (attenB)");
        }
        
        if (inputScalars.maskFP) {
            if (inputScalars.useBuffers) {
                int subset = 0;
                if (inputScalars.maskFPZ > 1) subset = osa_iter;
                SET_KERNEL_ARG_BUFFER(kernelBP, d_maskFPB[subset], 0, 6);
            } else {
                if (inputScalars.maskFPZ > 1) {
                    //kernelFP->setBuffer(d_maskFP3[osa_iter].get(), (NS::UInteger)0, 6);
                } else {
                    //kernelFP->setBuffer(d_maskFP[osa_iter].get(), (NS::UInteger)0, 6);
                }
            }
            if (DEBUG) mexPrint("backwardProjection: buffer 6 (maskFP) set");
        }

        if (inputScalars.maskBP) {
            if (inputScalars.useBuffers) {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_maskBPB, 0, 7);
            } else {
                if (inputScalars.maskBPZ > 1){
                    //status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP3);
                } else {
                    //status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP);
                }
            }
            if (DEBUG) mexPrint("backwardProjection: buffer 7 (maskBP) set");
        }

        if (compSens) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_xFull[0], 0, 8);
            SET_KERNEL_ARG_BUFFER(kernelBP, d_zFull[0], 0, 9);
        } else {
            if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[0][0], 0, 8);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[timestep][osa_iter], 0, 8);
            if (DEBUG) mexPrint("backwardProjection: buffer 8 set");

            if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
                SET_KERNEL_ARG_BUFFER(kernelBP, d_z[timestep][osa_iter], 0, 9);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_z[timestep][inputScalars.osa_iter0], 0, 9);
            if (DEBUG) mexPrint("backwardProjection: buffer 9 set");
        }

        if (inputScalars.normalization_correction) {
            if (compSens)
                SET_KERNEL_ARG_BUFFER(kernelBP, d_normFull[0], 0, 10);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_norm[osa_iter], 0, 10);
            if (DEBUG) mexPrint("backwardProjection: buffer 10 set (norm)");
        }

        if (inputScalars.scatter) {
            if (compSens) 
                SET_KERNEL_ARG_BUFFER(kernelBP, d_scatFull[0], 0, 11);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_scat[timestep][osa_iter], 0, 11);
            if (DEBUG) mexPrint("backwardProjection: buffer 11 set (scat)");
        }

        SET_KERNEL_ARG_BUFFER(kernelBP, d_Summ[uu], 0, 12);
        if (DEBUG) mexPrint("backwardProjection: buffer 12 set (Summ)");

        if ( (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0 ) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_xyindex[osa_iter], 0, 13);
            SET_KERNEL_ARG_BUFFER(kernelBP, d_zindex[osa_iter], 0, 14);
            if (DEBUG) mexPrint("backwardProjection: buffers 13 and 14 set (indices)");
        }

        if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
            if (!inputScalars.loadTOF) {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_trIndex[0][0], 0, 15);
                SET_KERNEL_ARG_BUFFER(kernelBP, d_axIndex[0][0], 0, 16);
            } else {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_trIndex[timestep][osa_iter], 0, 15);
                SET_KERNEL_ARG_BUFFER(kernelBP, d_axIndex[timestep][osa_iter], 0, 16);
            }
        }

        if (inputScalars.listmode > 0 && inputScalars.TOF) {
            if (!inputScalars.loadTOF) {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_TOFIndex[0][0], 0, 17);
            } else {
                SET_KERNEL_ARG_BUFFER(kernelBP, d_TOFIndex[timestep][osa_iter], 0, 17);
            }
        }

        if (inputScalars.raw)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_L[osa_iter], 0, 18);

        SET_KERNEL_ARG_BUFFER(kernelBP, d_output, 0, 19);
        if (DEBUG) mexPrint("backwardProjection: buffer 19 set (output)");

        SET_KERNEL_ARG_BUFFER(kernelBP, vec_opencl.d_rhs_os[uu], 0, 20);
        if (DEBUG) mexPrint("backwardProjection: buffer 20 set (rhs_os)");
    }
    if (inputScalars.CT && (inputScalars.BPType == 4 || inputScalars.BPType == 5)) { // TODO
        global[0] = inputScalars.Nx[ii] + erotusBP[0][ii];
        global[1] = inputScalars.Ny[ii] + erotusBP[1][ii];
        global[2] = inputScalars.Nz[ii];
        if (inputScalars.BPType == 4) {
            if (!inputScalars.largeDim) {
                if (!inputScalars.useHelical) {
                    global[2] = (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS;
                } else {
                    global[2] = (inputScalars.Nz[ii] + NVOXELSHELICAL - 1) / NVOXELSHELICAL;
                }
            }
        } else if (inputScalars.BPType == 5) {
            if (!inputScalars.pitch) {
                global[2] = (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5;
            }
        }
        if (!inputScalars.useBuffers) {
            NS::UInteger imX = inputScalars.nRowsD;
            NS::UInteger imY = inputScalars.nColsD;
            NS::UInteger imZ = length[osa_iter];
            if (inputScalars.BPType == 5) {
                imX++;
                imY++;
            }
            //cl::detail::size_t_array region = { imX, imY, imZ };
            //d_inputImage = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
            //status = CLCommandQueue[0].enqueueCopyBufferToImage(d_output, d_inputImage, 0, origin, region);
        }

        if (inputScalars.offset)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_T[osa_iter], 0, 1);
        if (inputScalars.useBuffers)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_output, 0, 2);
        //else
        //    status = kernelBP.setArg(kernelIndBPSubIter++, d_inputImage);

        if (inputScalars.CT && inputScalars.DSC > 0.f) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_angle, 0, 3);
        }
        SET_KERNEL_ARG_BUFFER(kernelBP, vec_opencl.d_rhs_os[uu], 0, 4);
        if (compSens)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_xFull[0], 0, 5);
        else
            if (!inputScalars.loadTOF && inputScalars.listmode > 0)
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[0][0], 0, 5);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[timestep][osa_iter], 0, 5);
        if (compSens)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_zFull[0], 0, 6);
        else
            SET_KERNEL_ARG_BUFFER(kernelBP, d_z[timestep][osa_iter], 0, 6);
        SET_KERNEL_ARG_BUFFER(kernelBP, d_Summ[ee], 0, 7);
        if (inputScalars.normalization_correction)
            SET_KERNEL_ARG_BUFFER(kernelBP, d_norm[osa_iter], 0, 8);
    }
    if (!inputScalars.CT && (inputScalars.BPType == 4 || inputScalars.BPType == 5)) { // TODO
        if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
            global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
            global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
            global[2] = length[osa_iter];
        } else if (inputScalars.listmode > 0 && compSens) {
            global[0] = static_cast<size_t>(inputScalars.det_per_ring + erotusSens[0]) / local[0];
            global[1] = (static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1]) / local[1];
            global[2] = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings);
        } else {
            erotus[0] = length[osa_iter] % local_size[0];

            if (erotus[0] > 0)
                erotus[0] = (local_size[0] - erotus[0]);
            global[0] = (length[osa_iter] + erotus[0]) / local[0];
            global[1] = 1;
            global[2] = 1;
        }
        
        SET_KERNEL_ARG_BUFFER(kernelBP, d_output, 0, 1);
        SET_KERNEL_ARG_BUFFER(kernelBP, vec_opencl.d_rhs_os[uu], 0, 2);
        if (compSens) {
            SET_KERNEL_ARG_BUFFER(kernelBP, d_xFull[0], 0, 3);
            SET_KERNEL_ARG_BUFFER(kernelBP, d_zFull[0], 0, 4);
        } else {
            if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[0][0], 0, 3);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_x[timestep][osa_iter], 0, 3);
            if (DEBUG) mexPrint("backwardProjection: buffer 3 set");

            if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
                SET_KERNEL_ARG_BUFFER(kernelBP, d_z[timestep][osa_iter], 0, 4);
            else if (inputScalars.indexBased && inputScalars.listmode > 0)
                SET_KERNEL_ARG_BUFFER(kernelBP, d_z[0][0], 0, 4);
            else
                SET_KERNEL_ARG_BUFFER(kernelBP, d_z[timestep][inputScalars.osa_iter0], 0, 4);
            if (DEBUG) mexPrint("backwardProjection: buffer 4 set");
        }
    }

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