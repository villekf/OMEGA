#import "ProjectorClassMetal.h"

// ---------- Box impls ----------
@implementation ScalarStructBox { void *_p; }
+ (instancetype)boxWithPointer:(void *)ptr { ScalarStructBox *b=[self new]; b->_p=ptr; return b; }
- (void *)ptr { return _p; }
@end
@implementation WeightingBox { void *_p; }
+ (instancetype)boxWithPointer:(void *)ptr { WeightingBox *b=[self new]; b->_p=ptr; return b; }
- (void *)ptr { return _p; }
@end
/*@implementation RecMethodsBox { void *_p; }
+ (instancetype)boxWithPointer:(void *)ptr { RecMethodsBox *b=[self new]; b->_p=ptr; return b; }
- (void *)ptr { return _p; }
@end*/

typedef struct {
    float global_factor;
	float d_epps;
	uint d_size_x;
	uint d_det_per_ring;
	float sigma_x;
	float coneOfResponseStdCoeffA;
    float coneOfResponseStdCoeffB;
    float coneOfResponseStdCoeffC;
	float crystalSizeX;
	float crystalSizeY;
	float orthWidth;
	float bmin;
	float bmax;
	float Vmax;
	uint d_sizey;
    long d_nProjections;
    uint d_Nx;
	uint d_Ny;
	uint d_Nz;
	float d_dx;
	float d_dy;
	float d_dz;
	float bx;
	float by;
	float bz;
	float d_bmaxx;
	float d_bmaxy;
	float d_bmaxz;
    unsigned char no_norm;
	unsigned long m_size;
	uint currentSubset;
	int aa;
} ParamsConst;

// -------- Projector class ---------
@interface ProjectorClass () {
    id<MTLLibrary> _libFP;
    id<MTLLibrary> _libBP;
    id<MTLFunction> _fnFP;
    id<MTLFunction> _fnBP;
    id<MTLComputePipelineState> _psoFP;
    id<MTLComputePipelineState> _psoBP;
    id<MTLCommandQueue> _queueFP;
    id<MTLCommandQueue> _queueBP;
    id<MTLCommandBuffer> _commandBufferFP;
    id<MTLCommandBuffer> _commandBufferBP;
    id<MTLComputeCommandEncoder> _encFP;
    id<MTLComputeCommandEncoder> _encBP;

    size_t _erotus[3];
    std::vector<std::vector<size_t>> _erotusBP;
    NSUInteger _localX;
    NSUInteger _localY;
    NSUInteger _localZ;
    NSUInteger _globalX;
    NSUInteger _globalY;
    NSUInteger _globalZ;
    std::vector<float> bx, by, bz, dx, dy, dz, bmaxx, bmaxy, bmaxz;
	std::vector<NSInteger> d_Nx, d_Ny, d_Nz;
    // Buffers
    id<MTLBuffer> _d_V;
    id<MTLBuffer> _d_rayShiftsDetector;
    id<MTLBuffer> _d_rayShiftsSource;
    std::vector<__strong id<MTLBuffer>> _d_x;
    std::vector<__strong id<MTLBuffer>> _d_z;
    std::vector<__strong id<MTLBuffer>> _d_xyindex;
    std::vector<__strong id<MTLBuffer>> _d_zindex;
    std::vector<__strong id<MTLBuffer>> _d_L;
}
@end

@implementation ProjectorClass

- (instancetype)init {
    if ((self = [super init])) {

    }
    return self;
}

static NSString *ReadUTF8(NSString *path, NSError **err) {
    return [NSString stringWithContentsOfFile:path encoding:NSUTF8StringEncoding error:err];
}

- (NSInteger)allocateOutput:(NSInteger)nBytes // Allocate and fill with zeroes
{
    _d_output = [_device newBufferWithLength:nBytes options:MTLResourceStorageModeShared];
    memset([_d_output contents], 0, nBytes); // if Shared/Managed
    [_d_output didModifyRange:NSMakeRange(0, nBytes)];
    return 0;
}

- (NSInteger)setOutput:(NSInteger)nBytes // Allocate and fill with input data
                    data:(const float*)d_output
{
    _d_output = [_device newBufferWithBytes:d_output length:nBytes options:MTLResourceStorageModeShared];
    return 0;
}

- (NSInteger)setImage:(NSInteger)nBytes
                    data:(const float*)d_im
{
    _d_im = [_device newBufferWithBytes:d_im length:nBytes options:MTLResourceStorageModeShared];
    return 0;
}

- (NSInteger)addProjector:(ScalarStructBox *)inputScalarsBox
                           weighting:(WeightingBox *)wVec
                              //method:(RecMethodsBox *)methodList
                    headerDirectory:(NSString *)headerDirectory
                                type:(NSInteger)type
{ 
    // Unbox to C++ refs (no copies).
    auto &inputScalars = *static_cast<scalarStruct *>(inputScalarsBox.ptr);
    auto &w_vec = *static_cast<Weighting *>(wVec.ptr);
    //auto const &methods = *static_cast<RecMethods const *>(methodList.ptr);
    mexPrintf("test1.1\n");
    // ---- Metal setup ----
    _device = MTLCreateSystemDefaultDevice();
    if (!_device) {
        NSLog(@"Metal device not found.");
        return -1;
    }
    
    // ---- Read kernel/header files from headerDirectory ----
    // TODO: projector types, different FP/BP programs
    NSString *base = [headerDirectory stringByStandardizingPath];
    NSString *headerPath = [base stringByAppendingPathComponent:@"general_opencl_functions.h"];
    NSString *kernelPath = [base stringByAppendingPathComponent:@"projectorType123.cl"];
    NSError *err = nil; 
    NSString *headerSrc = ReadUTF8(headerPath, &err);
    if (!headerSrc) { NSLog(@"Failed to read header: %@", err.localizedDescription); return -2; }
    NSString *kernelSrc = ReadUTF8(kernelPath, &err);
    if (!kernelSrc) { NSLog(@"Failed to read kernel: %@", err.localizedDescription); return -3; }
    NSString *src = [NSString stringWithFormat:@"%@\n\n%@", headerSrc, kernelSrc];

    // https://developer.apple.com/documentation/metal/mtlcompileoptions/preprocessormacros
    MTLCompileOptions *optsFP = [MTLCompileOptions new]; // TODO preprocessor macros
    MTLCompileOptions *optsBP = [MTLCompileOptions new]; // TODO preprocessor macros

    NSMutableDictionary *macros = [@{
        @"METAL": @""
    } mutableCopy];

    // TODO refactor to a function
    const bool siddonVal = (inputScalars.FPType == 1 || inputScalars.BPType == 1 || inputScalars.FPType == 4 || inputScalars.BPType == 4) ? true : false;

    if (inputScalars.SPECT) {
        macros[@"SPECT"] = @"";
    }
    macros[@"NBINS"] = @(inputScalars.nBins); 
    if ((siddonVal && ((inputScalars.n_rays * inputScalars.n_rays3D) > 1)) || inputScalars.SPECT) {;
        macros[@"N_RAYS"] = @(inputScalars.n_rays * inputScalars.n_rays3D); 
        macros[@"N_RAYS2D"] = @(inputScalars.n_rays); 
        macros[@"N_RAYS3D"] = @(inputScalars.n_rays3D); 
    }
    if (inputScalars.pitch)
        macros[@"PITCH"] = @"";
    if (((inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
        macros[@"SUBSETS"] = @"";

    if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
        macros[@"STYPE"] = @(inputScalars.subsetType); 
        macros[@"NSUBSETS"] = @(inputScalars.subsets);
    }

    NSMutableDictionary *macrosFP = macros;
    NSMutableDictionary *macrosBP = macros;
    
    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        macrosFP[@"SIDDON"] = @"";
        macrosFP[@"ATOMICF"] = @"";
        macrosBP[@"SIDDON"] = @"";
        macrosBP[@"ATOMICF"] = @"";
        if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)
            macrosFP[@"FP"] = @"";
        if (inputScalars.FPType == 3)
            macrosFP[@"VOL"] = @"";
        if (inputScalars.FPType == 2 || inputScalars.FPType == 3)
            macrosFP[@"ORTH"] = @"";

        if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3)
            macrosBP[@"BP"] = @"";
        if (inputScalars.BPType == 3)
            macrosBP[@"VOL"] = @"";
        if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
            macrosBP[@"ORTH"] = @"";
    }


    // END todo
    optsFP.preprocessorMacros = macrosFP;
    optsBP.preprocessorMacros = macrosBP;

    _libFP = [_device newLibraryWithSource:src options:optsFP error:&err];
    
    if (!_libFP) {
        NSString *errStr = err.localizedDescription;
        mexPrintf(errStr.UTF8String);
        NSLog(@"Metal shader compile failed: %@", err.localizedDescription);
        return -2;
    }
    
    _libBP = [_device newLibraryWithSource:src options:optsBP error:&err];
    if (!_libBP) {
        NSString *errStr = err.localizedDescription;
        mexPrintf(errStr.UTF8String);
        NSLog(@"Metal shader compile failed: %@", err.localizedDescription);
        return -2;
    }
    
    _fnFP = [_libFP newFunctionWithName:@"projectorType123"];
    _fnBP = [_libBP newFunctionWithName:@"projectorType123"];

    _psoFP = [_device newComputePipelineStateWithFunction:_fnFP error:&err];
    _psoBP = [_device newComputePipelineStateWithFunction:_fnBP error:&err];

    _queueFP = [_device newCommandQueue];
    _queueBP = [_device newCommandQueue];
    
    _localX = 64ULL;
    _localY = 1ULL;
    _localZ = 1ULL;
    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || (inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT)))
        _localX = 128ULL;
    if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
        if (inputScalars.nColsD > 1 && !(inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT))) {
            _localX = 16ULL;
            _localY = 16ULL;
        }
    }

    if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
        _erotus[0] = inputScalars.nRowsD % _localX;
        //if (inputScalars.FPType == 5)
        //    erotus[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP) % local_size[1];
        //else
            _erotus[1] = inputScalars.nColsD % _localY;
        if (_erotus[1] > 0)
            _erotus[1] = (_localY - _erotus[1]);
        if (_erotus[0] > 0)
            _erotus[0] = (_localX - _erotus[0]);
    }

    _erotusBP.resize(2);
    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
        _erotusBP[0].emplace_back(inputScalars.Nx[ii] % _localX);
        _erotusBP[1].emplace_back(inputScalars.Ny[ii] % _localY);
        if (_erotusBP[0][ii] > 0)
            _erotusBP[0][ii] = (_localX - _erotusBP[0][ii]);
        if (_erotusBP[1][ii] > 0)
            _erotusBP[1][ii] = (_localY - _erotusBP[1][ii]);
    }

    
    bx.resize(inputScalars.nMultiVolumes + 1);
    by.resize(inputScalars.nMultiVolumes + 1);
    bz.resize(inputScalars.nMultiVolumes + 1);
    dx.resize(inputScalars.nMultiVolumes + 1);
    dy.resize(inputScalars.nMultiVolumes + 1);
    dz.resize(inputScalars.nMultiVolumes + 1);
    d_Nx.resize(inputScalars.nMultiVolumes + 1);
    d_Ny.resize(inputScalars.nMultiVolumes + 1);
    d_Nz.resize(inputScalars.nMultiVolumes + 1);
    bmaxx.resize(inputScalars.nMultiVolumes + 1);
    bmaxy.resize(inputScalars.nMultiVolumes + 1);
    bmaxz.resize(inputScalars.nMultiVolumes + 1);

    
    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
        bx[ii] = inputScalars.bx[ii];
        by[ii] = inputScalars.by[ii];
        bz[ii] = inputScalars.bz[ii];
        dx[ii] = inputScalars.dx[ii];
        dy[ii] = inputScalars.dy[ii];
        dz[ii] = inputScalars.dz[ii];
        d_Nx[ii] = static_cast<NSInteger>(inputScalars.Nx[ii]);
        d_Ny[ii] = static_cast<NSInteger>(inputScalars.Ny[ii]);
        d_Nz[ii] = static_cast<NSInteger>(inputScalars.Nz[ii]);
        bmaxx[ii] = static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii];
        bmaxy[ii] = static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii];
        bmaxz[ii] = static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii];
    }
    //region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
    
    return 0; 
}

- (NSInteger)createBuffers:(ScalarStructBox *)inputScalarsBox
                           weighting:(WeightingBox *)wVec
                           x:(const float *)x
                       zDet:(const float *)z_det
                   xyIndex:(const uint32_t *)xy_index
                    zIndex:(const uint16_t *)z_index
                           L:(const uint16_t *)L
                    pituus:(const int64_t *)pituus
                     atten:(const float *)atten
                      norm:(const float *)norm
                 extraCorr:(const float *)extraCorr
                    length:(const int64_t *)length
                    type:(NSUInteger)type
{
    // Unbox to C++ refs (no copies).
    auto &inputScalars = *static_cast<scalarStruct *>(inputScalarsBox.ptr);
    auto &w_vec = *static_cast<Weighting *>(wVec.ptr);

    // Resize necessary buffer vectors
    if (inputScalars.raw)
        _d_L.resize(inputScalars.subsetsUsed);
    if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1) {
        _d_xyindex.resize(inputScalars.subsetsUsed);
        _d_zindex.resize(inputScalars.subsetsUsed);
    }
    if (inputScalars.projector_type != 6) {
        _d_x.resize(inputScalars.subsetsUsed);
        _d_z.resize(inputScalars.subsetsUsed);
    }

    // --- Buffers (from createAndWriteBuffers) --- TODO
    if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        NSUInteger bytesV = sizeof(float) * inputScalars.size_V;
        _d_V = [_device newBufferWithBytes:inputScalars.V length:bytesV options:MTLResourceStorageModeShared];
    }
    if (inputScalars.SPECT) {
        NSUInteger bytesRayShifts = sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections;
        _d_rayShiftsDetector = [_device newBufferWithBytes:w_vec.rayShiftsDetector length:bytesRayShifts options:MTLResourceStorageModeShared];
        _d_rayShiftsSource = [_device newBufferWithBytes:w_vec.rayShiftsSource length:bytesRayShifts options:MTLResourceStorageModeShared];
    }
    for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
        NSUInteger bytesZ = sizeof(float) * length[kk] * 2;
        NSUInteger bytesX = sizeof(float) * length[kk] * 6;
        _d_z[kk] = [_device newBufferWithBytes:&z_det[pituus[kk] * 2] length:bytesZ options:MTLResourceStorageModeShared];
        _d_x[kk] = [_device newBufferWithBytes:&x[pituus[kk] * 6] length:bytesX options:MTLResourceStorageModeShared];

        // Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
        if (inputScalars.raw && inputScalars.listmode != 1) {
            NSUInteger bytesL = sizeof(uint16_t) * length[kk] * 2;
            _d_L[kk] = [_device newBufferWithBytes:&L[pituus[kk] * 2] length:bytesL options:MTLResourceStorageModeShared];
        }
    }

    return 0;
}

- (NSInteger)forwardProjection:(ScalarStructBox *)inputScalarsBox
                           weighting:(WeightingBox *)wVec
                           osaIter:(uint32_t)osa_iter
                         length:(const int64_t *)length
                         m_size:(uint64_t)m_size
                              ii:(int32_t)ii
                              uu:(int)uu
{
    mexPrintf("init forwardProjection\n");
    // Unbox to C++ refs (no copies).
    auto &inputScalars = *static_cast<scalarStruct *>(inputScalarsBox.ptr);
    auto &w_vec = *static_cast<Weighting *>(wVec.ptr);

    _globalX = inputScalars.nRowsD + _erotus[0];
    _globalY = inputScalars.nColsD + _erotus[1];
    _globalZ = static_cast<size_t>(length[osa_iter]);
    mexPrintf("fp1.1\n");

    _commandBufferFP = [_queueFP commandBuffer];
    _encFP = [_commandBufferFP computeCommandEncoder];
    [_encFP setComputePipelineState:_psoFP];
    mexPrintf("fp1.2\n");

    ParamsConst params = {}; // Move params to struct
    params.global_factor = inputScalars.global_factor;
	params.d_epps = inputScalars.epps;
	params.d_size_x = inputScalars.nRowsD;
	params.d_det_per_ring = inputScalars.det_per_ring;
	params.sigma_x = inputScalars.sigma_x;
	params.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
    params.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
    params.coneOfResponseStdCoeffC = inputScalars.coneOfResponseStdCoeffC;
	params.crystalSizeX = w_vec.dPitchX;
	params.crystalSizeY = w_vec.dPitchY;
    if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        if (inputScalars.FPType == 2)
            params.orthWidth = inputScalars.tube_width;
        if (inputScalars.FPType == 3)
            params.orthWidth = inputScalars.cylRadiusProj3;
        params.bmin = inputScalars.bmin;
        params.bmax = inputScalars.bmax;
        params.Vmax = inputScalars.Vmax;
    }
	params.d_sizey = inputScalars.nColsD;
    params.d_nProjections = length[osa_iter];
    params.d_Nx = d_Nx[ii];
    params.d_Ny = d_Ny[ii];
    params.d_Nz = d_Nz[ii];
    params.d_dx = dx[ii];
    params.d_dy = dy[ii];
    params.d_dz = dz[ii];
    params.bx = bx[ii];
    params.by = by[ii];
    params.bz = bz[ii];
    params.d_bmaxx = bmaxx[ii];
    params.d_bmaxy = bmaxy[ii];
    params.d_bmaxz = bmaxz[ii];
    params.no_norm = _no_norm;
	params.m_size = m_size;
	params.currentSubset = osa_iter;
	params.aa = ii;
    mexPrintf("fp1.3\n");

    if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        [_encFP setBytes:&params length:sizeof(params) atIndex:0];
        if (inputScalars.SPECT) {
            [_encFP setBuffer:_d_rayShiftsDetector offset:0 atIndex:1];
            [_encFP setBuffer:_d_rayShiftsSource offset:0 atIndex:2];
        }
        if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
            [_encFP setBuffer:_d_V offset:0 atIndex:4];
        }
    }
    [_encFP setBuffer:_d_x[osa_iter] offset:0 atIndex:5];
    [_encFP setBuffer:_d_z[osa_iter] offset:0 atIndex:6];
    mexPrintf("fp1.4\n");
    [_encFP setBuffer:_d_Summ[uu] offset:0 atIndex:7];
    mexPrintf("fp1.5\n");
    if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
        //[_enc setBuffer:_d_xyindex[osa_iter] offset:0 atIndex:_kernelIndFP++];
        //[_enc setBuffer:_d_zindex[osa_iter] offset:0 atIndex:_kernelIndFP++];
    }

    //if (inputScalars.raw)
        //kernelFP.setArg(kernelIndFPSubIter++, d_L[osa_iter]);
    [_encFP setBuffer:_d_im offset:0 atIndex:8];    
    [_encFP setBuffer:_d_output offset:0 atIndex:9];

    // --- Build Metal sizes that exactly match CL's parameters ---
    MTLSize threadsPerThreadgroup = MTLSizeMake(_localX, _localY, _localZ);

    // Validate the local size against the pipeline’s limit
    //NSUInteger maxTG = _pso.maxTotalThreadsPerThreadgroup;
    /*NSAssert(_localX * _localY * _localZ <= maxTG,
            @"Local work-group size (%lu) exceeds Metal's maxTotalThreadsPerThreadgroup (%lu). "
            @"Choose a smaller local size to match CL.",
            (unsigned long)(_localX*_localY*_localZ), (unsigned long)maxTG);*/

    // In OpenCL, global must be a multiple of local (with zero offset).
    // Mirror that here and compute the number of threadgroups per grid.
    /*NSAssert(_globalX % _localX == 0 && _globalY % _localY == 0 && _globalZ % _localZ == 0,
            @"global must be an exact multiple of local in each dimension (as in OpenCL).");*/

    MTLSize threadgroupsPerGrid = MTLSizeMake(_globalX / _localX, _globalY / _localY, _globalZ / _localZ);
    
    [_encFP dispatchThreadgroups:threadgroupsPerGrid threadsPerThreadgroup:threadsPerThreadgroup];
    
    [_encFP endEncoding];
    mexPrintf("fp1.6\n");
    [_commandBufferFP commit];
    mexPrintf("fp1.7\n");
    [_commandBufferFP waitUntilCompleted];

    return 0;
}

- (NSInteger)backwardProjection:(ScalarStructBox *)inputScalarsBox
                           weighting:(WeightingBox *)wVec
                           osaIter:(uint32_t)osa_iter
                         length:(const int64_t *)length
                         m_size:(uint64_t)m_size
                         computeSens:(bool)computeSens
                              ii:(int32_t)ii
                              uu:(int)uu
{
    mexPrintf("init backwardProjection\n");
    // Unbox to C++ refs (no copies).
    auto &inputScalars = *static_cast<scalarStruct *>(inputScalarsBox.ptr);
    auto &w_vec = *static_cast<Weighting *>(wVec.ptr);

    _globalX = inputScalars.nRowsD + _erotus[0];
    _globalY = inputScalars.nColsD + _erotus[1];
    _globalZ = static_cast<size_t>(length[osa_iter]);
    mexPrintf("bp1.1\n");

    _commandBufferBP = [_queueBP commandBuffer];
    mexPrintf("bp1.11\n");
    _encBP = [_commandBufferBP computeCommandEncoder];
    mexPrintf("bp1.12\n");
    [_encBP setComputePipelineState:_psoBP];
    mexPrintf("bp1.2\n");

    ParamsConst params = {}; // Move params to struct
    params.global_factor = inputScalars.global_factor;
	params.d_epps = inputScalars.epps;
	params.d_size_x = inputScalars.nRowsD;
	params.d_det_per_ring = inputScalars.det_per_ring;
	params.sigma_x = inputScalars.sigma_x;
	params.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
    params.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
    params.coneOfResponseStdCoeffC = inputScalars.coneOfResponseStdCoeffC;
	params.crystalSizeX = w_vec.dPitchX;
	params.crystalSizeY = w_vec.dPitchY;
    if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
        if (inputScalars.FPType == 2)
            params.orthWidth = inputScalars.tube_width;
        if (inputScalars.FPType == 3)
            params.orthWidth = inputScalars.cylRadiusProj3;
        params.bmin = inputScalars.bmin;
        params.bmax = inputScalars.bmax;
        params.Vmax = inputScalars.Vmax;
    }
	params.d_sizey = inputScalars.nColsD;
    params.d_nProjections = length[osa_iter];
    params.d_Nx = d_Nx[ii];
    params.d_Ny = d_Ny[ii];
    params.d_Nz = d_Nz[ii];
    params.d_dx = dx[ii];
    params.d_dy = dy[ii];
    params.d_dz = dz[ii];
    params.bx = bx[ii];
    params.by = by[ii];
    params.bz = bz[ii];
    params.d_bmaxx = bmaxx[ii];
    params.d_bmaxy = bmaxy[ii];
    params.d_bmaxz = bmaxz[ii];
    params.no_norm = _no_norm;
	params.m_size = m_size;
	params.currentSubset = osa_iter;
	params.aa = ii;
    mexPrintf("bp1.3\n");

    if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
        [_encBP setBytes:&params length:sizeof(params) atIndex:0];
        if (inputScalars.SPECT) {
            [_encBP setBuffer:_d_rayShiftsDetector offset:0 atIndex:1];
            [_encBP setBuffer:_d_rayShiftsSource offset:0 atIndex:2];
        }
        if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
            [_encBP setBuffer:_d_V offset:0 atIndex:4];
        }
    }
    [_encBP setBuffer:_d_x[osa_iter] offset:0 atIndex:5];
    [_encBP setBuffer:_d_z[osa_iter] offset:0 atIndex:6];
    mexPrintf("bp1.4\n");
    [_encBP setBuffer:_d_Summ[uu] offset:0 atIndex:7];
    mexPrintf("bp1.5\n");
    if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
        //[_enc setBuffer:_d_xyindex[osa_iter] offset:0 atIndex:_kernelIndFP++];
        //[_enc setBuffer:_d_zindex[osa_iter] offset:0 atIndex:_kernelIndFP++];
    }

    //if (inputScalars.raw)
        //kernelFP.setArg(kernelIndFPSubIter++, d_L[osa_iter]); 
    [_encBP setBuffer:_d_output offset:0 atIndex:8];
    [_encBP setBuffer:_d_rhs_os[uu] offset:0 atIndex:9];   

    // --- Build Metal sizes that exactly match CL's parameters ---
    MTLSize threadsPerThreadgroup = MTLSizeMake(_localX, _localY, _localZ);
    MTLSize threadgroupsPerGrid = MTLSizeMake(_globalX / _localX, _globalY / _localY, _globalZ / _localZ);
    
    [_encBP dispatchThreadgroups:threadgroupsPerGrid threadsPerThreadgroup:threadsPerThreadgroup];
    
    [_encBP endEncoding];
    mexPrintf("bp1.6\n");
    [_commandBufferBP commit];
    mexPrintf("bp1.7\n");
    [_commandBufferBP waitUntilCompleted];

    return 0;
}

@end