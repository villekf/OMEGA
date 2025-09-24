#pragma once
#include "ProjectorClassMetal.mm"

template <typename T, typename C>
inline void reconstruction_metal(const float* z_det, const float* x,
    scalarStruct& inputScalars, Weighting& w_vec, RecMethods& MethodList,
    const int64_t* pituus, const char* header_directory,
    const float* meas /* =sinograms */, const float* im /* =image volume */, T* output, C* sensIm,
    const int type = 0, const int no_norm = 1, const float* rand = nullptr,
    const float* atten = nullptr, const float* norm = nullptr, const float* extraCorr = nullptr,
    const size_t size_gauss = 0, const uint32_t* xy_index = nullptr,
    const uint16_t* z_index = nullptr, const uint16_t* L = nullptr
) {
    const C tyyppi = (C)0;
	std::vector<int64_t> length(inputScalars.subsetsUsed); // Number of measurements in each subset

    for (uint32_t kk = 0; kk < inputScalars.subsetsUsed; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];
	uint64_t m_size = length[inputScalars.osa_iter0];

    ProjectorClass *proj = [[ProjectorClass alloc] init];
    NSInteger status = 0;

    status = [proj addProjector:[ScalarStructBox boxWithPointer:&inputScalars]
                                       weighting:[WeightingBox boxWithPointer:&w_vec]
                                          //method:[RecMethodsBox boxWithPointer:&methods]
                                headerDirectory:[NSString stringWithUTF8String:header_directory]
                                            type:-1];

    if (status != 0)
		return;
	proj.no_norm = no_norm;

    status = [proj createBuffers:[ScalarStructBox boxWithPointer:&inputScalars]
                                weighting:[WeightingBox boxWithPointer:&w_vec]
                                x:x
                            zDet:z_det
                        xyIndex:xy_index
                            zIndex:z_index
                                L:L
                            pituus:pituus
                            atten:atten
                            norm:norm
                        extraCorr:extraCorr
                            length:length.data()
                            type:type];

    if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[inputScalars.osa_iter0];
    
    NSUInteger bytesOutput = sizeof(float) * m_size * inputScalars.nBins;
    if (type == 1) { // Forward projection A*x
        status = [proj allocateOutput:bytesOutput]; // TODO replace with proj.d_output = ...
	}
	if (type == 2) { // Backward projection A'*y
        id<MTLBuffer> buf = [proj.device newBufferWithBytes:meas length:bytesOutput options:MTLResourceStorageModeShared];
        proj.d_output = buf;

        if (!proj.d_rhs_os) {
            proj.d_rhs_os = [NSMutableArray arrayWithCapacity:inputScalars.nMultiVolumes + 1];
        }

		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
            NSUInteger len = sizeof(float) * inputScalars.im_dim[ii];
            id<MTLBuffer> buf = [proj.device newBufferWithLength:len options:MTLResourceStorageModeShared];
            [proj.d_rhs_os addObject:buf];
            memset(buf.contents, 0, len);
			
            //proj.vec_opencl.d_rhs_os.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_WRITE_ONLY, sizeof(T) * inputScalars.im_dim[ii], NULL, &status));
			//status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.vec_opencl.d_rhs_os[ii], (T)0, 0, sizeof(T) * inputScalars.im_dim[ii]);

			if (proj.no_norm == 0) {} // Sensitivity image stuff (TODO)
			else {
                id<MTLBuffer> buf = [proj.device newBufferWithLength:sizeof(float) options:MTLResourceStorageModeShared]; // or Private
                [proj.d_Summ addObject:buf];
            }
		}
	}

    for (uint32_t iter = 0; iter < inputScalars.Niter; iter++) {
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
			uint32_t subIter = osa_iter;
			m_size = length[osa_iter];
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];
            if (type == 0) {}
            else if (type == 1) {
                size_t uu = 0;
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
                    NSInteger bytesIm = inputScalars.im_dim[ii] * sizeof(float);
                    const float *src = im + uu;
                    proj.d_im = [proj.device newBufferWithBytes:src length:bytesIm options:MTLResourceStorageModeShared];

                    status = [proj forwardProjection:[ScalarStructBox boxWithPointer:&inputScalars]
                                weighting:[WeightingBox boxWithPointer:&w_vec]
                                osaIter:osa_iter
                                length:length.data()
                                m_size:m_size
                                ii:ii
                                uu:0];
                    uu += inputScalars.im_dim[ii];
                }
            }
            if (type == 2 || type == 0) {
                for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					int uu = ii;
					if (type == 0) {}
                    else {
                        status = [proj backwardProjection:[ScalarStructBox boxWithPointer:&inputScalars]
                                weighting:[WeightingBox boxWithPointer:&w_vec]
                                osaIter:osa_iter
                                length:length.data()
                                m_size:m_size
                                computeSens:false
                                ii:ii
                                uu:uu];
                    }
                    if (type == 0) {}
                }
            }
        }
        if (type == 0) {}
    }
    
    if (type == 0) {}
    if (type == 1) { // Forward projection A*x
        size_t byteCount = sizeof(float) * m_size * inputScalars.nBins;
        void *src = proj.d_output.contents;            // id<MTLBuffer>
        memcpy(output, src, byteCount);                // copy into your CPU array
    }
    if (type == 2) { // Backward projection A'*y
        size_t uuElems = 0;
        for (int ii = 0; ii <= inputScalars.nMultiVolumes; ++ii) {
            const size_t elemCount = inputScalars.im_dim[ii];
            const size_t byteCount = elemCount * sizeof(float);

            const void *src = proj.d_rhs_os[ii].contents; // Shared buffer
            memcpy(output + uuElems, src, byteCount);

            uuElems += elemCount;
        }
    }
    
}