#include "AF_opencl_functions.hpp"

void forwardProjection(cl::CommandQueue& af_queue, cl_uint& kernelIndFPSubIter, const int64_t length, const scalarStruct& inputScalars,
	const Weighting& w_vec, const uint64_t& m_size, af::array& outputFP, const size_t local_size[], cl::Kernel& kernelFP) {
	cl_int status = CL_SUCCESS;
	size_t erotus[3];
	af::sync();
	af_queue.finish();
	cl::Buffer d_output = cl::Buffer(*outputFP.device<cl_mem>(), true);
	//cl::Buffer d_output = cl::Buffer(*outputFP.device<cl_mem>(), true);
	//d_output = cl::Buffer(af_context, CL_MEM_READ_WRITE, sizeof(cl_float) * m_size, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	//const size_t global_size = length[i] + erotus;
	size_t global_size = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;

	cl::NDRange local(local_size[0]);
	cl::NDRange global(global_size);
	global = { w_vec.size_x + erotus[0], w_vec.size_y + erotus[1], static_cast<size_t>(length) };
	local = { local_size[0] , local_size[1] };

	if (DEBUG) {
		mexPrintf("global[0] = %u\n", global[0]);
		mexPrintf("local[0] = %u\n", local[0]);
		mexPrintf("local[1] = %u\n", local[1]);
		mexPrintf("global[1] = %u\n", global[1]);
		mexPrintf("global[2] = %u\n", global[2]);
		mexPrintf("erotus[0] = %u\n", erotus[0]);
		mexPrintf("erotus[1] = %u\n", erotus[1]);
		mexPrintf("global.dimensions() = %u\n", global.dimensions());
		mexPrintf("local.dimensions() = %u\n", local.dimensions());
		mexPrintf("m_size = %u\n", m_size);
		mexPrintf("size_x = %u\n", w_vec.size_x);
		mexPrintf("size_y = %u\n", w_vec.size_y);
		mexPrintf("length[osa_iter] = %u\n", length);
		mexPrintf("st = %u\n", st);
		mexPrintf("listmode = %u\n", listmode);
		mexPrintf("maskBP = %u\n", inputScalars.maskBP);
		mexPrintf("st = %u\n", st);
		mexPrintf("im_dim = %u\n", im_dim);
		mexPrintf("no_norm = %u\n", no_norm);
		mexPrintf("NVOXELS = %u\n", NVOXELS);
		mexEvalString("pause(.0001);");
		//mexEvalString("pause(2);");
	}

	af::sync();
	af_queue.finish();
	if (inputScalars.projector_type == 4) {
		status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, d_output);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else if (inputScalars.projector_type == 5) {
		status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os_int);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = kernelFP.setArg(kernelIndFPSubIter++, d_output);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (inputScalars.meanFP) {

		}
	}
	else if (inputScalars.projector_type == 1 || inputScalars.projector_type == 14) {

	}
	kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
	if (iter > 0 || !MethodList.LSQR)
		status = af_queue.enqueueNDRangeKernel(kernelFP, cl::NDRange(), global, local, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		//outputFP.unlock();
		if (compute_norm_matrix == 1u) {
			Summ[0].unlock();
		}
		else {
			if (no_norm == 0u) {
				Summ[osa_iter].unlock();
			}
			else
				apu_sum.unlock();
		}
		if (inputScalars.use_psf)
			vec.im_os_blurred.unlock();
		else
			vec.im_os.unlock();
		vec.rhs_os.unlock();
		outputFP.unlock();
		return;
	}
	else if (DEBUG) {
		mexPrintf("Forward projection kernel launched successfully\n");
		mexEvalString("pause(.0001);");
	}
	status = af_queue.finish();
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Queue finish failed after FP kernel\n");
		mexEvalString("pause(.0001);");
	}

	outputFP.unlock();
	af::sync();
}