/**************************************************************************
* Implements either the improved Siddon's algorithm, orthogonal 
* distance based projector or volume-based projector for OMEGA.
*
* Implementation 1 supports only Siddon's algorithm and outputs the sparse
* system matrix. Implementation 4 supports all projectors and is computed
* matrix-free.
* 
* This file contains the oct-functions for both the implementation type 1
* and type 4.
* 
* A precomputed vector that has the number of voxels a LOR/ray passes, as 
* well as the number of voxels that have been traversed in previous 
* LORs/rays is required if options.precompute_lor = true, but is optional 
* otherwise. If the precompute_lor options is set to false and implementation 
* type 4 has been selected then all the measurements are investigated, but 
* only those intercepting the FOV will be included.
* 
* Both implementations use OpenMP for parallellization (if available, 
* otherwise the code will be sequential with no parallelization).
* 
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include "projector_functions.h"
 
using namespace std;


DEFUN_DLD(projector_oct, prhs, nargout, "projector_oct") {

	int ind = 0;
	// Load the input arguments
	// Image size in y-direction
	const uint32_t Ny = prhs(ind).uint32_scalar_value();
	ind++;

	// Image size in x-direction
	const uint32_t Nx = prhs(ind).uint32_scalar_value();
	ind++;

	// Image size in z-direction
	const uint32_t Nz = prhs(ind).uint32_scalar_value();
	ind++;

	// Distance between adjacent pixels in x-direction
	const double d = prhs(ind).scalar_value();
	ind++;

	// Distance between adjacent pixels in z-direction
	const double dz = prhs(ind).scalar_value();
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	const double by = prhs(ind).scalar_value();
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	const double bx = prhs(ind).scalar_value();
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	const double bz = prhs(ind).scalar_value();
	ind++;

	// Coordinates of the detectors in z-direction
	const NDArray z_det_ = prhs(ind).array_value();

	vector<double> z_det_vec(z_det_.numel(), 0.);
	std::copy_n(z_det_.data(), z_det_.numel(), z_det_vec.begin());
	ind++;

	const double* z_det = z_det_.fortran_vec();

	// Coordinates of the detectors in x-direction
	const NDArray x_ = prhs(ind).array_value();
	ind++;

	const double* x = x_.fortran_vec();

	// Coordinates of the detectors in y-direction
	const NDArray y_ = prhs(ind).array_value();
	ind++;

	const double* y = y_.fortran_vec();

	// Distance between adjacent pixels in y-direction
	const double dy = prhs(ind).scalar_value();
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in y-direction
	const NDArray yy = prhs(ind).array_value();

	// from array to std::vector
	vector<double> yy_vec(yy.numel(), 0.);
	std::copy_n(yy.data(), yy.numel(), yy_vec.begin());
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in x-direction
	const NDArray xx = prhs(ind).array_value();

	vector<double> xx_vec(yy.numel(), 0.);
	std::copy_n(xx.data(), xx.numel(), xx_vec.begin());
	ind++;

	// Number of sinograms used
	const uint32_t NSinos = prhs(ind).uint32_scalar_value();
	ind++;

	// Number of slices included
	const uint32_t NSlices = prhs(ind).uint32_scalar_value();
	ind++;

	// Number of detector indices
	const uint32_t size_x = prhs(ind).uint32_scalar_value();
	ind++;

	// Maximum value of the z-direction detector coordinates
	const double zmax = prhs(ind).scalar_value();
	ind++;

	// attenuation values (attenuation images)
	const NDArray atten_ = prhs(ind).array_value();
	ind++;

	const double* atten = atten_.fortran_vec();

	// Normalization coefficients
	const NDArray norm_coef_ = prhs(ind).array_value();
	ind++;

	const double* norm_coef = norm_coef_.fortran_vec();

	// Randoms
	const NDArray randoms_ = prhs(ind).array_value();
	ind++;

	const double* randoms = randoms_.fortran_vec();

	// Number of measurements/LORs
	const uint32_t pituus = prhs(ind).uint32_scalar_value();
	ind++;

	// Is the attenuation correction included
	const bool attenuation_correction = prhs(ind).bool_value();
	ind++;

	// Is the normalization correction included
	const bool normalization = prhs(ind).bool_value();
	ind++;

	// Is the randoms/scatter correction included
	const bool randoms_correction = prhs(ind).bool_value();
	ind++;

	// Is scatter correction included as multiplication (system matrix)
	const bool scatter = prhs(ind).bool_value();
	ind++;

	// Scatter data
	const NDArray scatter_coef_ = prhs(ind).array_value();
	ind++;

	const double* scatter_coef = scatter_coef_.fortran_vec();

	// Global correction factor
	const double global_factor = prhs(ind).scalar_value();
	ind++;

	// Number of voxels the current LOR/ray traverses, precomputed data ONLY
	uint16NDArray lor1_ = prhs(ind).uint16_array_value();
	ind++;

	const uint16_t* lor1 = reinterpret_cast<uint16_t*>(lor1_.fortran_vec());

	// For sinogram data, the indices of the detectors corresponding to the current sinogram bin
	uint32NDArray xy_index_ = prhs(ind).uint32_array_value();
	ind++;

	const uint32_t* xy_index = reinterpret_cast<uint32_t*>(xy_index_.fortran_vec());

	// Same as above, but for z-direction
	uint16NDArray z_index_ = prhs(ind).uint16_array_value();
	ind++;

	const uint16_t* z_index = reinterpret_cast<uint16_t*>(z_index_.fortran_vec());

	// Total number of sinograms
	const uint32_t TotSinos = prhs(ind).uint32_scalar_value();
	ind++;

	// Detector pair numbers, for raw list-mode data
	uint16NDArray L_ = prhs(ind).uint16_array_value();
	const size_t numRows = L_.numel();
	ind++;

	const uint16_t* L = reinterpret_cast<uint16_t*>(L_.fortran_vec());

	// Location (ring numbers) of pseudo rings, if present
	uint32NDArray pseudos_ = prhs(ind).uint32_array_value();
	const uint32_t pRows = pseudos_.numel();
	ind++;

	const uint32_t* pseudos = reinterpret_cast<uint32_t*>(pseudos_.fortran_vec());

	// Number of detectors per ring
	const uint32_t det_per_ring = prhs(ind).uint32_scalar_value();
	ind++;

	// Is TOF data used?
	const bool TOF = prhs(ind).bool_scalar_value();
	ind++;

	// Size of single TOF-subset
	const int64_t TOFSize = prhs(ind).int64_t_scalar_value();
	ind++;

	// Variance of the Gaussian TOF
	const double sigma_x = prhs(ind).scalar_value();
	ind++;

	// Centers of the TOF-bins
	const double* TOFCenter = prhs(ind).array_value();
	ind++;

	// Index offset for TOF subsets
	const int64_t nBins = prhs(ind).int64_t_scalar_value();
	ind++;

	const uint32_t dec_v = prhs(ind).uint32_scalar_value();
	ind++;

	// Are status messages displayed
	const bool verbose = prhs(ind).bool_value();
	ind++;

	// Is raw list-mode data used
	const bool raw = prhs(ind).bool_value();
	ind++;

	const uint32_t type = prhs(ind).uint32_scalar_value();
	ind++;

	// Number of measurements/LORs
	int64_t loop_var_par = pituus;

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	octave_value_list retval(nargout);

	// Implementation 1, with precomputed_lor = true
	if (type == 0u) {

		// Total number of voxels traversed by the LORs at the specific LOR
		// e.g. if the first LOR traverses through 10 voxels then at the second lor the value is 10
		const uint64NDArray lor2_ = prhs(ind).uint64_array_value();
		ind++;

		const uint64_t* lor2 = reinterpret_cast<const uint64_t*>(lor2_.fortran_vec());

		// Total number of non-zero values
		const uint64_t summa = prhs(ind).uint64_scalar_value();
		ind++;

		// Is this a transmission image reconstruction (true) or emission image reconstruction (false)
		const bool attenuation_phase = prhs(ind).bool_value();
		ind++;

		const uint32_t projector_type = prhs(ind).uint32_scalar_value();
		ind++;

		// output sizes
		const octave_idx_type N = static_cast<octave_idx_type>(Nx) * static_cast<octave_idx_type>(Ny) * static_cast<octave_idx_type>(Nz);
		const octave_idx_type nzmax = summa;
		const octave_idx_type rows = pituus;

		//NDArray sm(dim_vector(1, 1));

		// Create the MATLAB sparse matrix
		SparseMatrix sm(N, rows, nzmax);

		// Non-zero elements of the matrix
		double* elements = sm.data();

		//// Row indices
		size_t* indices = reinterpret_cast<size_t*>(sm.ridx());

		//// Column indices
		size_t* lor = reinterpret_cast<size_t*>(sm.cidx());

		for (size_t kk = 0; kk <= loop_var_par; kk++)
			lor[kk] = lor2[kk];

		NDArray Ll;
		// If doing Inveon attenuation
		if (attenuation_phase)
			Ll.resize(dim_vector(pituus, 1));
		else
			Ll.resize(dim_vector(1, 1));

		double* ll = Ll.fortran_vec();

		// Timing
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (projector_type == 2u) {

			//if (nrhs != 52)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 52.");

			// Width of the TOR
			const double crystal_size = prhs(ind).scalar_value();
			ind++;

			// Coordinates of the pixel centers in y-direction
			const NDArray x_center_ = prhs(ind).array_value();
			ind++;

			const double* x_center = x_center_.fortran_vec();

			// Coordinates of the pixel centers in x-direction
			const NDArray y_center_ = prhs(ind).array_value();
			ind++;

			const double* y_center = y_center_.fortran_vec();

			// Coordinates of the pixel centers in z-direction
			const NDArray z_center_ = prhs(ind).array_value();
			ind++;

			const double* z_center = z_center_.fortran_vec();

			const double crystal_size_z = prhs(ind).scalar_value();
			ind++;

			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			orth_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, scatter, scatter_coef);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				octave_stdout << "Orthogonal distance based ray tracer took " << (float)time_span.count() << " seconds\n";
			//	mexPrintf("Orthogonal distance based ray tracer took %f seconds\n", ((float)time_span.count()));
			//	mexEvalString("pause(.001);");
			}
		}
		else if (projector_type == 1u) {

			// run the Improved Siddon's algorithm, precomputed_lor = true
			improved_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, global_factor, scatter, scatter_coef);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				octave_stdout << "Improved Siddon took " << (float)time_span.count() << " seconds\n";
			//	mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
			//	mexEvalString("pause(.001);");
			}
		}
		else if ((projector_type == 3u)) {

			//if (nrhs != 52)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 52.");

			// Width of the TOR
			const double crystal_size = prhs(ind).scalar_value();
			ind++;

			// Coordinates of the pixel centers in y-direction
			const NDArray x_center_ = prhs(ind).array_value();
			ind++;

			const double* x_center = x_center_.fortran_vec();

			// Coordinates of the pixel centers in x-direction
			const NDArray y_center_ = prhs(ind).array_value();
			ind++;

			const double* y_center = y_center_.fortran_vec();

			// Coordinates of the pixel centers in z-direction
			const NDArray z_center_ = prhs(ind).array_value();
			ind++;

			const double* z_center = z_center_.fortran_vec();

			const double crystal_size_z = prhs(ind).scalar_value();
			ind++;

			const double bmin = prhs(ind).scalar_value();
			ind++;

			const double bmax = prhs(ind).scalar_value();
			ind++;

			const double Vmax = prhs(ind).scalar_value();
			ind++;

			const NDArray V_ = prhs(ind).array_value();
			ind++;

			const double* V = V_.fortran_vec();

			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			vol_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, bmin, 
				bmax, Vmax, V, scatter, scatter_coef);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				octave_stdout << "Volume-based ray tracer took " << (float)time_span.count() << " seconds\n";
			//	mexPrintf("Volume-based ray tracer took %f seconds\n", ((float)time_span.count()));
			//	mexEvalString("pause(.001);");
			}
		}
		//else
		//	mexErrMsgTxt("Unsupported projector");

		retval(0) = octave_value(sm);
		retval(1) = octave_value(Ll);

		return retval;

	}
	// Implementation 4
	else if (type == 1u) {

		//if (nrhs < 47)
		//	mexErrMsgTxt("Too few input arguments.  There must be at least 47.");
		//else if (nrhs > 56)
		//	mexErrMsgTxt("Too many input arguments.  There can be at most 56.");

		//if (nlhs != 2)
		//	mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Small constant to prevent division by zero
		const double epps = prhs(ind).scalar_value();
		ind++;

		// Measurement data
		const NDArray Sino_ = prhs(ind).array_value();
		ind++;

		const double* Sino = Sino_.fortran_vec();

		// Current estimates
		NDArray osem_apu_ = prhs(ind).array_value();
		ind++;

		double* osem_apu = osem_apu_.fortran_vec();

		// Projector used
		const uint32_t projector_type = prhs(ind).uint32_scalar_value();
		ind++;

		const size_t N = static_cast<size_t>(Nx) * static_cast<size_t>(Ny) * static_cast<size_t>(Nz);

		// If 1, do not compute the normalization constant anymore
		const bool no_norm = prhs(ind).bool_value();
		ind++;

		// precomputed_lor = false
		const bool precompute = prhs(ind).bool_value();
		ind++;

		const bool fp = prhs(ind).bool_value();
		ind++;

		const bool list_mode_format = prhs(ind).bool_value();
		ind++;

		NDArray Summ_(dim_vector(N, 1));

		// Normalization constants
		double* Summ = Summ_.fortran_vec();

		NDArray rhs_;
		if (fp)
			rhs_.resize(dim_vector(pituus, 1));
		else
			rhs_.resize(dim_vector(N, 1));

		double* rhs = rhs_.fortran_vec();

		//clock_t time = clock();
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		// Orthogonal
		if (projector_type == 2u) {

			//if (nrhs < 51)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 51.");

			// Width of the strip in 2D case
			const double crystal_size = prhs(ind).scalar_value();
			ind++;

			// Coordinates of the pixel centers in y-direction
			const NDArray x_center_ = prhs(ind).array_value();
			ind++;

			const double* x_center = x_center_.fortran_vec();

			// Coordinates of the pixel centers in x-direction
			const NDArray y_center_ = prhs(ind).array_value();
			ind++;

			const double* y_center = y_center_.fortran_vec();

			// Coordinates of the pixel centers in z-direction
			const NDArray z_center_ = prhs(ind).array_value();
			ind++;

			const double* z_center = z_center_.fortran_vec();

			// Width of the TOR in 3D case
			const double crystal_size_z = prhs(ind).scalar_value();
			ind++;

			if (precompute) {
				sequential_orth_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det, 
					NSlices, Nx, Ny, Nz, d, dz,	bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef);
			}
			else {
				sequential_orth_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef);
			}
		}
		// Improved Siddon
		else if (projector_type == 1u) {

			//if (nrhs < 48)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 48.");

			// Number of rays in Siddon (transaxial)
			uint16_t n_rays = prhs(ind).uint16_scalar_value();
			ind++;

			// Number of rays in Siddon (axial)
			uint16_t n_rays3D = prhs(ind).uint16_scalar_value();
			ind++;

			// Crystal pitch in z-direction (for multi-ray)
			const double cr_pz = prhs(ind).scalar_value();
			ind++;

			if (precompute) {
				sequential_improved_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, no_norm, global_factor, fp, scatter, scatter_coef, TOF, TOFSize,
					sigma_x, TOFCenter, nBins, dec_v);
			}
			else {
				sequential_improved_siddon_no_precompute(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, TotSinos,
					epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, cr_pz, no_norm, n_rays, n_rays3D, global_factor, fp, list_mode_format, 
					scatter, scatter_coef, TOF, TOFSize, sigma_x, TOFCenter, nBins, dec_v);
			}


		}
		else if ((projector_type == 3u)) {
			//if (nrhs < 54)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 54.");

			// Coordinates of the pixel centers in y-direction
			const NDArray x_center_ = prhs(ind).array_value();
			ind++;

			const double* x_center = x_center_.fortran_vec();

			// Coordinates of the pixel centers in x-direction
			const NDArray y_center_ = prhs(ind).array_value();
			ind++;

			const double* y_center = y_center_.fortran_vec();

			// Coordinates of the pixel centers in z-direction
			const NDArray z_center_ = prhs(ind).array_value();
			ind++;

			const double* z_center = z_center_.fortran_vec();

			// Width of the TOR in 3D case
			const double bmin = prhs(ind).scalar_value();
			ind++;

			// Width of the TOR in 3D case
			const double bmax = prhs(ind).scalar_value();
			ind++;

			// Width of the TOR in 3D case
			const double Vmax = prhs(ind).scalar_value();
			ind++;

			// Coordinates of the pixel centers in y-direction
			const NDArray V_ = prhs(ind).array_value();
			ind++;

			const double* V = V_.fortran_vec();

			if (precompute) {
				sequential_volume_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det,
					NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef);
			}
			else {
				sequential_volume_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef);
			}
		}

		retval(0) = octave_value(Summ_);
		retval(1) = octave_value(rhs_);

		return retval;
	}
	// Implementation 1, precomputed_lor = false
	else if (type == 2u) {

		//if (nrhs < 43)
		//	mexErrMsgTxt("Too few input arguments.  There must be at least 43.");
		//else if (nrhs > 49)
		//	mexErrMsgTxt("Too many input arguments.  There can be at most 49.");

		//if (nlhs != 3)
		//	mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// How many elements are preallocated in memory
		const uint32_t ind_size = prhs(ind).uint32_scalar_value();
		ind++;

		// Starting ring
		const uint32_t block1 = prhs(ind).uint32_scalar_value();
		ind++;

		// End ring
		const uint32_t blocks = prhs(ind).uint32_scalar_value();
		ind++;

		// Subset indices
		const uint32NDArray index_ = prhs(ind).uint32_array_value();
		const uint32_t index_size = index_.numel();
		ind++;

		const uint32_t* index = reinterpret_cast<const uint32_t*>(index_.fortran_vec());

		// Projector
		const uint32_t projector_type = prhs(ind).uint32_scalar_value();
		ind++;

		// Number of LORs
		if (index_size > 1ULL && !raw) {
			loop_var_par = index_size;
		}
		else if (!raw) {
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(size_x);
		}
		else {
			loop_var_par = numRows / 2ULL;
		}

		uint16NDArray lor_(dim_vector(loop_var_par, 1));

		uint16_t* lor = reinterpret_cast<uint16_t*>(lor_.fortran_vec());

		vector<uint32_t> indices;

		vector<double> elements;

		// Reserve some memory
		indices.reserve(ind_size);
		elements.reserve(ind_size);

		// The maximum elements of the pixel space in both x- and y-directions
		const double maxyy = yy_vec.back();
		const double maxxx = xx_vec.back();

		//clock_t time = clock();
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		uint32_t lj = 0U;

		if (projector_type == 1u) {

			// run the Improved Siddon's algorithm
			lj = improved_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
				det_per_ring, blocks, block1, L, pseudos, pRows, global_factor, scatter, scatter_coef);
		}
		else if (projector_type == 2u) {

			//if (nrhs != 49)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 49.");

			const double crystal_size = prhs(ind).scalar_value();
			ind++;

			const NDArray x_center_ = prhs(ind).array_value();
			ind++;

			const double* x_center = x_center_.fortran_vec();

			const NDArray y_center_ = prhs(ind).array_value();
			ind++;

			const double* y_center = y_center_.fortran_vec();

			// Coordinates of the pixel centers in z-direction
			const NDArray z_center_ = prhs(ind).array_value();
			ind++;

			const double* z_center = z_center_.fortran_vec();

			const double crystal_size_z = prhs(ind).scalar_value();
			ind++;

			// run the Orthogonal Siddon algorithm
			//lj = orth_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			//	yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
			//	det_per_ring, blocks, block1, L, pseudos, pRows, crystal_size, crystal_size_z, y_center, x_center, z_center, dec_v);
		}
		// Original Siddon's ray tracer
		else if (projector_type == 0u) {

			//if (nrhs < 46)
			//	mexErrMsgTxt("Too few input arguments.  There must be at least 46.");

			// Voxel numbers in x-direction
			const NDArray iij = prhs(ind).array_value();
			ind++;

			// Voxel numbers in y-direction
			const NDArray jji = prhs(ind).array_value();
			ind++;

			// Voxel numbers in z-direction
			const NDArray kkj = prhs(ind).array_value();
			ind++;

			vector<double> iij_vec(iij.numel(), 0.);
			std::copy_n(iij.data(), iij.numel(), iij_vec.begin());

			vector<double> jji_vec(jji.numel(), 0.);
			std::copy_n(jji.data(), jji.numel(), jji_vec.begin());

			vector<double> kkj_vec(kkj.numel(), 0.);
			std::copy_n(kkj.data(), kkj.numel(), kkj_vec.begin());

			// run the original Siddon's algorithm
			lj = original_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw,
				det_per_ring, blocks, block1, L, pseudos, pRows, iij_vec, jji_vec, kkj_vec, global_factor, scatter, scatter_coef);
		}

		const size_t outSize1 = static_cast<size_t>(lj) * 2ULL;
		const size_t outSize2 = 1ULL;

		// Create the Octave output vectors (row and column indices, elements)

		uint32NDArray outputMatrix2(dim_vector(indices.size(), outSize2));

		//for (size_t ll = 0ULL; ll < indices.size(); ll++)
		//	outputMatrix2(ll) = static_cast<octave_uint32>(indices[ll]);

		std::copy_n(indices.begin(), indices.size(), outputMatrix2.fortran_vec());
		indices.erase(indices.begin(), indices.end());
		indices.shrink_to_fit();

		NDArray outputMatrix3(dim_vector(elements.size(), outSize2));

		//for (size_t ll = 0ULL; ll < elements.size(); ll++)
		//	outputMatrix3(ll) = elements[ll];

		std::copy_n(elements.begin(), elements.size(), outputMatrix3.fortran_vec());

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		if (verbose) {
			octave_stdout << "Function elapsed time is " << (float)time_span.count() << " seconds\n";
		//	mexPrintf("Function elapsed time is %f seconds\n", ((float)time_span.count()));
		//	mexEvalString("pause(.001);");
		}

		retval(0) = octave_value(lor_);
		retval(1) = octave_value(outputMatrix2);
		retval(2) = octave_value(outputMatrix3);

		return retval;

	}
	// Precomputation phase
	else if (type == 3u) {
		//if (nrhs < 51)
		//	mexErrMsgTxt("Too few input arguments.  There must be at least 51.");
		//else if (nrhs > 51)
		//	mexErrMsgTxt("Too many input arguments.  There can be at most 51.");

		//if (nlhs != 3)
		//	mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// Starting ring
		const uint32_t block1 = prhs(ind).uint32_scalar_value();
		ind++;

		// End ring
		const uint32_t blocks = prhs(ind).uint32_scalar_value();
		ind++;

		// Projector
		const uint32_t projector_type = prhs(ind).uint32_scalar_value();
		ind++;

		const double crystal_size = prhs(ind).scalar_value();
		ind++;

		const NDArray x_center_ = prhs(ind).array_value();
		ind++;

		const double* x_center = x_center_.fortran_vec();

		const NDArray y_center_ = prhs(ind).array_value();
		ind++;

		const double* y_center = y_center_.fortran_vec();

		// Coordinates of the pixel centers in z-direction
		const NDArray z_center_ = prhs(ind).array_value();
		ind++;

		const double* z_center = z_center_.fortran_vec();

		const double crystal_size_z = prhs(ind).scalar_value();
		ind++;

		// Width of the TOR in 3D case
		const double bmin = prhs(ind).scalar_value();
		ind++;

		// Width of the TOR in 3D case
		const double bmax = prhs(ind).scalar_value();
		ind++;

		// Width of the TOR in 3D case
		const double Vmax = prhs(ind).scalar_value();
		ind++;

		// Coordinates of the pixel centers in y-direction
		const NDArray V_ = prhs(ind).array_value();
		ind++;

		const double* V = V_.fortran_vec();

		const uint32_t tyyppi = prhs(ind).uint32_scalar_value();
		ind++;

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(size_x);

		uint16NDArray lor_(dim_vector(loop_var_par, 1));
		uint16NDArray lor_orth_;
		uint16NDArray lor_vol_;

		if (tyyppi == 0) {
			lor_orth_.resize(dim_vector(1, 1));
			lor_vol_.resize(dim_vector(1, 1));
		}
		else if (tyyppi == 1 || tyyppi == 2) {

			if (tyyppi == 2u)
				lor_orth_.resize(dim_vector(loop_var_par * 2LL, 1));
			else
				lor_orth_.resize(dim_vector(loop_var_par, 1));
			lor_vol_.resize(dim_vector(1, 1));

		}
		else {
			lor_orth_.resize(dim_vector(1, 1));
			lor_vol_.resize(dim_vector(loop_var_par, 1));
		}

		uint16_t* lor = reinterpret_cast<uint16_t*>(lor_.fortran_vec());
		uint16_t* lor_orth = reinterpret_cast<uint16_t*>(lor_orth_.fortran_vec());
		uint16_t* lor_vol = reinterpret_cast<uint16_t*>(lor_vol_.fortran_vec());

		improved_siddon_precomputation_phase(loop_var_par, size_x, zmax, TotSinos, lor, maxyy, maxxx, xx_vec, z_det_vec, dy, yy_vec, x, y, z_det, NSlices, Nx, Ny, Nz,
			d, dz, bx, by, bz, block1, blocks, L, pseudos, raw, pRows, det_per_ring, tyyppi, lor_orth, lor_vol, crystal_size, crystal_size_z, x_center, y_center, z_center, 
			bmin, bmax, Vmax, V);


		retval(0) = octave_value(lor_);
		retval(1) = octave_value(lor_orth_);
		retval(2) = octave_value(lor_vol_);

		return retval;

	}
}