__kernel void Convolution3D(__global float* input,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z, 
	const uint Nx, const uint Ny, const uint Nz, const uint Nyx) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (0, 0, 0, 0);
	float result = 0.f;
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		for (int j = -window_size_y; j <= window_size_y; j++) {
			for (int i = -window_size_x; i <= window_size_x; i++) {
				ind_uus.x = ind.x + i;
				ind_uus.y = ind.y + j;
				ind_uus.z = ind.z + k;
				if (ind_uus.x >= Nx)
					ind_uus.x = ind.x - i + 1;
				else if (ind_uus.x < 0)
					ind_uus.x = ind.x - (i + 1);
				if (ind_uus.y >= Ny)
					ind_uus.y = ind.y - j + 1;
				else if (ind_uus.y < 0)
					ind_uus.y = ind.y - (j + 1);
				if (ind_uus.z >= Nz)
					ind_uus.z = ind.z - k + 1;
				else if (ind_uus.z < 0)
					ind_uus.z = ind.z - (k + 1);
				int indeksi = ind_uus.x + ind_uus.y * Nx + ind_uus.z * Nyx;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[ind.x + ind.y * Nx + ind.z * Nyx] = result;
}