
/******************************************************************************
 * Matrix free projectors for forward and backward projections. This is the
 *branchless distance-driven projector. Based on:
 *https://doi.org/10.1109/TCI.2017.2675705
 *
 * Used by implementations 2, 3 and 5.
 *
 * The forward and backward projections are separate kernels.
 *
 *
 * INPUTS:
 * d_nRows = the number of detector elements (rows),
 * d_nCols = the number of detector elements (columns),
 * d_dPitch = Either a vector of float2 or two floats if PYTHON is defined, the
 *detector size/pitch in both "row" and "column" directions maskFP = 2D/3D
 *Forward projection mask, i.e. LORs/measurements with 0 will be skipped d_N =
 *image size in x/y/z- dimension, float3 or three floats (if PYTHON is defined),
 * d_b = distance from the pixel space to origin (z/x/y-dimension), float3 or
 *three floats (if PYTHON is defined), d_Size = precomputed scaling value,
 *float2 or two floats (if PYTHON is defined), see
 *computeProjectorScalingValues.m d_d = distance between adjecent voxels in
 *z/x/y-dimension, float3 or three floats (if PYTHON is defined), d_scale =
 *precomputed scaling value, float3 or three floats (if PYTHON is defined), see
 *computeProjectorScalingValues.m d_xy/z = detector x/y/z-coordinates, d_uv =
 *Direction coordinates for the detector panel, d_IImageY = Integral image of
 *the image volume for XZ-plane, d_IImageX = Integral image of the image volume
 *for YZ-plane, d_forw = forward projection, d_meanV = mean values for each
 *integral image slice, first for d_IImageY d_nProjections = Number of
 *projections/sinograms,
 *
 * OUTPUTS:
 * d_forw = forward projection,
 *
 * Copyright (C) 2022-2024 Ville-Veikko Wettenhovi
 *
 * This program is free software: you can redistribute it and/or modify it under
 *the terms of the GNU General Public License as published by the Free Software
 *Foundation, either version 3 of the License, or (at your option) any later
 *version.
 *
 * This program is distributed in the hope that it wiL be useful, but WITHOUT
 *ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *details.
 *
 * You should have received a copy of the GNU General Public License along with
 *this program. If not, see <https://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifdef OPENCL
CONSTANT sampler_t sampler2 =
    CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif

// Forward projection for the branchless distance-driven ray tracer
#ifdef FP
#ifdef OPENCL
__kernel
    __attribute__((vec_type_hint(float)))
    __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#else
KERNEL
#endif
void projectorType5Forward(
#if defined(METAL)
    SCALAR_PARAMS(scalarParams) BUF0
#else
    const uint d_nRows, const uint d_nCols
#ifdef PYTHON
    , const float d_dPitchX, const float d_dPitchY
#else
    , const float2 d_dPitch
#endif
#ifdef PYTHON
    , const uint d_Nx, const uint d_Ny, const uint d_Nz
    , const float bx, const float by, const float bz
    , const float d_SizeX, const float d_SizeY
    , const float d_dx, const float d_dy, const float d_dz
    , const float d_scalex, const float d_scaley, const float d_scalez
#else
    , const uint3 d_N
    , const float3 b
    , const float2 d_Size
    , const float3 d_d
    , const float3 d_scale
#endif
#endif
#if defined(USEGLOBAL)
    , const CLGLOBAL float *CLRESTRICT d_xyz BUF1
    , const CLGLOBAL float *CLRESTRICT d_uv BUF2
#else
    , CONSTANT float *d_xyz BUF1
    , CONSTANT float *d_uv BUF2
#endif
    , IMAGE3D d_IImageY TEX3
    , IMAGE3D d_IImageX TEX4
    , CLGLOBAL float *d_forw BUF5
#ifdef MEANDISTANCEFP
    , CONSTANT float *d_meanV BUF6
#endif
#ifdef MASKFP
#ifdef MASKFP3D
    , IMAGE3D maskFP TEX7
#else
    , IMAGE2D maskFP TEX7 
#endif
#endif
#ifdef NORM
    , const CLGLOBAL float *CLRESTRICT d_norm BUF8
#endif
#if !defined(METAL)
    , const LONG d_nProjections
#else
    , uint3 temp_i [[thread_position_in_grid]] // global id
#endif
) {
#if defined(METAL)
    UNPACK_SCALAR_PARAMS_5_FP(scalarParams);
    int GID0 = temp_i.x;
    int GID1 = temp_i.y;
    int GID2 = temp_i.z;
#endif
  const int3 i = MINT3(GID0, GID1 * NVOXELSFP, GID2);

  if (i.x >= d_nRows || i.y >= d_nCols || i.z >= d_nProjections)
    return;
  size_t idx = GID0 + GID1 * NVOXELSFP * d_nRows + GID2 * d_nCols * d_nRows;

#ifdef MASKFP
#ifdef CUDA
#ifdef MASKFP3D
  const int maskVal = tex3D<unsigned char>(maskFP, i.x, i.y, i.z);
#else
  const int maskVal = tex2D<unsigned char>(maskFP, i.x, i.y);
#endif
#else
#ifdef MASKFP3D
  const int maskVal =
      read_imageui(maskFP, sampler_MASK, (int4)(i.x, i.y, i.z, 0)).w;
#else
  const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
#endif
  if (maskVal == 0)
    return;
#endif
#ifdef PYTHON
  const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
  const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
  const float3 d_d = make_float3(d_dx, d_dy, d_dz);
  const float3 b = make_float3(bx, by, bz);
  const float3 d_scale = make_float3(d_scalex, d_scaley, d_scalez);
  const float2 d_Size = make_float2(d_SizeX, d_SizeY);
#endif
  float3 s, d, dL, dR, dU, dD;
  float temp[NVOXELSFP];
  for (int zz = 0; zz < NVOXELSFP; zz++)
    temp[zz] = 0.f;
#ifdef CT
  // Get source and detector coordinates
  getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_nRows, d_nCols, d_dPitch,
                           &dL, &dR, &dU, &dD);
  // Weight values
  const float3 uVector = FABS(NORMALIZE(d - s));
  float kerroin = d_d.x * d_d.y * d_d.z;
  const float apuZ = (b.z - s.z - d_d.z / 2.f);
  if (FABS(s.x) <= FABS(s.y)) {
    // Coordinates of the corners
    const float2 X =
        MFLOAT2((dL.x - s.x) / (dL.y - s.y), (dR.x - s.x) / (dR.y - s.y));
    const float2 Z =
        MFLOAT2((dU.z - s.z) / (dU.y - s.y), (dD.z - s.z) / (dD.y - s.y));
    kerroin /= uVector.y;
    const float apuX = (b.x - s.x - d_d.x / 2.f);
    const float apuV = b.y - s.y;
    for (uint jj = 0; jj < d_N.y; jj++) {
      float dy = CFLOAT(jj) * d_d.y + d_d.y / 2.f;
      const float v = dy + apuV;
      float2 xLR = v * X;
      float2 zUD = v * Z;
      // Shift the origin to the corner of the slice
      xLR -= apuX;
      zUD -= apuZ;
      if (xLR.x > xLR.y) {
        const float apu = xLR.x;
        xLR.x = xLR.y;
        xLR.y = apu;
      }
      if (zUD.x > zUD.y) {
        const float apu = zUD.x;
        zUD.x = zUD.y;
        zUD.y = apu;
      }
      const float area = FABS((xLR.x - xLR.y) * (zUD.x - zUD.y));
      // Scale between 0 and 1
      xLR *= d_scale.x;
      dy *= d_Size.y;
      zUD *= d_scale.z;
#ifdef CUDA
      float A = tex3D<float>(d_IImageY, xLR.y, zUD.x, dy);
      float B = tex3D<float>(d_IImageY, xLR.x, zUD.y, dy);
      float C = tex3D<float>(d_IImageY, xLR.y, zUD.y, dy);
      float D = tex3D<float>(d_IImageY, xLR.x, zUD.x, dy);
#elif defined(OPENCL)
      float A =
          read_imagef(d_IImageY, sampler2, (float4)(xLR.y, zUD.x, dy, 0.f)).w;
      float B =
          read_imagef(d_IImageY, sampler2, (float4)(xLR.x, zUD.y, dy, 0.f)).w;
      float C =
          read_imagef(d_IImageY, sampler2, (float4)(xLR.y, zUD.y, dy, 0.f)).w;
      float D =
          read_imagef(d_IImageY, sampler2, (float4)(xLR.x, zUD.x, dy, 0.f)).w;
#elif defined(METAL)
      float A = d_IImageY.sample(sampler2, (float3)(xLR.y, zUD.x, dy)).r;
      float B = d_IImageY.sample(sampler2, (float3)(xLR.x, zUD.y, dy)).r;
      float C = d_IImageY.sample(sampler2, (float3)(xLR.y, zUD.y, dy)).r;
      float D = d_IImageY.sample(sampler2, (float3)(xLR.x, zUD.x, dy)).r;
#endif
      float apu = C + D - A - B;
#ifdef MEANDISTANCEFP
#ifdef CUDA
      const float area2 =
          d_meanV[jj + d_N.x] * fabs((xLR.x - xLR.y) * (zUD.x - zUD.y)) *
          CFLOAT(get_image_width(d_IImageY) * get_image_height(d_IImageY));
#else
      const float area2 = d_meanV[jj + d_N.x] *
                          FABS((xLR.x - xLR.y) * (zUD.x - zUD.y)) *
                          CFLOAT((d_N.y + 1) * (d_N.z + 1));
#endif
      apu += area2;
#endif
      temp[0] += apu / area;
      const float xInterval = FABS(zUD.x - zUD.y);
      for (int zz = 1; zz < NVOXELSFP; zz++) {
        const uint ind = i.y + zz;
        if (ind >= d_nCols)
          break;
        D = B;
        A = C;
        zUD.y += xInterval;
#if defined(CUDA)
        B = tex3D<float>(d_IImageY, xLR.x, zUD.y, dy);
        C = tex3D<float>(d_IImageY, xLR.y, zUD.y, dy);
#elif defined(OPENCL)
        B = read_imagef(d_IImageY, sampler2, (float4)(xLR.x, zUD.y, dy, 0.f)).w;
        C = read_imagef(d_IImageY, sampler2, (float4)(xLR.y, zUD.y, dy, 0.f)).w;
#elif defined(METAL)
        B = d_IImageY.sample(sampler2, (float3)(xLR.x, zUD.y, dy)).r;
        C = d_IImageY.sample(sampler2, (float3)(xLR.y, zUD.y, dy)).r;
#endif
        apu = C + D - A - B;
#ifdef MEANDISTANCEFP
        apu += (area2);
#endif
        temp[zz] += apu / area;
      }
    }
  } else {
    const float2 Y =
        MFLOAT2((dL.y - s.y) / (dL.x - s.x), (dR.y - s.y) / (dR.x - s.x));
    const float2 Z =
        MFLOAT2((dU.z - s.z) / (dU.x - s.x), (dD.z - s.z) / (dD.x - s.x));
    kerroin /= uVector.x;
    const float apuY = (b.y - s.y - d_d.y / 2.f);
    const float apuV = b.x - s.x;
    for (uint ii = 0; ii < d_N.x; ii++) {
      float dx = CFLOAT(ii) * d_d.x + d_d.x / 2.f;
      const float v = dx + apuV;
      float2 yLR = v * Y;
      float2 zUD = v * Z;
      yLR -= apuY;
      zUD -= apuZ;
      if (yLR.x > yLR.y) {
        const float apu = yLR.x;
        yLR.x = yLR.y;
        yLR.y = apu;
      }
      if (zUD.x > zUD.y) {
        const float apu = zUD.x;
        zUD.x = zUD.y;
        zUD.y = apu;
      }
      const float area = FABS((yLR.x - yLR.y) * (zUD.x - zUD.y));
      yLR *= d_scale.y;
      dx *= d_Size.x;
      zUD *= d_scale.z;
#if defined(CUDA)
      float A = tex3D<float>(d_IImageX, yLR.y, zUD.x, dx);
      float B = tex3D<float>(d_IImageX, yLR.x, zUD.y, dx);
      float C = tex3D<float>(d_IImageX, yLR.y, zUD.y, dx);
      float D = tex3D<float>(d_IImageX, yLR.x, zUD.x, dx);
#elif defined(OPENCL)
      float A =
          read_imagef(d_IImageX, sampler2, (float4)(yLR.y, zUD.x, dx, 0.f)).w;
      float B =
          read_imagef(d_IImageX, sampler2, (float4)(yLR.x, zUD.y, dx, 0.f)).w;
      float C =
          read_imagef(d_IImageX, sampler2, (float4)(yLR.y, zUD.y, dx, 0.f)).w;
      float D =
          read_imagef(d_IImageX, sampler2, (float4)(yLR.x, zUD.x, dx, 0.f)).w;
#elif defined(METAL)
      float A = d_IImageX.sample(sampler2, (float3)(yLR.y, zUD.x, dx)).r;
      float B = d_IImageX.sample(sampler2, (float3)(yLR.x, zUD.y, dx)).r;
      float C = d_IImageX.sample(sampler2, (float3)(yLR.y, zUD.y, dx)).r;
      float D = d_IImageX.sample(sampler2, (float3)(yLR.x, zUD.x, dx)).r;
#endif
      float apu = C + D - A - B;
#ifdef MEANDISTANCEFP
#ifdef CUDA
      const float area2 =
          d_meanV[ii] * fabs((yLR.x - yLR.y) * (zUD.x - zUD.y)) *
          CFLOAT(get_image_width(d_IImageX) * get_image_height(d_IImageX));
#else
      const float area2 = d_meanV[ii] *
                          FABS((yLR.x - yLR.y) * (zUD.x - zUD.y)) *
                          CFLOAT((d_N.x + 1) * (d_N.z + 1));
#endif
      apu += (area2);
#endif
      temp[0] += apu / area;
      const float yInterval = FABS(zUD.x - zUD.y);
      for (int zz = 1; zz < NVOXELSFP; zz++) {
        const uint ind = i.y + zz;
        if (ind >= d_nCols)
          break;
        D = B;
        A = C;
        zUD.y += yInterval;
#if defined(CUDA)
        B = tex3D<float>(d_IImageX, yLR.x, zUD.y, dx);
        C = tex3D<float>(d_IImageX, yLR.y, zUD.y, dx);
#elif defined(OPENCL)
        B = read_imagef(d_IImageX, sampler2, (float4)(yLR.x, zUD.y, dx, 0.f)).w;
        C = read_imagef(d_IImageX, sampler2, (float4)(yLR.y, zUD.y, dx, 0.f)).w;
#elif defined(METAL)
        B = d_IImageX.sample(sampler2, (float3)(yLR.x, zUD.y, dx)).r;
        C = d_IImageX.sample(sampler2, (float3)(yLR.y, zUD.y, dx)).r;
#endif
        apu = C + D - A - B;
#ifdef MEANDISTANCEFP
        apu += (area2);
#endif
        temp[zz] += apu / area;
      }
    }
  }
  for (int zz = 0; zz < NVOXELSFP; zz++) {
    if (idx >= d_nCols * d_nRows * d_nProjections)
      break;
#if defined(NORM)
    temp[zz] *= local_norm;
#endif
    d_forw[idx] += temp[zz] * kerroin;
    idx += d_nRows;
  }
#else

#endif
}
#endif

#ifdef BP
#ifndef NVOXELS5
#define NVOXELS5 1
#endif

/******************************************************************************
 * BDD backprojection.
 *
 * INPUTS:
 * d_nRows = the number of detector elements (rows),
 * d_nCols = the number of detector elements (columns),
 * d_dPitch = Either a vector of float2 or two floats if PYTHON is defined, the
 *detector size/pitch in both "row" and "column" directions maskBP = 2D/3D
 *backward projection mask, i.e. voxels with 0 will be skipped T = redundancy
 *weights for offset imaging, d_N = image size in x/y/z- dimension, float3 or
 *three floats (if PYTHON is defined), d_b = distance from the pixel space to
 *origin (z/x/y-dimension), float3 or three floats (if PYTHON is defined), d_d =
 *distance between adjecent voxels in z/x/y-dimension, float3 or three floats
 *(if PYTHON is defined), d_scale = precomputed scaling value, float3 or three
 *floats (if PYTHON is defined), see computeProjectorScalingValues.m d_Size =
 *precomputed scaling value, float2 or two floats (if PYTHON is defined), see
 *computeProjectorScalingValues.m d_xy/z = detector x/y/z-coordinates, d_uv =
 *Direction coordinates for the detector panel, d_IImage = Integral image of the
 *forward projections, d_forw = backward projection, d_Summ = buffer for d_Summ
 *(sensitivity image), d_meanV = mean values for each integral image slice,
 * d_nProjections = Number of projections/sinograms,
 * ii = The current volume, for multi-resolution reconstruction, 0 can be used
 *when not using multi-resolution
 *
 * OUTPUTS:
 * d_forw = backprojection,
 ******************************************************************************/

#ifdef OPENCL
__kernel
    __attribute__((vec_type_hint(float)))
    __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#else
extern "C" __global__
#endif
    void
    projectorType5Backward(const uint d_nRows, const uint d_nCols,
#ifdef PYTHON
                           const float d_dPitchX, const float d_dPitchY,
#else
                           const float2 d_dPitch,
#endif
#ifdef OFFSET
                           CONSTANT float *T,
#endif
#ifdef PYTHON
                           const uint d_Nx, const uint d_Ny, const uint d_Nz,
                           const float bx, const float by, const float bz,
                           const float d_dx, const float d_dy, const float d_dz,
                           const float d_scalex, const float d_scaley,
                           const float d_scalez, const float d_SizeX,
                           const float d_SizeY,
#else
                           const uint3 d_N, const float3 b, const float3 d_d,
                           const float3 d_scale, const float2 d_Size,
#endif
#if defined(USEGLOBAL)
                           const CLGLOBAL float *CLRESTRICT d_xyz,
#else
                           CONSTANT float *d_xyz,
#endif
                           CONSTANT float *d_uv, IMAGE3D d_IImage,
                           CLGLOBAL float *d_forw, CLGLOBAL float *d_Summ,
#ifdef MEANDISTANCEBP
                           CONSTANT float *d_meanV,
#endif
#ifdef NORM
                           const CLGLOBAL float *CLRESTRICT d_norm,
#endif
                           const uchar no_norm,
#ifdef MASKBP
#ifdef MASKBP3D
                           IMAGE3D maskBP,
#else
                           IMAGE2D maskBP,
#endif
#endif
                           const LONG d_nProjections, const int ii) {
  const int3 i = MINT3(GID0, GID1, GID2 * NVOXELS5);
#ifdef PYTHON
  const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
#endif
  if (i.x >= d_N.x || i.y >= d_N.y || i.z >= d_N.z)
    return;
  size_t idx = GID0 + GID1 * d_N.x + GID2 * NVOXELS5 * d_N.y * d_N.x;
#ifdef MASKBP
  if (ii == 0) {
#ifdef CUDA
#ifdef MASKBP3D
    const int maskVal = tex3D<unsigned char>(maskBP, i.x, i.y, i.z);
#else
    const int maskVal = tex2D<unsigned char>(maskBP, i.x, i.y);
#endif
#else
#ifdef MASKBP3D
    const int maskVal =
        read_imageui(maskBP, sampler_MASK, (int4)(i.x, i.y, i.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
#endif
    if (maskVal == 0)
      return;
  }
#endif
#ifdef PYTHON
  const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
  const float3 d_d = make_float3(d_dx, d_dy, d_dz);
  const float3 b = make_float3(bx, by, bz);
  const float3 d_scale = make_float3(d_scalex, d_scaley, d_scalez);
  const float2 d_Size = make_float2(d_SizeX, d_SizeY);
#endif
  float temp[NVOXELS5];
  float wSum[NVOXELS5];
  for (int zz = 0; zz < NVOXELS5; zz++) {
    temp[zz] = 0.f;
    if (no_norm == 0u)
      wSum[zz] = 0.f;
  }
#ifdef USEMAD
  const float3 dV = FMAD3(CFLOAT3(i), d_d, d_d / 2.f + b);
#else
  const float3 dV = CFLOAT3(i) * d_d + d_d / 2.f + b;
#endif
  const float2 koko = MFLOAT2(CFLOAT(d_nRows + 1) * d_dPitch.x,
                              CFLOAT(d_nCols + 1) * d_dPitch.y);
  const float2 indeksi = MFLOAT2(CFLOAT(d_nRows) / 2.f, CFLOAT(d_nCols) / 2.f);
#ifdef CT
  float3 dV2 = dV;
  for (int kk = 0; kk < d_nProjections; kk++) {
    float3 vLU, vRD;
    const int id = kk * 6;
    const float3 s = CMFLOAT3(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
    const float3 d = CMFLOAT3(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
#if defined(PITCH)
    const float3 apuX =
        CMFLOAT3(d_uv[kk * NA], d_uv[kk * NA + 1], d_uv[kk * NA + 2]) *
        indeksi.x;
    const float3 apuY =
        CMFLOAT3(d_uv[kk * NA + 3], d_uv[kk * NA + 4], d_uv[kk * NA + 5]) *
        indeksi.y;
#else
    const float3 apuX =
        CMFLOAT3(indeksi.x * d_uv[kk * NA], indeksi.x * d_uv[kk * NA + 1], 0.f);
    const float3 apuY = CMFLOAT3(0.f, 0.f, indeksi.y * d_dPitch.y);
#endif
    const float3 d2 = apuX - apuY;
    const float3 d3 = d - apuX - apuY;
    const float3 normX = normalize(apuX);
    const float3 normY = normalize(apuY);
    const float3 v = normalize(dV - s);
    const float3 aVar = d3 - d;
    const float3 crossP = cross(d2, d3 - d);
    const float upperPart = dot(crossP, s - d);
    if (fabs(s.x) <= fabs(s.y)) {
      if (s.y > 0.f) {
        vLU = CMFLOAT3(dV.x - d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s; // A
        vRD = CMFLOAT3(dV.x + d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s; // D
      } else {
        vLU = CMFLOAT3(dV.x + d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s;
        vRD = CMFLOAT3(dV.x - d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s;
      }
    } else {
      if (s.x > 0.f) {
        vLU = CMFLOAT3(dV.x, dV.y + d_d.y / 2.f, dV.z - d_d.z / 2.f) - s; // A
        vRD = CMFLOAT3(dV.x, dV.y - d_d.y / 2.f, dV.z + d_d.z / 2.f) - s; // D
      } else {
        vLU = CMFLOAT3(dV.x, dV.y - d_d.y / 2.f, dV.z - d_d.z / 2.f) - s;
        vRD = CMFLOAT3(dV.x, dV.y + d_d.y / 2.f, dV.z + d_d.z / 2.f) - s;
      }
    }
    float tLU = upperPart / (dot(-vLU, crossP));
    ;
    float tRD = upperPart / (dot(-vRD, crossP));
#ifdef USEMAD
    float3 pLU = FMAD3(vLU, tLU, s);
    float3 pRD = FMAD3(vRD, tRD, s);
#else
    float3 pLU = vLU * tLU + s;
    float3 pRD = vRD * tRD + s;
#endif
    pLU -= d3;
    pRD -= d3;
    float Ax = dot(pLU, normX) + d_dPitch.x / 2.f;
    float Ay = dot(pLU, normY) + d_dPitch.y / 2.f;
    float Dx = dot(pRD, normX) + d_dPitch.x / 2.f;
    float Dy = dot(pRD, normY) + d_dPitch.y / 2.f;
    if (Dx > Ax) {
      const float xA = Ax;
      Ax = Dx;
      Dx = xA;
    }
    if (Dy > Ay) {
      const float yA = Ay;
      Ay = Dy;
      Dy = yA;
    }
    const float dz = (CFLOAT(kk) + 0.5f) / CFLOAT(d_nProjections);
    float3 coordA = MFLOAT3(Ax / koko.x, Ay / koko.y, dz);
    float3 coordB = MFLOAT3(Ax / koko.x, Dy / koko.y, dz);
    float3 coordC = MFLOAT3(Dx / koko.x, Ay / koko.y, dz);
    float3 coordD = MFLOAT3(Dx / koko.x, Dy / koko.y, dz);
#ifdef CUDA
    float A = tex3D<float>(d_IImage, coordA.x, coordA.y, coordA.z);
    float B = tex3D<float>(d_IImage, coordB.x, coordB.y, coordB.z);
    float C = tex3D<float>(d_IImage, coordC.x, coordC.y, coordC.z);
    float D = tex3D<float>(d_IImage, coordD.x, coordD.y, coordD.z);
#else
    float A = read_imagef(d_IImage, sampler2, (float4)(coordA, 0.f)).w;
    float B = read_imagef(d_IImage, sampler2, (float4)(coordB, 0.f)).w;
    float C = read_imagef(d_IImage, sampler2, (float4)(coordC, 0.f)).w;
    float D = read_imagef(d_IImage, sampler2, (float4)(coordD, 0.f)).w;
#endif
#ifdef NORM
    LONG indX = CLONG_rtz(coordA.x * CFLOAT(d_size_x));
    LONG indY = CLONG_rtz(coordA.y * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
    LONG indZ = CLONG_rtz(coordA.z * CFLOAT(d_nProjections)) *
                CLONG_rtz(d_sizey) * CLONG_rtz(d_size_x);
    A *= d_norm[indX + indY + indZ];
    indX = CLONG_rtz(coordB.x * CFLOAT(d_size_x));
    indY = CLONG_rtz(coordB.y * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
    indZ = CLONG_rtz(coordB.z * CFLOAT(d_nProjections)) * CLONG_rtz(d_sizey) *
           CLONG_rtz(d_size_x);
    B *= d_norm[indX + indY + indZ];
    indX = CLONG_rtz(coordB.x * CFLOAT(d_size_x));
    indY = CLONG_rtz(coordB.y * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
    indZ = CLONG_rtz(coordB.z * CFLOAT(d_nProjections)) * CLONG_rtz(d_sizey) *
           CLONG_rtz(d_size_x);
    C *= d_norm[indX + indY + indZ];
    indX = CLONG_rtz(coordB.x * CFLOAT(d_size_x));
    indY = CLONG_rtz(coordB.y * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
    indZ = CLONG_rtz(coordB.z * CFLOAT(d_nProjections)) * CLONG_rtz(d_sizey) *
           CLONG_rtz(d_size_x);
    D *= d_norm[indX + indY + indZ];
#endif
    float kerroin;
    if (fabs(s.x) <= fabs(s.y))
      kerroin = d_d.y / fabs(v.y);
    else
      kerroin = d_d.x / fabs(v.x);
#ifdef OFFSET
    float OffTT;
    float OffT = T[kk];
    if (OffT > koko.x / 2.f) {
      coordA.x = fabs(coordA.x - 1.f);
      coordB.x = fabs(coordB.x - 1.f);
      coordC.x = fabs(coordC.x - 1.f);
      coordD.x = fabs(coordD.x - 1.f);
      OffTT = koko.x - OffT;
    } else
      OffTT = OffT;
    coordA.x *= koko.x;
    coordB.x *= koko.x;
    coordC.x *= koko.x;
    coordD.x *= koko.x;
    coordA.x -= OffTT;
    coordB.x -= OffTT;
    coordC.x -= OffTT;
    coordD.x -= OffTT;
    float wA = 1.f, wB = 1.f, wC = 1.f, wD = 1.f;
    if (coordA.x <= OffTT && coordA.x >= -OffTT) {
      wA = SINF(M_PI_2_F * ((coordA.x + OffTT) / (2.f * OffTT)));
      wA *= wA;
    } else if (coordA.x < -OffTT)
      wA = 0.f;
    if (coordB.x <= OffTT && coordB.x >= -OffTT) {
      wB = SINF(M_PI_2_F * ((coordB.x + OffTT) / (2.f * OffTT)));
      wB *= wB;
    } else if (coordB.x < -OffTT)
      wB = 0.f;
    if (coordC.x <= OffTT && coordC.x >= -OffTT) {
      wC = SINF(M_PI_2_F * ((coordC.x + OffTT) / (2.f * OffTT)));
      wC *= wC;
    } else if (coordC.x < -OffTT)
      wC = 0.f;
    if (coordD.x <= OffTT && coordD.x >= -OffTT) {
      wD = SINF(M_PI_2_F * ((coordD.x + OffTT) / (2.f * OffTT)));
      wD *= wD;
    } else if (coordD.x < -OffTT)
      wD = 0.f;

    const float w = (wA + wB + wC + wD) / 4.f;
    float apu = (A + D - C - B) * w;
#else
    float apu = (A + D - C - B);
#endif
#ifdef MEANDISTANCEBP
    float area = d_meanV[kk] * CFLOAT((d_nRows + 1) * (d_nCols + 1));
#endif
    if (apu != 0.f) {
#ifdef MEANDISTANCEBP
      apu += area * fabs(coordA.y - coordD.y) * fabs(coordA.x - coordD.x);
#endif
#ifdef FDK
      kerroin /= 6.f;
#endif
      temp[0] += apu * kerroin;
      if (no_norm == 0u)
        wSum[0] += kerroin;
    }
    for (int zz = 1; zz < NVOXELS5; zz++) {
      const uint ind = i.z + zz;
      if (ind >= d_N.z)
        break;
      vRD.z += d_d.z;
      vLU.z += d_d.z;
      dV2.z = CFLOAT(i.z + zz) * d_d.z + d_d.z / 2.f + b.z;
      float3 v2 = normalize(dV2 - s);
      if (fabs(s.x) <= fabs(s.y))
        kerroin = d_d.y / fabs(v2.y);
      else
        kerroin = d_d.x / fabs(v2.x);
      tLU = upperPart / (dot(-vLU, crossP));
      tRD = upperPart / (dot(-vRD, crossP));
#ifdef USEMAD
      pLU = FMAD3(vLU, tLU, s);
      pRD = FMAD3(vRD, tRD, s);
#else
      pLU = vLU * tLU + s;
      pRD = vRD * tRD + s;
#endif
      pRD -= d3;
      pLU -= d3;
      Ay = dot(pLU, normY) + d_dPitch.y / 2.f;
      Dy = dot(pRD, normY) + d_dPitch.y / 2.f;
      if (Dy > Ay) {
        const float yA = Ay;
        Ay = Dy;
        Dy = yA;
        B = A;
        D = C;
        coordA = CMFLOAT3(Ax / koko.x, Ay / koko.y, dz);
        coordC = CMFLOAT3(Dx / koko.x, Ay / koko.y, dz);
#ifdef CUDA
        A = tex3D<float>(d_IImage, coordA.x, coordA.y, coordA.z);
        C = tex3D<float>(d_IImage, coordC.x, coordC.y, coordC.z);
#else
        A = read_imagef(d_IImage, sampler2, (float4)(coordA, 0.f)).w;
        C = read_imagef(d_IImage, sampler2, (float4)(coordC, 0.f)).w;
#endif
#ifdef MEANDISTANCEBP
        coordD = (float4)(Dx / koko.x, Dy / koko.y, dz, 0.f);
#endif
      } else {
        A = B;
        C = D;
        coordB = CMFLOAT3(Ax / koko.x, Dy / koko.y, dz);
        coordD = CMFLOAT3(Dx / koko.x, Dy / koko.y, dz);
#ifdef CUDA
        B = tex3D<float>(d_IImage, coordB.x, coordB.y, coordB.z);
        D = tex3D<float>(d_IImage, coordD.x, coordD.y, coordD.z);
#else
        B = read_imagef(d_IImage, sampler2, (float4)(coordB, 0.f)).w;
        D = read_imagef(d_IImage, sampler2, (float4)(coordD, 0.f)).w;
#endif
#ifdef MEANDISTANCEBP
        coordA = (float4)(Ax / koko.x, Ay / koko.y, dz, 0.f);
#endif
      }
#ifdef OFFSET
      apu = (A * wA + D * wD - C * wC - B * wB);
#else
      apu = (A + D - C - B);
#endif
      if (apu != 0.f) {
#ifdef MEANDISTANCEBP
        apu += area * fabs(coordA.y - coordD.y) * fabs(coordA.x - coordD.x);
#endif
#ifdef FDK
        kerroin /= 6.f;
#endif
        temp[zz] += apu * kerroin;
        if (no_norm == 0u)
          wSum[zz] += kerroin;
      }
    }
  }
  for (int zz = 0; zz < NVOXELS5; zz++) {
    const uint ind = i.z + zz;
    if (ind >= d_N.z)
      break;
    d_forw[idx] = temp[zz];
    if (no_norm == 0u)
      d_Summ[idx] = wSum[zz];
    idx += d_N.y * d_N.x;
  }
#endif
}
#endif
