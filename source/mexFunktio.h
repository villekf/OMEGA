#pragma once
#include <cstdio>
#include <cstdint>
#ifndef OCTAVE
#include "mex.h"
//#undef MX_HAS_INTERLEAVED_COMPLEX

const bool getScalarBool(const mxArray* mx, int ind);

const int8_t getScalarInt8(const mxArray* mx, int ind);

const uint8_t getScalarUInt8(const mxArray* mx, int ind);

const int16_t getScalarInt16(const mxArray* mx, int ind);

const uint16_t getScalarUInt16(const mxArray* mx, int ind);

const int32_t getScalarInt32(const mxArray* mx, int ind);

const uint32_t getScalarUInt32(const mxArray* mx, int ind);

const int64_t getScalarInt64(const mxArray* mx, int ind);

const uint64_t getScalarUInt64(const mxArray* mx, int ind);

const float getScalarFloat(const mxArray* mx, int ind);

const double getScalarDouble(const mxArray* mx, int ind);
#endif