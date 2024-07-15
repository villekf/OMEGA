#pragma once
//#define DEBUG false
#include "mex.h"
#include <string>
#include <cstring>
#ifndef TARGET_API_VERSION
#define TARGET_API_VERSION 0
#endif
//#undef MX_HAS_INTERLEAVED_COMPLEX

inline mxArray* getField(const mxArray* in, mwIndex i, const char* name) {
    const int field_num = mxGetFieldNumber(in, name);
    std::string output = name;
    output = "Error with " + output;
    mxArray* out = nullptr;
    if (field_num >= 0)
        out = mxGetField(in, i, name);
    else
        mexErrMsgTxt(output.c_str());
    return out;
}

inline const bool getScalarBool(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const bool scalar = (bool)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar bool is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const bool scalar = (bool)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar bool is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const int8_t getScalarInt8(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const int8_t scalar = (int8_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int8 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const int8_t scalar = (int8_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int8 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const uint8_t getScalarUInt8(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const uint8_t scalar = (uint8_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint8 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const uint8_t scalar = (uint8_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint8 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const int16_t getScalarInt16(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const int16_t scalar = (int16_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int16 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const int16_t scalar = (int16_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int16 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const uint16_t getScalarUInt16(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const uint16_t scalar = (uint16_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint16 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const uint16_t scalar = (uint16_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint16 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const int32_t getScalarInt32(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const int32_t scalar = (int32_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int32 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const int32_t scalar = (int32_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int32 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const uint32_t getScalarUInt32(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const uint32_t scalar = (uint32_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint32 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const uint32_t scalar = (uint32_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint32 is not of size == [1 1].\n");
            return false;
        }
    }
}


inline const int64_t getScalarInt64(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const int64_t scalar = (int64_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int64 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const int64_t scalar = (int64_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar int64 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const uint64_t getScalarUInt64(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const uint64_t scalar = (uint64_t)mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint64 is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const uint64_t scalar = (uint64_t)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar uint64 is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const float getScalarFloat(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const float scalar = (float)mxGetScalar(mx);
            return scalar;
        }
        else if (mxGetNumberOfElements(mx) == 0)
            return 0.f;
        else {
            mexErrMsgTxt("Scalar float is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const float scalar = (float)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar float is not of size == [1 1].\n");
            return false;
        }
    }
}

inline const double getScalarDouble(const mxArray* mx, int ind = 0, const char* name = "") {
    if (std::strcmp(name, "") == 0) {
        if (mxGetNumberOfElements(mx) == 1) {
            const double scalar = mxGetScalar(mx);
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar double is not of size == [1 1].\n");
            return false;
        }
    }
    else {
        if (mxGetNumberOfElements(getField(mx, ind, name)) == 1) {
            const double scalar = (double)mxGetScalar(getField(mx, 0, name));
            return scalar;
        }
        else {
            mexErrMsgTxt("Scalar double is not of size == [1 1].\n");
            return false;
        }
    }
}

inline float* getSingles(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (float*)mxGetSingles(mxGetCell(mx, i));
#else
        return (float*)mxGetData(mxGetCell(mx, i));
#endif
    }
    else if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (float*)mxGetSingles(mx);
#else
        return (float*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (float*)mxGetSingles(getField(mx, i, name));
#else
        return (float*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline double* getDoubles(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (double*)mxGetDoubles(mxGetCell(mx, i));
#else
        return (double*)mxGetData(mxGetCell(mx, i));
#endif
    }
    else if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (double*)mxGetDoubles(mx);
#else
        return (double*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (double*)mxGetDoubles(getField(mx, i, name));
#else
        return (double*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline uint32_t* getUint32s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint32_t*)mxGetUint32s(mx);
#else
        return (uint32_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint32_t*)mxGetUint32s(getField(mx, i, name));
#else
        return (uint32_t*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline int32_t* getInt32s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (int32_t*)mxGetInt32s(mx);
#else
        return (int32_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (int32_t*)mxGetInt32s(getField(mx, i, name));
#else
        return (int32_t*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline uint64_t* getUint64s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint64_t*)mxGetUint64s(mx);
#else
        return (uint64_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint64_t*)mxGetUint64s(getField(mx, i, name));
#else
        return (uint64_t*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline int64_t* getInt64s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (int64_t*)mxGetInt64s(mx);
#else
        return (int64_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (int64_t*)mxGetInt64s(getField(mx, i, name));
#else
        return (int64_t*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline uint16_t* getUint16s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint16_t*)mxGetUint16s(mx);
#else
        return (uint16_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint16_t*)mxGetUint16s(getField(mx, i, name));
#else
        return (uint16_t*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline bool* getBools(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (bool*)mxGetLogicals(mx);
#else
        return (bool*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (bool*)mxGetLogicals(getField(mx, i, name));
#else
        return (bool*)mxGetData(getField(mx, i, name));
#endif
    }
}

inline uint8_t* getUint8s(const mxArray* mx, const char* name, mwIndex i = 0) {
    if (std::strcmp(name, "solu") == 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint8_t*)mxGetUint8s(mx);
#else
        return (uint8_t*)mxGetData(mx);
#endif
    }
    else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
        return (uint8_t*)mxGetUint8s(getField(mx, i, name));
#else
        return (uint8_t*)mxGetData(getField(mx, i, name));
#endif
    }
}