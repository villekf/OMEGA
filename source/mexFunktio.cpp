#include "mexFunktio.h"
#ifndef OCTAVE

const bool getScalarBool(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsInt8(mx) && !mxIsLogical) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type logical.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const bool scalar = (bool)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const int8_t getScalarInt8(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsInt8(mx) && !mxIsLogical) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type int8.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const int8_t scalar = (int8_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const uint8_t getScalarUInt8(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsUint8(mx) && !mxIsLogical) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type uint8.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const uint8_t scalar = (uint8_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const int16_t getScalarInt16(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsInt16(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type int16.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const int16_t scalar = (int16_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const uint16_t getScalarUInt16(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsUint16(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type uint16.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const uint16_t scalar = (uint16_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const int32_t getScalarInt32(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsInt32(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type int32.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const int32_t scalar = (int32_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const uint32_t getScalarUInt32(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsUint32(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type uint32.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const uint32_t scalar = (uint32_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}


const int64_t getScalarInt64(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsInt64(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type int64.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const int64_t scalar = (int64_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const uint64_t getScalarUInt64(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsUint64(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type uint64.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const uint64_t scalar = (uint64_t)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const float getScalarFloat(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsSingle(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type single.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const float scalar = (float)mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

const double getScalarDouble(const mxArray* mx, const int ind) {
    // Check input
    //if (!mxIsDouble(mx)) {
    //    char ch[40];
    //    std::snprintf(ch, 40, "Scalar number %d is not of type double.\n", ind);
    //    mexErrMsgTxt(ch);
    //}

    if (mxGetNumberOfElements(mx) == 1) {
        const double scalar = mxGetScalar(mx);
        return scalar;
    }
    else {
        mexErrMsgTxt("Scalar is not of size == [1 1].\n");
        return false;
    }
}

#endif