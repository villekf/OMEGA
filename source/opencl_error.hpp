/**************************************************************************
* Header for OpenCL error list
*
* Copyright(C) 2019  Ville - Veikko Wettenhovi
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
//#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_CL_1_2_DEFAULT_BUILD
//#if defined(__APPLE__) || defined(__MACOSX)
//#include <OpenCL/cl.h>
//#else
//#include <CL/cl.h>
//#endif
//#include <CL/cl.hpp>
#include "cl2.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "mexFunktio.h"

#define TH 100000000000.f
#define TH32 100000.f

#define getErrorString(ans) { gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cl_int code, const char* file, int line);

//const char *getErrorString(cl_int error);

std::string header_to_string(const char* header_directory);