// -*-c++-*-
#ifndef STORMM_HPC_TUTORIAL_I_H
#define STORMM_HPC_TUTORIAL_I_H

#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"

using stormm::card::GpuDetails;

//-------------------------------------------------------------------------------------------------
// Wrapper for summation kernels of various data types
//
// Arguments:
//   vdata:    The array of data to sum
//   n:        The trusted length of vdata
//   vresult:   Separate array to hold the result of the sum
//   ct_data:  Type ID code to indicat ethe data type, e.g. int, of both vdata and vresult
//-------------------------------------------------------------------------------------------------
void wrapTheSummationLaunch(const void* vdata, const size_t n, void* vresult, const size_t ct_data,
                            const GpuDetails &gpu);

#endif
