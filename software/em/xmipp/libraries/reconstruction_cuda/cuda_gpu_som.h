
#ifndef CUDA_GPU_SOM_H
#define CUDA_GPU_SOM_H
#include "cuda_xmipp_utils.h"

void cuda_generate_8rotations(GpuMultidimArrayAtGpu<float> &images);
void cuda_calculate_correlations(GpuMultidimArrayAtGpu<float> &Iexp, GpuMultidimArrayAtGpu<float> &Iref, GpuMultidimArrayAtGpu<float> &cc);
#endif
