
//Host includes
#include "cuda_gpu_som.h"
#include <iostream>

__global__ void kernel_generate_8rotations (float *images, size_t xdim, size_t ydim, size_t ndim)
{
	size_t xydim = xdim*ydim;
	size_t ndim_8 = ndim/8;
	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx>=ndim_8*xydim)
		return;
	size_t myImgN  = idx / xydim;
	size_t myPixel = idx % xydim;
	size_t myY = myPixel / xdim;
	size_t myX = myPixel % xdim;
	size_t myYm = ydim-1-myY;
	size_t myXm = xdim-1-myX;

	size_t skipOneImageRow = ndim_8*xydim;
	float *myImg = images + myImgN*xydim;
	float *ptrOut = myImg + myPixel + skipOneImageRow;

	// -x, y
	*ptrOut = * (myImg + (myY*xdim) + myXm);
	ptrOut+=skipOneImageRow;

	// x, -y
	*ptrOut = * (myImg + (myYm*xdim) + myX);
	ptrOut+=skipOneImageRow;

	// -x, -y
	*ptrOut = * (myImg + (myYm*xdim) + myXm);
	ptrOut+=skipOneImageRow;

	// y, x
	*ptrOut = * (myImg + (myX*xdim) + myY);
	ptrOut+=skipOneImageRow;

	// -y, x
	*ptrOut = * (myImg + (myX*xdim) + myYm);
	ptrOut+=skipOneImageRow;

	// y, -x
	*ptrOut = * (myImg + (myXm*xdim) + myY);
	ptrOut+=skipOneImageRow;

	// -y, -x
	*ptrOut = * (myImg + (myXm*xdim) + myYm);
}

void cuda_generate_8rotations(GpuMultidimArrayAtGpu<float> &images)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    images.calculateGridSizeVectorized(blockSize, gridSize);
    gridSize.x = (size_t) ceil(gridSize.x/8.0);
//    std::cout << "GridSize " << gridSize.x << std::endl;
    kernel_generate_8rotations <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
    		(images.d_data, images.xdim, images.ydim, images.ndim);
    waitForKernelToFinish();
}

// ==========================================================
__global__ void kernel_calculate_correlations (float *Iexp, float *Iref, float *cc, size_t xdim, size_t ydim, size_t nexp, size_t nref)
{
	size_t totalNcc = nexp*nref*8;
	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx>=totalNcc)
		return;
	size_t xydim = xdim*ydim;
	size_t nexp8 = nexp*8;
	size_t myIrefIdx = idx/nexp8;
	size_t remaining = idx%nexp8;
	size_t myRotIdx = remaining/nexp;
	size_t myIexpIdx = remaining%nexp;
	size_t skipOneImageRow = nexp*xydim;

	float *myIref=Iref+myIrefIdx*xydim;
	float *myIexp=Iexp+myIexpIdx*xydim+myRotIdx*skipOneImageRow;
	float localcc=0.0;
	for (size_t ij=0; ij<xydim; ++ij)
		localcc+=(*myIref++)*(*myIexp++);
	*(cc+idx)=localcc/xydim;
}

void cuda_calculate_correlations(GpuMultidimArrayAtGpu<float> &Iexp, GpuMultidimArrayAtGpu<float> &Iref, GpuMultidimArrayAtGpu<float> &cc)
{
	cc.resize(Iexp.ndim/8,8,Iref.ndim);
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    cc.calculateGridSizeVectorized(blockSize, gridSize);
    kernel_calculate_correlations<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
    		(Iexp.d_data, Iref.d_data, cc.d_data, Iexp.xdim, Iexp.ydim, Iexp.ndim/8, Iref.ndim);
    waitForKernelToFinish();
}
