
#include <cuda_runtime_api.h>
#include "reconstruction_cuda/cuda_utils.h" // cannot be in header as it includes cuda headers
#include "cuda_gpu_reconstruct_fourier.h"
#include "reconstruction_cuda/cuda_basic_math.h"

#define BLOCK_DIM_X 32

// run per each pixel of the dest array
__global__
void kernel2(const float2* __restrict__ src, float2* dest, int noOfImages, size_t oldX, size_t oldY, size_t newX, size_t newY,
		const float* __restrict__ filter) {
	// assign pixel to thread
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;

	if (idx == 0 && idy ==0) {
		printf("kernle2 called %p %p %d old:%lu %lu new:%lu %lu filter %p\n", src, dest, noOfImages, oldX, oldY, newX, newY, filter);
	}
	if (idx >= newX || idy >= newY ) return;

	int yhalf = (newY+1)/2;

	size_t fIndex = idy*newX + idx; // index within single image
	size_t origY = (idy <= yhalf) ? idy : (oldY - (newY-idy)); // take top N/2+1 and bottom N/2 lines
	for (int n = 0; n < noOfImages; n++) {
		size_t iIndex = n*oldX*oldY + origY*oldX + idx; // index within consecutive images
		size_t oIndex = n*newX*newY + fIndex; // index within consecutive images
//		if (iIndex >= 16785408l || iIndex < 0 || oIndex >= 5177900l || oIndex < 0) {
//			printf("problem: %p %p old:%lu %lu new:%lu %lu : i:%lu o:%lu\nyhalf: %d origY %lu thread %d %d \n", src, dest, oldX, oldY, newX, newY, iIndex, oIndex,
//					yhalf, origY, idx, idy);
//		}
//		if (fIndex >= 2588950) {
//			printf("problem: %p %p old:%lu %lu new:%lu %lu : i:%lu o:%lu f:%lu \nyhalf: %d origY %lu thread %d %d \n", src, dest, oldX, oldY, newX, newY, iIndex, oIndex, fIndex,
//								yhalf, origY, idx, idy);
//		}
		dest[oIndex] = src[iIndex] * filter[fIndex];
	}

//	int halfY = iSizeY / 2;
//	float normFactor = iSizeY*iSizeY;
//	int oSizeX = oBuffer->fftSizeX;
//
//	// input is an image in Fourier space (not normalized)
//	// with low frequencies in the inner corners
//	for (int n = 0; n < iLength; n++) {
//		float2 freq;
//		if ((idy < iSizeY) // for all input lines
//				&& (idx < oSizeX)) { // for all output pixels in the line
//			// process line only if it can hold sufficiently high frequency, i.e. process only
//			// first and last N lines
//			if (idy < oSizeX || idy >= (iSizeY - oSizeX)) {
//				// check the frequency
//				freq.x = FFT_IDX2DIGFREQ(idx, iSizeY);
//				freq.y = FFT_IDX2DIGFREQ(idy, iSizeY);
//				if ((freq.x * freq.x + freq.y * freq.y) > maxResolutionSqr) {
//					continue;
//				}
//				// do the shift (lower line will move up, upper down)
//				int newY = (idy < halfY) ? (idy + oSizeX) : (idy - iSizeY + oSizeX);
//				int oIndex = newY*oSizeX + idx;
//
//				int iIndex = n*iSizeY*iSizeX + idy*iSizeX + idx;
//				float* iValue = (float*)&(iFouriers[iIndex]);
//
//				// copy data and perform normalization
//				oBuffer->getNthItem(oBuffer->FFTs, n)[2*oIndex] = iValue[0] / normFactor;
//				oBuffer->getNthItem(oBuffer->FFTs, n)[2*oIndex + 1] = iValue[1] / normFactor;
//			}
//		}
//	}
}




void kernel1(float* imgs, size_t oldX, size_t oldY, int noOfImages, size_t newX, size_t newY,
		float* filter,
		std::complex<float>*& result) {
//		float*& result) {

//	FIXME Assert newX <= oldX. same with Y


	size_t noOfFloats = noOfImages * std::max(oldX*oldY, (oldX/2+1) * oldY * 2);
	float* d_imgs;
	gpuMalloc((void**) &d_imgs,noOfFloats*sizeof(float));
	gpuErrchk(cudaMemcpy(d_imgs, imgs, noOfFloats*sizeof(float), cudaMemcpyHostToDevice));
	// store to proper structure
	GpuMultidimArrayAtGpu<float> imagesGPU(oldX, oldY, 1, noOfImages, d_imgs);
//	imagesGPU.copyToGpu(imgs);

//	************
//	IN-OF-PLACE
//	************
	GpuMultidimArrayAtGpu<std::complex<float> > resultingFFT(imagesGPU.Xdim / 2 + 1,
			imagesGPU.Ydim,
			imagesGPU.Zdim,
			imagesGPU.Ndim,
			(std::complex<float>*)imagesGPU.d_data);

//	************
//	OUT-OF-PLACE
//	************
//	GpuMultidimArrayAtGpu<std::complex<float> > resultingFFT;

// perform FFT
	mycufftHandle myhandle;
	std::cout << "about to do FFT" << std::endl;
	imagesGPU.fft(resultingFFT, myhandle);
//	myhandle.clear(); // release unnecessary l || oIndex < 0) {
	//			printf("problem: %p %p old:%lu %lu new:%lu %lu : i:%lu o:%lu\nyhalf: %d origY %lu thread %d %d \n", src, dest, oldX, oldY, newX, newY, iIndex, oIndex,
	//					yhalf, origY, idx, idy);
	//		}memory
	std::cout << "FFT done" << std::endl;

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );


	// crop FFT
	float* d_cropped;
	size_t newFFTX = newX / 2 + 1;
	size_t noOfCroppedFloats = noOfImages * newFFTX * newY * 2; // complex

	// copy filter
	float* d_filter;
	gpuMalloc((void**) &d_filter,newFFTX * newY*sizeof(float));
	gpuErrchk(cudaMemcpy(d_filter, filter, newFFTX * newY*sizeof(float), cudaMemcpyHostToDevice));

	gpuMalloc((void**) &d_cropped,noOfCroppedFloats*sizeof(float));
	cudaMemset(d_cropped, 0.f, noOfCroppedFloats*sizeof(float));
	dim3 dimBlock(BLOCK_DIM_X, BLOCK_DIM_X);
	dim3 dimGrid(ceil(newFFTX/(float)dimBlock.x), ceil(newY/(float)dimBlock.y));
	kernel2<<<dimGrid, dimBlock>>>((float2*)resultingFFT.d_data,(float2*) d_cropped, noOfImages, resultingFFT.Xdim, resultingFFT.Ydim, newFFTX, newY, d_filter);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	gpuErrchk( cudaPeekAtLastError() );

// copy out results
	std::cout << "about to copy to host" << std::endl;
	result = new std::complex<float>[noOfImages*newFFTX*newY]();
//	printf("result: %p\nFFTs: %p\n", result, resultingFFT.d_data );
//	resultingFFT.copyToCpu(result);
	printf ("about to copy to host: %p %p %d\n", result, d_cropped, noOfCroppedFloats*sizeof(float));
	gpuErrchk(cudaMemcpy((void*)result, (void*)d_cropped, noOfCroppedFloats*sizeof(float), cudaMemcpyDeviceToHost));
	std::cout << "copy to host done" << std::endl;
	resultingFFT.d_data = NULL; // unbind connection
//	std::cout << "No of elems: " << resultingFFT.nzyxdim  << " X:" << resultingFFT.Xdim << " Y:" << resultingFFT.Ydim<< std::endl;
}
