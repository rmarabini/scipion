
#include <cuda_runtime_api.h>
#include "reconstruction_cuda/cuda_utils.h" // cannot be in header as it includes cuda headers
#include "cuda_gpu_reconstruct_fourier.h"






void kernel1(float* imgs, size_t oldX, size_t oldY, int noOfImages, int newX, int newY,
		std::complex<float>*& result) {

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
	myhandle.clear(); // release unnecessary memory
	std::cout << "FFT done" << std::endl;

	std::cout << "about to copy to host" << std::endl;
	result = new std::complex<float>[resultingFFT.nzyxdim]();
	printf("result: %p\nFFTs: %p\n", result, resultingFFT.d_data );
	resultingFFT.copyToCpu(result);
	std::cout << "copy to host done" << std::endl;
	resultingFFT.d_data = NULL; // unbind connection
	std::cout << "No of elems: " << resultingFFT.nzyxdim  << " X:" << resultingFFT.Xdim << " Y:" << resultingFFT.Ydim<< std::endl;
}
