
#include <cuda_runtime_api.h>
#include "reconstruction_cuda/cuda_utils.h" // cannot be in header as it includes cuda headers
#include "cuda_gpu_reconstruct_fourier.h"






void kernel1(float* imgs, size_t oldX, size_t oldY, int noOfImages, int newX, int newY,
		std::complex<float>*& result) {
	// store to proper structure
	GpuMultidimArrayAtGpu<float> imagesGPU(oldX, oldY, 1, noOfImages);
	imagesGPU.copyToGpu(imgs);
	// perform FFT
	GpuMultidimArrayAtGpu<std::complex<float> > resultingFFT;
	mycufftHandle myhandle;
	std::cout << "about to do FFT" << std::endl;
	imagesGPU.fft(resultingFFT, myhandle);
	myhandle.clear(); // release unnecessary memory
	std::cout << "FFT done" << std::endl;

	std::cout << "about to copy to host" << std::endl;
	result = new std::complex<float>[resultingFFT.nzyxdim]();
	printf("result: %p\nFFTs: %p", result, resultingFFT.d_data );
	resultingFFT.copyToCpu(result);
	std::cout << "copy to host done" << std::endl;
	std::cout << "No of elems: " << resultingFFT.nzyxdim  << " X:" << resultingFFT.Xdim << " Y:" << resultingFFT.Ydim<< std::endl;
}
