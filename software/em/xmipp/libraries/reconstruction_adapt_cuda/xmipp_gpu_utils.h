/***************************************************************************
 *
 * Authors:    Amaya Jimenez      ajimenez@cnb.csic.es (2017)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef _PROG_GPU_UTILS
#define _PROG_GPU_UTILS

#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <reconstruction_cuda/cuda_gpu_correlation.h>
#include <reconstruction_cuda/cuda_xmipp_utils.h>

template<typename T>
class GpuMultidimArrayAtCpu
{
public:
	size_t xdim, ydim, zdim, ndim, yxdim, zyxdim, nzyxdim;
    T* data;

	GpuMultidimArrayAtCpu()
    {
		xdim=ydim=zdim=ndim=yxdim=zyxdim=nzyxdim=0;
		data=NULL;
    }

	GpuMultidimArrayAtCpu(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {
		resize(_Xdim, _Ydim, _Zdim, _Ndim);
    }

	void resize(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {
		xdim=_Xdim;
		ydim=_Ydim;
		zdim=_Zdim;
		ndim=_Ndim;
        yxdim=(size_t)_Ydim*_Xdim;
        zyxdim=yxdim*_Zdim;
        nzyxdim=zyxdim*_Ndim;
        if (nzyxdim>0)
        	data=new T[nzyxdim];
        else
        	data=NULL;
    }

	void initZeros(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
	{
		resize(_Xdim, _Ydim, _Zdim, _Ndim);
		initZeros();
	}

	void initZeros()
	{
		memset(data,0,nzyxdim*sizeof(T));
	}

	void fillImage(size_t n, const MultidimArray<T> &from)
	{
		memcpy(data+n*zyxdim, MULTIDIM_ARRAY(from), MULTIDIM_SIZE(from)*sizeof(T));
	}

	void fillThisWithImage(size_t n, const MultidimArray<double> &Ifrom)
	{
		T *ptrTo=data+n*zyxdim;
		double *ptrFrom=MULTIDIM_ARRAY(Ifrom);
		for (size_t i=0; i<zyxdim; i++)
			*ptrTo++ = (T) *ptrFrom++;
	}

	void fillImageWithThis(size_t n, MultidimArray<double> &Ito)
	{
		T *ptrFrom=data+n*zyxdim;
		double *ptrTo=MULTIDIM_ARRAY(Ito);
		for (size_t i=0; i<zyxdim; i++)
			*ptrTo++ = (double) *ptrFrom++;
	}

	void addImage(size_t n, const MultidimArray<double> &ItoAdd)
	{
		T *ptrTo=data+n*zyxdim;
		double *ptrFrom=MULTIDIM_ARRAY(ItoAdd);
		for (size_t i=0; i<zyxdim; i++)
			*ptrTo++ += (T) *ptrFrom++;
	}

	double sum()
	{
		double retval=0.0;
		T *ptr=data;
		for (size_t i=0; i<nzyxdim; i++)
			retval+=*ptr++;
		return retval;
	}

	bool isEmpty()
	{
		return data==NULL;
	}

	void copyToGpu(GpuMultidimArrayAtGpu<T> &gpuArray)
	{
		if (gpuArray.isEmpty())
			gpuArray.resize(xdim,ydim,zdim,ndim);

		gpuCopyFromCPUToGPU(data, gpuArray.d_data, nzyxdim*sizeof(T));
	}

	void copyFromGpu(GpuMultidimArrayAtGpu<T> &gpuArray)
	{
		gpuCopyFromGPUToCPU(gpuArray.d_data, data, nzyxdim*sizeof(T));
	}

	void copyToGpuMultiple(GpuMultidimArrayAtGpu<T> &gpuArray, int numCopy)
	{
		if (gpuArray.isEmpty())
			gpuArray.resize(xdim,ydim,zdim,numCopy);

		int index = 0;
		for(int i=0; i<numCopy; i++){
			gpuCopyFromCPUToGPU(data, &gpuArray.d_data[index], nzyxdim*sizeof(T));
			index=index+nzyxdim;
		}
	}

	void clear()
	{
		delete []data;
	}

	~GpuMultidimArrayAtCpu()
	{
		clear();
	}

	void write(const FileName &fnOut)
	{
		Image<double> I;
		I().resize(ndim,zdim,ydim,xdim);
		T *ptrFrom=data;
		double *ptrTo=MULTIDIM_ARRAY(I());
		for (size_t i=0; i<nzyxdim; i++)
			*ptrTo++ = (double) *ptrFrom++;
		I.write(fnOut);
	}
};

template<class T>
void fillImage(GpuMultidimArrayAtGpu<T> &to, const MultidimArray<T> &from, size_t n=0)
{
	gpuCopyFromCPUToGPU(MULTIDIM_ARRAY(from), to.d_data+n*MULTIDIM_SIZE(from), MULTIDIM_SIZE(from)*sizeof(T));
}

//@}
#endif






