/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#include <data/xmipp_image.h>
#include <data/normalize.h>
#include <data/metadata_extension.h>

#include "xmipp_gpu_utils.h"
#include "xmipp_gpu_som.h"
#include <reconstruction_cuda/cuda_gpu_som.h>

// Read arguments ==========================================================
void ProgGpuSOM::readParams()
{
    fn_exp = getParam("-i");
    fn_out = getParam("-o");

    Niter = getIntParam("--iter");
    somXdim = getIntParam("--somdim",0);
    somYdim = getIntParam("--somdim",1);
    normalizeImages = !checkParam("--dontNormalize");
}

// Show ====================================================================
void ProgGpuSOM::show()
{
    std::cout
	<< "Input experimental:             " << fn_exp    << std::endl
	<< "Output directory:               " << fn_out   << std::endl
	<< "Iterations:                     " << Niter     << std::endl
	<< "SOM Xdim:                       " << somXdim   << std::endl
	<< "SOM Ydim:                       " << somYdim   << std::endl
	<< "Normalize images:               " << normalizeImages << std::endl
    ;
}

// usage ===================================================================
void ProgGpuSOM::defineParams()
{

	addParamsLine("   -i <md_exp>                : Metadata file with input experimental images");
    addParamsLine("   [-o <out=\"classes.stk\">]     : Output file");
	addParamsLine("   [--iter <N=10>]            : Number of iterations");
	addParamsLine("   [--somdim <xdim=10> <ydim=10>]  : Size of the SOM map");
	addParamsLine("   [--dontNormalize]          : Don't normalize input images");
    addUsageLine("Computes a SOM map of size xdim x ydim. A very rough alignment is performed");
}

// Produce side info ======================================================
void ProgGpuSOM::readImage(Image<double> &I, size_t objId, bool applyGeo) const
{
    if (applyGeo)
        I.readApplyGeo(SFexp, objId);
    else
    {
        FileName fnImg;
        SFexp.getValue(MDL_IMAGE, fnImg, objId);
        I.read(fnImg);
    }
    if (normalizeImages)
    	I().statisticsAdjust(0, 1);
}

void ProgGpuSOM::produceSideInfo()
{
	SFexp.read(fn_exp);
	size_t Nimgs = SFexp.size();
    size_t somdim = somXdim*somYdim;
    Image<double> I;
	if (fn_out.exists())
	{
		I.read(fn_out);
		MultidimArray<double> Iq;
	    Iref.initZeros(XSIZE(I()), YSIZE(I()), 1, NSIZE(I()));
	    for (size_t q=0; q<somdim; ++q)
	    {
	    	Iq.aliasImageInStack(I(),q);
	    	Iref.fillThisWithImage(q,Iq);
	    }
	}
	else
	{
		std::cout << "Initializing classes ...\n";
		init_progress_bar(Nimgs);
	    MultidimArray<double> Irefq;
	    size_t q=0, idx=0;
	    size_t xdim, ydim, zdim, ndim;
	    getImageSize(fn_exp, xdim, ydim, zdim, ndim);
	    Iref.initZeros(xdim, ydim, 1, somdim);
	    FOR_ALL_OBJECTS_IN_METADATA(SFexp)
	    {
	    	readImage(I, __iter.objId, true);
	    	normalize_OldXmipp(I());
	    	Iref.fillThisWithImage(q,I());
	    	q=(q+1)%somdim;
	    	if (idx%100==0)
	    		progress_bar(idx);
	    	idx++;
	    }
		progress_bar(Nimgs);
	}

	// Normalize the reference images
	I().resizeNoCopy(Iref.ydim,Iref.xdim);
    for (size_t q=0; q<somdim; ++q)
    {
    	Iref.fillImageWithThis(q,I());
    	normalize_OldXmipp(I());
    	Iref.fillThisWithImage(q,I());
    }

    // Transfer to the GPU
	Iref.copyToGpu(Iref_gpu);
	float gpuMemory[3];
	cuda_check_gpu_memory(gpuMemory);

	Nblock = (size_t) floor(gpuMemory[1]/(8*MULTIDIM_SIZE(I())*sizeof(float)));
	Nblock = std::min(Nblock,Nimgs);

	// *** Limit Nblock by the number of gridsize
	std::cout << "GPU Experimental block size: " << Nblock << std::endl;
	Iexp.resize(Iref.xdim,Iref.ydim,1,8*Nblock);
}

// Generate SOM --------------------------------------------------------
//#define DEBUG
void ProgGpuSOM::run()
{
	produceSideInfo();

	// Read the experimental images
	size_t n=0;
	Image<double> I;
    FOR_ALL_OBJECTS_IN_METADATA(SFexp)
    {
    	readImage(I, __iter.objId, true);
    	normalize_OldXmipp(I());
    	Iexp.fillThisWithImage(n,I());
    	n++;
    }
    Iexp.copyToGpu(Iexp_gpu);
    cuda_generate_8rotations(Iexp_gpu);
	Iexp.copyFromGpu(Iexp_gpu);
#ifdef DEBUG
	Iexp.write("PPPgpu.xmp");
#endif

	for (int iter=0; iter<Niter; iter++)
	{
		std::cout << "Iteration " << iter << std::endl;
		if (iter>0)
			Iref.copyToGpu(Iref_gpu);

		// Calculate correlations
		cuda_calculate_correlations(Iexp_gpu,Iref_gpu,cc_gpu);
		size_t nexp=cc_gpu.xdim;
		size_t nrot=cc_gpu.ydim;
		size_t nref=cc_gpu.zdim;
		cc.resize(nexp,nrot,nref,1);
		cc.copyFromGpu(cc_gpu);
#ifdef DEBUG
		cc.write("PPPcc.vol");
#endif

		// For each experimental image calculate the winning node
		winnerRot.resize(nexp);
		winnerRot.initConstant(-1);
		winnerRef=winnerRot;
		winnerCC.resize(nexp);
		winnerCC.initConstant(-2.0);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(cc)
			if (DIRECT_A3D_ELEM(cc,k,i,j)>DIRECT_A1D_ELEM(winnerCC,j) && DIRECT_A3D_ELEM(cc,k,i,j)>0)
			{
				DIRECT_A1D_ELEM(winnerCC,j)=DIRECT_A3D_ELEM(cc,k,i,j);
				DIRECT_A1D_ELEM(winnerRot,j)=(int)i;
				DIRECT_A1D_ELEM(winnerRef,j)=(int)k;
			}

		// Update the winning nodes
		Iref.initZeros();
		refWeights.initZeros(nref);
		size_t xydim=XSIZE(Iref)*YSIZE(Iref);
		size_t skipOneImageRow = nexp*xydim;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(winnerRef)
		{
			int refIdx=DIRECT_A1D_ELEM(winnerRef,i);
			float *ptrFrom=MULTIDIM_ARRAY(Iexp)+DIRECT_A1D_ELEM(winnerRot,i)*skipOneImageRow+i*xydim;
			float *ptrTo=MULTIDIM_ARRAY(Iref)+refIdx*xydim;
			for (size_t ij=0; ij<xydim; ++ij)
				*ptrTo++ += *ptrFrom++;
			DIRECT_A1D_ELEM(refWeights,refIdx)+=DIRECT_A1D_ELEM(winnerCC,i);
		}

		// Normalize by the weight
		for (size_t k=0; k<nref; ++k)
			if (DIRECT_A1D_ELEM(refWeights,k)>0)
			{
				float iWeight=1.0f/DIRECT_A1D_ELEM(refWeights,k);
				float *ptrTo=MULTIDIM_ARRAY(Iref)+k*xydim;
				for (size_t ij=0; ij<xydim; ++ij)
					*ptrTo++ *= iWeight;
			}
	}
    Iref.write(fn_out);
}

