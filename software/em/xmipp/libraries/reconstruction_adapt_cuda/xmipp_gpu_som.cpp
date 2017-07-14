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
#include <data/metadata_extension.h>

#include "xmipp_gpu_utils.h"
#include "xmipp_gpu_som.h"

// Read arguments ==========================================================
void ProgGpuSOM::readParams()
{
    fn_exp = getParam("-i");
    fn_odir = getParam("--odir");

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
	<< "Output directory:               " << fn_odir   << std::endl
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
    addParamsLine("   [--odir <outdir=\".\"	>]        : Output directory");
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
	FileName fnClasses = fn_odir+"/classes.xmd";
	int Nimgs = SFexp.size();
	if (fnClasses.exists())
	{

	}
	else
	{
		std::cout << "Initializing classes ...\n";
		init_progress_bar(Nimgs);
	    Image<double> I;
	    MultidimArray<double> Irefq;
	    size_t q=0, idx=0;
	    size_t somdim = somXdim*somYdim;
	    size_t xdim, ydim, zdim, ndim;
	    getImageSize(fn_exp, xdim, ydim, zdim, ndim);
	    Iref.initZeros(xdim, ydim, 1, somdim);
	    FOR_ALL_OBJECTS_IN_METADATA(SFexp)
	    {
	    	readImage(I, __iter.objId, true);
	    	float *ptrq=Iref.data+Iref.yxdim*q;
	    	double *ptrI=&I(0,0);
	    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
	    	{
	    		*ptrq += (float)*ptrI;
	    		ptrq++;
	    		ptrI++;
	    	}

	    	q=(q+1)%somdim;
	    	if (idx%100==0)
	    		progress_bar(idx);
	    	idx++;
	    }
		progress_bar(Nimgs);
		Iref.copyToGpu(Iref_gpu);
//		Iref.write(fn_odir+"/classes.stk");
	}


}

// Compute correlation --------------------------------------------------------
void ProgGpuSOM::run()
{
	produceSideInfo();
}

