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
#include <data/mask.h>
#include <data/metadata_extension.h>

#include "xmipp_gpu_utils.h"
#include "xmipp_gpu_som.h"

#include <math.h>

// Read arguments ==========================================================
void ProgGpuSOM::readParams()
{
    fn_exp = getParam("-i");
    fn_odir = getParam("--odir");

    Niter = getIntParam("--iter");
    somXdim = getIntParam("--somdim",0);
    somYdim = getIntParam("--somdim",1);
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
    ;
}

// usage ===================================================================
void ProgGpuSOM::defineParams()
{

	addParamsLine("   -i <md_exp>                : Metadata file with input experimental images");
    addParamsLine("   [--odir <outdir=\".\"	>]        : Output directory");
	addParamsLine("   [--iter <N=10>]            : Number of iterations");
	addParamsLine("   [--somdim <xdim=10> <ydim=10>]  : Size of the SOM map");
    addUsageLine("Computes a SOM map of size xdim x ydim. A very rough alignment is performed");

}

// Compute correlation --------------------------------------------------------
void ProgGpuSOM::run()
{
}

