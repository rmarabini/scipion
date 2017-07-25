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
#ifndef _PROG_GPU_SOM
#define _PROG_GPU_SOM

#include <data/xmipp_program.h>
#include "xmipp_gpu_utils.h"

class ProgGpuSOM: public XmippProgram
{
public:
	FileName fn_exp, fn_out;
	int somXdim, somYdim;
	int Niter;
	bool normalizeImages;
	double sideWeight;

public:
    //Input metadata file
    MetaData SFref, SFexp;

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Read image */
    void readImage(Image<double> &I, size_t objId, bool applyGeo) const;

    /** Produce side info */
    void produceSideInfo();

    /** processImage */
    void run();

public:
    GpuMultidimArrayAtCpu<float> Iref, Iexp, cc;
    GpuMultidimArrayAtGpu<float> Iref_gpu, Iexp_gpu, cc_gpu;
    size_t Nblock;
    MultidimArray<int> winnerRot, winnerRef;
    MultidimArray<float> winnerCC, refWeights;
};
//@}
#endif
