/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#ifndef _PROG_CLASSIFY_SIGNIFICANT
#define _PROG_CLASSIFY_SIGNIFICANT

#include <data/xmipp_program.h>
#include "fourier_projection.h"
#include "fourier_filter.h"

/**@defgroup ClassifySignificant Classify a set of images into a discrete set of classes
   @ingroup ReconsLibrary */
//@{

/** Classify Significant Parameters. */
class ProgClassifySignificant: public XmippProgram
{
public:
    /** Filename of the reference volumes */
    FileName fnVols;
    /** Filename of indexes to study */
    FileName fnIds;
    /** Filename of angles assigned */
    FileName fnAngles;
    /** Output file */
    FileName fnOut;
    /** Padding factor */
    int pad;
public:
    // Fourier projector
    std::vector<FourierProjector *> projector;
	// Theoretical projection
	std::vector<Projection *> P;
	// Set of Ids
	std::vector<size_t> setIds;
	// Set of Angles
	std::vector<VMetaData> setAngles;
	// Current row
	std::vector<size_t> currentRowIdx;
	// Set of Angles for a particular image
	std::vector<VMetaData> subsetAngles;
public:
    /// Destructor
    ~ProgClassifySignificant();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void produceSideInfo();

    /** Predict angles and shift.
        At the input the pose parameters must have an initial guess of the
        parameters. At the output they have the estimated pose.*/
    void run();

    void selectSubset(size_t particleId);
};
//@}
#endif
