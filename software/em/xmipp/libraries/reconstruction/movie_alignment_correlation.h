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

#ifndef _PROG_MOVIE_ALIGNMENT_CORRELATION
#define _PROG_MOVIE_ALIGNMENT_CORRELATION

#include <data/xmipp_program.h>

/**@defgroup MovieAlignmentCorrelation Movie alignment by correlation
   @ingroup ReconsLibrary */
//@{

/** Movie alignment correlation Parameters. */
class ProgMovieAlignmentCorrelation: public XmippProgram
{
public:
    /** Filename of movie metadata */
    FileName fnMovie;
    /** Correction images */
    FileName fnDark, fnGain;
    /** Max shift */
    double maxShift;
    /** Sampling rate */
    double Ts;
    /** Max freq. */
    double maxFreq;
    /** Solver iterations */
    int solverIterations;
    /** Aligned movie */
    FileName fnAligned;
    /** Aligned micrograph */
    FileName fnAvg;
    /** Aligned micrograph */
    FileName fnInitialAvg;
    /** Metadata with shifts */
    FileName fnOut;
    /** First and last frame*/
    int nfirst, nlast;
    /** First and last frame*/
    int nfirstSum, nlastSum;
    /** Do not calculate and use the input shifts */
    bool useInputShifts;
    /** Binning factor */
    double bin;
    /** Bspline order */
    int BsplineOrder;
    /** Outside mode */
    int outsideMode;
    /** Outside value */
    double outsideValue;

    /*****************************/
    /** crop corner **/
    /*****************************/
    /** x left top corner **/
    int xLTcorner;
    /** y left top corner **/
    int yLTcorner;
    /** x right down corner **/
    int xDRcorner;
    /** y right down corner **/
    int yDRcorner;

public:
    // Fourier transforms of the input images
	std::vector< MultidimArray<std::complex<double> > * > frameFourier;

	// Target sampling rate
	double newTs;

	// Target size of the frames
	int newXdim, newYdim;
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Run
    void run();


private:
	int findReferenceImage(size_t N, const Matrix1D<double>& shiftX,
			const Matrix1D<double>& shiftY);
	void loadData(const MetaData& movie, const Image<double>& dark,
			const Image<double>& gain,
			double targetOccupancy,
			const MultidimArray<double>& lpf);
	void solveEquationSystem(Matrix1D<double>& bX, Matrix1D<double>& bY,
			Matrix2D<double>& A, Matrix1D<double>& shiftX,
			Matrix1D<double>& shiftY);
	void loadDarkCorrection(Image<double>& dark);
	void loadGainCorrection(Image<double>& gain);
	void computeShifts(size_t N, const Matrix1D<double>& bX,
			const Matrix1D<double>& bY, const Matrix2D<double>& A);
	void constructLPF(double targetOccupancy, const MultidimArray<double>& lpf);
	void setNewDimensions(double& targetOccupancy, const MetaData& movie,
			double& sizeFactor);
	void readMovie(MetaData& movie);
	void storeRelativeShifts(int bestIref, const Matrix1D<double>& shiftX,
			const Matrix1D<double>& shiftY, double sizeFactor, MetaData& movie);
	void setZeroShift(MetaData& movie);
	int findShiftsAndStore(MetaData& movie, Image<double>& dark,
			Image<double>& gain);
	void applyShiftsComputeAverage(const MetaData& movie,
			const Image<double>& dark, const Image<double>& gain,
			Image<double>& initialMic, size_t& Ninitial,
			Image<double>& averageMicrograph, size_t& N);
	void storeResults(Image<double> initialMic, size_t Ninitial,
			Image<double> averageMicrograph, size_t N, const MetaData& movie,
			int bestIref);
	void correctLoopIndices(const MetaData& movie);
};
//@}
#endif
