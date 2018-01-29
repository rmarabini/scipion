/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#include "reconstruction_adapt_cuda/movie_alignment_correlation_gpu.h"


// FIXME: REMOVE
#include <sstream>
#include "data/filters.h"
#include "data/xmipp_fftw.h"
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
// FIXME

void ProgMovieAlignmentCorrelationGPU::loadFrame(const MetaData& movie, size_t objId, bool crop, Image<float>& out) {
	FileName fnFrame;
	movie.getValue(MDL_IMAGE, fnFrame, objId);
	if (crop) {
		Image<double>tmp;
		tmp.read(fnFrame);
		tmp().window(out(), yLTcorner, xLTcorner, yDRcorner, xDRcorner);
	} else {
		out.read(fnFrame);
	}
}

void ProgMovieAlignmentCorrelationGPU::loadData(const MetaData& movie,
		const Image<double>& dark, const Image<double>& gain,
		double targetOccupancy, const MultidimArray<double>& lpf) {
	Image<float> frame, gainF, darkF;
	// FIXME consider loading imgs in full size and cropping them on GPU
	bool cropInput = (yDRcorner != -1);
	// find input image size
	gainF.data.resize(gain());
	darkF.data.resize(dark());

	loadFrame(movie, movie.firstObject(), cropInput, frame);
	int noOfImgs = nlast - nfirst + 1;
	size_t noOfFloats = noOfImgs * std::max(frame.data.yxdim, (frame.data.xdim/2+1) * frame.data.ydim * 2);
	float* imgs = new float[noOfFloats]();
	std::cout << "noOfFloats: " << noOfFloats << std::endl;
	int counter = -1;
	int paddedLineLength = (frame.data.xdim/2+1)*2;
	FOR_ALL_OBJECTS_IN_METADATA(movie) {
		counter++;
		if (counter < nfirst ) continue;
		if (counter > nlast) break;

		loadFrame(movie, __iter.objId, cropInput, frame);
		// FIXME optimize if necessary
		if (XSIZE(darkF()) > 0)
			frame() -= darkF();
		if (XSIZE(gainF()) > 0)
			frame() *= gainF();


//		************
//		IN-PLACE
//		************
		// copy line by line, adding offset at the end of each line
		// result is the same image, padded in the X dimension to (N/2+1)*2
		float* dest = imgs + ((counter-nfirst) * paddedLineLength * frame.data.ydim); // points to first float in the image
		for (int i = 0; i < frame.data.ydim; i++) {
			memcpy(dest + (paddedLineLength * i),
					frame.data.data + i*frame.data.xdim,
					frame.data.xdim * sizeof(float));
		}

//		************
//		OUT-OF-PLACE
//		************
//		// add image at the end of the stack (that is already long enough)
//		memcpy(imgs + ((counter-nfirst) * ((frame.data.xdim/2+1) * frame.data.ydim * 2)),
//				frame.data.data,
//				frame.data.yxdim * sizeof(float));
	}

	Image<float> aaaa((frame.data.xdim/2+1)*2, frame.data.ydim, 1, noOfImgs);
	aaaa.data.data = imgs;
	aaaa.write("images.vol");

	std::complex<float>* result;
	kernel1(imgs, frame.data.xdim, frame.data.ydim, noOfImgs, newXdim, newYdim, result);
	// 16785408 X:2049 Y:4096
	Image<float> tmp(2049, 4096, 1, noOfImgs);
	for (size_t i = 0; i < 16785408L; i++) {
//	for (size_t i = 0; i < 8388608L; i++) {
		float val = result[i].real() / 16785408.f;
		if (val < 1) tmp.data[i] = val;
		else std::cout << "skipping " << val << " at position " << i << std::endl;

	}
	tmp.write("fftFromGPU" + SSTR(counter) + ".vol");
}

void ProgMovieAlignmentCorrelationGPU::computeShifts(size_t N,
		const Matrix1D<double>& bX, const Matrix1D<double>& bY,
		const Matrix2D<double>& A) {
	return;
	// FIXME refactor

	int idx = 0;
	MultidimArray<double> Mcorr;
	Mcorr.resizeNoCopy(newYdim, newXdim);
	Mcorr.setXmippOrigin();
	CorrelationAux aux;
	for (size_t i = 0; i < N - 1; ++i) {
		for (size_t j = i + 1; j < N; ++j) {
			bestShift(*frameFourier[i], *frameFourier[j], Mcorr, bX(idx),
					bY(idx), aux, NULL, maxShift);
			if (verbose)
				std::cerr << "Frame " << i + nfirst << " to Frame "
						<< j + nfirst << " -> (" << bX(idx) << "," << bY(idx)
						<< ")\n";
			for (int ij = i; ij < j; ij++)
				A(idx, ij) = 1;

			idx++;
		}
		delete frameFourier[i];
	}
}
