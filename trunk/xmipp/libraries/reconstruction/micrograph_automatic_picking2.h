/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
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

#ifndef MICROGRAPH_AUTOMATIC_PICKING_H
#define MICROGRAPH_AUTOMATIC_PICKING_H

#include <data/micrograph.h>
#include <data/mask.h>
#include <data/xmipp_program.h>
#include <reconstruction/fourier_filter.h>
#include <reconstruction/transform_geometry.h>
#include <classification/naive_bayes.h>

/// @defgroup AutomaticPicking Image denoising
/// @ingroup ReconsLibrary
//@{
/* Particle ---------------------------------------------------------------- */
class Particle2
{
public:
    FileName micrograph;       // Micrograph
    int x, y;                  // position in micrograph
    char status;               // rejected=0, selected=1 or moved=2
    Matrix1D<double> vec;      // vector of that particle
    double cost;               // Associated cost

    // Print
    friend std::ostream & operator << (std::ostream &_out, const Particle2 &_p);

    // Read
    void read(std::istream &_in, int _vec_size);
};

/* Automatic particle picking ---------------------------------------------- */
/** Class to perform the automatic particle picking */
class AutoParticlePicking2
{
public:
    static const double __penalization = 10;
    static const double __gray_bins = 8;
    static const double __radial_bins = 16;
    static const double __highpass_cutoff = 0.02;
    static const int    __reduction=2; // Of the piece with respect to the micrograph

    Micrograph                *__m;
    Image<double>              microImage;
    FileName                   fn_micrograph;
    int                        __numThreads;
    Mask                       __mask;
    int                        piece_xsize;
    int                        particle_radius;
    int                        __piece_overlap;
    int                        __scan_overlap;
    int                        __learn_overlap;
    double                       scaleRate;
    bool                       __fast;
    bool                       __incore;
    std::vector<Particle2>      __auto_candidates;
    std::vector<Particle2>      __rejected_particles;
    std::vector < MultidimArray<double> > filterBank;

public:

    /// Empty constructor
    AutoParticlePicking2(const FileName &fn, Micrograph *_m, bool __fast);

    /// Destructor
    ~AutoParticlePicking2();

    /// Read the micrograph in memory
    void readMicrograph();

    /// load models with a name
    void loadModels(const FileName &fn_root);

    /// Save models
    void saveModels(const FileName &fn_root) const;

    /** Save the feature vectors of the automatically selected particles.
        * Nvectors is the number of vectors to save, given by saveAutoParticles.
        */
    void saveAutoFeatureVectors(const FileName &fn, int Nvectors) const;

    /// Load the feature vectors of the automatically selected particles
    void loadAutoFeatureVectors(const FileName &fn);

    /// Training
    void learnParticles();

    /// Build vectors
    void  buildPositiveVectors();
};

class ProgMicrographAutomaticPicking2: public XmippProgram
{
public:
    /// Micrograph filename
    FileName fn_micrograph;
    /// Model rootname
    FileName fn_model;
    /// Training coordinates
    FileName fn_train;
    /// Particle size
    int size;
    /// Mode
    String mode;
    /// Number of threads
    int Nthreads;
    /// Output rootname
    FileName fn_root;
    /// Fast
    bool fast;
    /// In core
    bool incore;
public:
    /// Read parameters
    void readParams();

    /// Show parameters
    void show();

    /// Define Parameters
    void defineParams();

    /** Run */
    void run();
};

//@}
#endif
