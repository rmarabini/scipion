# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import pyworkflow.em as em

import simple
from convert import rowToAlignment, writeSetOfParticles, writeSetOfClasses2D


class ProtSimpleReconstruct(em.ProtReconstruct3D):
    """
    This protocols wraps *simple_recvol* program for reconstructing volumes
    given input orientations and state assignments.

    The algorithm is based on direct Fourier inversion with a Kaiser-Bessel (KB)
    interpolation kernel. This window function reduces the real-space ripple
    artifacts associated with direct moving windowed-sinc interpolation.
    The feature sought when implementing this algorithm was to enable quick,
    reliable reconstruction from aligned individual particle images.
    """
    _label = 'reconstruct'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='')

        form.addParam('maskRadius', params.IntParam, default=-1,
                      label='Particle mask radius (px)',
                      help='')

        form.addParam('lowPassFilter', params.IntParam, default=20,
                      label='Low pass filter (A)')

        form.addParam('symmetry', params.StringParam, default='c1',
                      label='Symmetry group',
                      help='Possibilities are: cn|dn|t|o|i{c1}')

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.getUniqueId())
        self._insertFunctionStep('reconstructStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        stackFn = self.getParticlesStack()
        docFn = stackFn.replace('.mrcs', '.txt')
        ctfFn = self._getExtraPath('ctfparams.txt')
        writeSetOfParticles(self.inputParticles.get(), stackFn, docFn, ctfFn,
                            alignType=em.ALIGN_PROJ)

    def reconstructStep(self):
        inputParticles = self.inputParticles.get()

        # simple_recvol stk=<ptcls.ext> smpd=<sampling distance(in A)>
        # oritab=<algndoc.txt> msk=<mask radius(in pixels)>
        # [lp=<low-pass limit{20}>] [frac=<fraction of ptcls to include{1.}>]
        # [nthr=<nr of openMP threads{1}>] [pgrp=<cn|dn|t|o|i{c1}>]
        args = self._getCommonArgs()
        args += ' pgrp=%s' % self.symmetry.get().lower()
        args += ' oritab=particles.txt'

        self.runJob("simple_recvol", args, cwd=self._getExtraPath())

    def _getCommonArgs(self):
        """ Return common command line argument for programs:
        - simple_prime2D_init
        - simple_prime2D
        """
        inputParticles = self.inputParticles.get()

        # We will run simple_prime2d in the extra folder, so 'particles.mrcs'
        # should be there
        args = " stk=particles.mrcs"
        args += " smpd=%f" % inputParticles.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " lp=%d" % self.lowPassFilter
        args += " nthr=%d" % self.numberOfThreads

        return args

    def createOutputStep(self):
        volFile = self._getExtraPath('recvol_state1.mrc')
        vol = em.Volume(volFile)
        vol.setSamplingRate(self.inputParticles.get().getSamplingRate())

        self._defineOutputs(outputVol=vol)
        self._defineSourceRelation(self.inputParticles, vol)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []

    # --------------------------- UTILS functions ------------------------------
    def getMaskRadius(self):
        # If mask radius is -1, use half of the particle size
        xdim, _, _ = self.inputParticles.get().getDimensions()

        return xdim / 2 if self.maskRadius < 0 else self.maskRadius.get()

    def getParticlesStack(self):
        return self._getExtraPath('particles.mrcs')

