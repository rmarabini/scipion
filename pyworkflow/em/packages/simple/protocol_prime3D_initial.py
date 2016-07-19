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


class ProtPrime3DInitial(em.ProtInitialVolume):
    """
    This protocol wraps *simple_prime3D*, an ab inito reconstruction/refinement
    program based on probabilistic projection matching.
    """
    _label = 'prime 3d initial'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfParticles,SetOfClasses,SetOfAverages',
                      label="Input set", important=True,
                      help='You can selected the following type of sets:'
                           'SetOfPartices, SetOfClasses or SetOfAverages '
                           'in order to produce an initial volume. ')

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
                                 self.inputSet.getUniqueId())
        self._insertFunctionStep('reconstructStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        inputSet = self.inputSet.get()
        stkFn = self._getExtraPath('particles.mrcs')

        if isinstance(inputSet, em.SetOfAverages):
            inputSet.writeStack(stkFn)
        elif isinstance(inputSet, em.SetOfClasses2D):
            writeSetOfClasses2D(inputSet, stkFn,
                                stackFn=None, docFn=None, ctfFn=None)
        elif isinstance(inputSet, em.SetOfParticles):
            writeSetOfParticles(inputSet, stkFn, docFn=None, ctfFn=None)
        else:
            raise Exception('Unexpected input type: %s' % type(inputSet))

    def reconstructStep(self):
        # simple_prime3D_init stk=selected_cavgs.spi smpd=1.62 msk=60 nthr=8 lp=20
        args = self._getCommonArgs()

        self.runJob("simple_prime3D_init", args, cwd=self._getExtraPath())

    def _getCommonArgs(self):
        inputSet = self.inputSet.get()

        # We will run simple_prime2d in the extra folder, so 'particles.mrcs'
        # should be there
        args = " stk=particles.mrcs"
        args += " smpd=%f" % inputSet.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " lp=%d" % self.lowPassFilter
        args += " nthr=%d" % self.numberOfThreads

        return args

    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary
    
    def _citations(self):
        return ['Elmlund2013']
    
    def _methods(self):
        return []

    # --------------------------- UTILS functions ------------------------------
    def getMaskRadius(self):
        # If mask radius is -1, use half of the particle size
        xdim, _, _ = self.inputSet.get().getDimensions()

        return self.maskRadius.get() if self.maskRadius < 0 else xdim / 2

