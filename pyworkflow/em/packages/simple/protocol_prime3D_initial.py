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

        form.addParam('symmetry', params.StringParam, default='c1',
                      label='Symmetry group',
                      help='Possibilities are: cn|dn|t|o|i{c1}')

        form.addParam('searchSymAxis', params.BooleanParam, default=False,
                      label='Search for the symmetry axis?',
                      help='')

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputSet.getUniqueId())

        sym = self.symmetry.get().lower()
        initSym = 'c1' if self.searchSymAxis else sym
        self._insertFunctionStep('prime3DInitStep', initSym)
        # The following files should be output from the prime 3d - init step
        initVol, initOritab = 'startvol_state1.mrc', 'prime3D_startdoc.txt'
        self._insertFunctionStep('prime3DStep', initSym, initVol, initOritab)

        if self.searchSymAxis:
            self._insertStep('findSymAxisStep', initVol, initOritab)
            # The following files should be produced after finding symmetry axis
            vol, oritab = "", ""
            self._insertStep('prime3DStep', sym, vol, oritab)


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

    def prime3DInitStep(self, initSym):
        args = self._getCommonArgs()
        args += " pgrp=%s" % initSym

        self.runJob("simple_prime3D_init", args, cwd=self._getExtraPath())

    def prime3DStep(self, sym, vol, oritab):
        args = self._getCommonArgs()
        args += " vol1=%s" % vol
        args += " pgrp=%s" % sym
        args += " oritab=%s" % oritab

        self.runJob("simple_prime3D", args, cwd=self._getExtraPath())

    def findSymAxisStep(self, vol, oritab):
        # $ simple_symsrch vol1=prime3D_round_16/recvol_state1.spi smpd=1.62 msk=60
        # oritab=prime3D_round_16/prime3Ddoc_16.txt pgrp=d7 outfile=sym_d7.txt nthr=8
        # lp=20 > SYMOUT
        args = self._getCommonArgs(stk=False)
        args += " vol1=%s" % vol
        args += " oritab=%s" % oritab
        args += " lp=%d" % 20 # FIXME use lowpass filter output from prime3d
        args += " "

    def _getCommonArgs(self, stk=True):
        inputSet = self.inputSet.get()

        args = " stk=particles.mrcs" if stk else ""
        args += " smpd=%f" % inputSet.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " nthr=%d" % self.numberOfThreads

        return args

    def createOutputStep(self):
        return
        vol = em.Volume()
        vol.setFileName()
        vol.setSamplingRate(self.inputSet.get().getSamplingRate())

        self._defineOutputs(outputVol=vol)
        self._defineSourceRelation(self.inputSet, vol)

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

    def getInitOutputs(self):
        return ['startvol_state1.mrc', 'prime3D_startdoc.txt']

