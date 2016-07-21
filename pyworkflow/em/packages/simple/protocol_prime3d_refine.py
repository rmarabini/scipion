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
from convert import rowToAlignment, writeSetOfParticles


class ProtPrime3DRefine(em.ProtRefine3D):
    """
    """
    _label = 'prime 3d refine'

    #--------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignment',
                      label="Input particles", important=True,
                      help='')

        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label='Input volume', important=True,
                      help='')

        form.addParam('maskRadius', params.IntParam, default=-1,
                      label='Particle mask radius (px)',
                      help='')

        form.addParam('lowPassFilter', params.IntParam, default=20,
                      label='Low pass filter (A)')

        form.addParam('symmetry', params.StringParam, default='c1',
                      label='Symmetry group',
                      help='Possibilities are: cn|dn|t|o|i{c1}')

        form.addParam('splitData', params.BooleanParam, default=True,
                      label='Split data in halves?',
                      help='')

        form.addParam('extraParams', params.StringParam,
                      label='Extra parameters',
                      help='')

        form.addParallelSection(threads=4, mpi=0)

    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.getUniqueId())
        self._insertFunctionStep('refineStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        stackFn = self.getParticlesStack()
        docFn = stackFn.replace('.mrcs', '.txt')
        ctfFn = stackFn.replace('.mrcs', '_ctf.txt')
        writeSetOfParticles(self.inputParticles.get(), stackFn, docFn, ctfFn,
                            alignType=em.ALIGN_PROJ, applyTransform=False)
        em.ImageHandler().convert(self.inputVolume.get(),
                                  self._getExtraPath('volume.mrc'))

    def refineStep(self):
        # $ nohup distr_simple.pl prg=prime3D stk=groel-stk.spi
        # vol1=sym_recvol_state1msk.spi smpd=1.62 msk=60
        # eo=yes oritab=mapped_ptcls_params.txt npart=8
        args = self._getCommonArgs()
        args += ' vol1=volume.mrc'
        args += ' pgrp=%s' % self.symmetry.get().lower()
        args += ' oritab=particles.txt'
        args += ' eo=%s' % ('yes' if self.splitData else 'no')
        if not self.extraParams.empty():
            args += ' %s' % self.extraParams

        self.runJob("simple_prime3D", args, cwd=self._getExtraPath())

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
        vol = em.Volume(self.getVolFile())
        vol.setSamplingRate(self.inputParticles.get().getSamplingRate())

        partSet = self._createSetOfParticles()
        self._fillParticlesFromIter(partSet)

        self._defineOutputs(outputVol=vol, outputParticles=partSet)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineSourceRelation(self.inputVolume, vol)
        self._defineSourceRelation(self.inputParticles, partSet)
        self._defineSourceRelation(self.inputVolume, partSet)

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

    def getDocFile(self, iteration=None):
        """ Return the document file with alignment parameters for the given
        iteration, if None passed, return the last iteration.
        """
        if iteration is None:
            iteration = self.getLastIteration()

        return self._getExtraPath("prime3D_doc%d.txt" % iteration)

    def getVolFile(self, iteration=None):
        """ Return the document file with alignment parameters for the given
        iteration, if None passed, return the last iteration.
        """
        if iteration is None:
            iteration = self.getLastIteration()

        return self._getExtraPath("recvol_state1_iter%d.mrc" % iteration)

    def getLastIteration(self):
        lastIter = 1
        while os.path.exists(self.getVolFile(lastIter)):
            lastIter += 1
        return lastIter - 1

    def _fillParticlesFromIter(self, partSet, iteration=None):
        """ Create the SetOfClasses2D from a given iteration. """
        inputParticles = self.inputParticles.get()

        doc = simple.SimpleDocFile(self.getDocFile(iteration))

        partSet.copyInfo(inputParticles)
        partSet.setAlignmentProj()
        partSet.copyItems(inputParticles,
                          updateItemCallback=self._updateParticle,
                          itemDataIterator=doc.iterValues())

        doc.close()

    def getIterParticles(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        partSqlite = self._getExtraPath('particles_it%03d.sqlite' % it)

        if clean:
            pwutils.cleanPath(partSqlite)

        if not os.path.exists(partSqlite):
            partSet = em.SetOfParticles(filename=partSqlite)
            self._fillParticlesFromIter(partSet, it)
            partSet.write()
            partSet.close()

        return partSqlite

    def _updateParticle(self, item, row):
        item.setTransform(rowToAlignment(row, em.ALIGN_PROJ))
