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


class ProtPrime2D(em.ProtClassify2D):
    """
    This protocol wraps the *prime2D* program,
    which is a reference-free 2D alignment/clustering algorithm adopted
    from the prime3D probabilistic ab initio 3D reconstruction algorithm.

    Do not search the origin shifts initially, when the cluster centers are
    of low quality.
    """
    _label = 'prime 2d'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('generateReferences', params.BooleanParam, default=True,
                      label='Generate initial references?',
                      help='If you select *Yes*, you should provide as input '
                           'a SetOfClasses2D to use the averages as reference '
                           'and reclassify the images.\n'
                           'If *No*, then you should provide the as input '
                           'a SetOfParticles and the number of classes.')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      condition='generateReferences',
                      label="Input particles", important=True,
                      help='')

        form.addParam('doCTF', params.BooleanParam, default=True,
                      label='Do CTF correction?',
                      help='')

        form.addParam('numberOfClasses', params.IntParam,
                      condition='generateReferences',
                      label="Number of classes",
                      help="")

        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D',
                      condition='not generateReferences',
                      label="Input 2D classes", important=True,
                      help='')

        form.addParam('maskRadius', params.IntParam, default=-1,
                      label='Particle mask radius (px)',
                      help='')

        form.addParam('originShift', params.FloatParam, default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Search origin shifts (px)',
                      help='')

        form.addParam('lowPassFilter', params.IntParam, default=20,
                      label='Low pass filter (A)')


        form.addParam('doAlign', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Align input particles?',
                      help="")

        form.addParallelSection(threads=4, mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.getUniqueId())
        self._insertFunctionStep('runPrime2DStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        stackFn = self.getParticlesStack()

        # Set to None if doCTF is False, then the ctf information will not be
        # written, neither if the input does not have CTF info.
        ctfFn = self._getExtraPath('ctfparams.txt') if self.doCTF else None

        if self.generateReferences:
            writeSetOfParticles(self.inputParticles.get(), stackFn, None, ctfFn)
        else:
            inputClasses = self.inputClasses.get()
            clsStack = self._getExtraPath('startcavgs.mrc')
            docFn = self._getExtraPath('prime2D_startdoc.txt')
            writeSetOfClasses2D(inputClasses, clsStack, stackFn, docFn, ctfFn)

    def _getCommonArgs(self):
        """ Return common command line argument for programs:
        - simple_prime2D
        """
        inputParticles = self.inputParticles.get()

        # We will run simple_prime2d in the extra folder, so 'particles.mrcs'
        # should be there
        args = " stk=particles.mrcs"
        args += " smpd=%0.3f" % inputParticles.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " lp=%d" % self.lowPassFilter
        args += " ncls=%d" % self.numberOfClasses
        args += " nparts=%d" % self.numberOfMpi
        args += " nthr=%d" % self.numberOfThreads


        return args

    def runPrime2DStep(self):
        args = self._getCommonArgs()

        if self.originShift > 0:
            args += " trs=%d" % self.originShift

        if not self.generateReferences:
            args += " refs=startcavgs.mrc"
            args += " oritab=prime2D_startdoc.txt"

        args += " ctf=%s " % ('yes' if self.doCTF else 'no')

        if self.doCTF:
            args += ' deftab=ctfparams.txt'

        self.runJob(simple.getProgram("prime2D", distr=True), args,
                    cwd=self._getExtraPath())

    def createOutputStep(self):
        classes2D = self._createSetOfClasses2D(self.inputParticles.get())
        self._fillClassesFromIter(classes2D, self.getLastIteration())

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary
    
    def _citations(self):
        return ['Reboul2016']
    
    def _methods(self):
        return []

    # --------------------------- UTILS functions ------------------------------
    def getMaskRadius(self):
        # If mask radius is -1, use half of the particle size
        xdim, _, _ = self.inputParticles.get().getDimensions()

        return xdim / 2 if self.maskRadius < 0 else self.maskRadius.get()

    def getNumberOfClasses(self):
        return self.numberOfClasses.get()

    def getDocFile(self, iteration=None):
        """ Return the document file with alignment parameters for the given
        iteration, if None passed, return the last iteration.
        """
        if iteration is None:
            iteration = self.getLastIteration()

        return self._getExtraPath("prime2D_doc%d.txt" % iteration)

    def getClassesStack(self, iteration=None):
        if iteration is None:
            iteration = self.getLastIteration()

        return self._getExtraPath("cavgs_iter%d.mrc" % iteration)

    def getParticlesStack(self):
        return self._getExtraPath('particles.mrcs')

    def getLastIteration(self):
        lastIter = 1
        while os.path.exists(self.getClassesStack(lastIter)):
            lastIter += 1
        return lastIter - 1

    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        self._iter = iteration
        self._count = 0
        doc = simple.SimpleDocFile(self.getDocFile(iteration))
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=doc.iterValues())
        doc.close()

    def getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getExtraPath('classes_it%03d.sqlite' % it)

        if clean:
            pwutils.cleanPath(data_classes)

        if not os.path.exists(data_classes):
            clsSet = em.SetOfClasses2D(filename=data_classes)
            clsSet.setImages(self.inputParticles.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _updateParticle(self, item, row):
        self._count += 1
        item.setClassId(int(float(row['class'])))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        item.setLocation((self._count, self.getParticlesStack()))

    def _updateClass(self, item):
        item.setAlignment2D()
        fn = self.getClassesStack(self._iter) + ":mrcs"
        item.getRepresentative().setLocation(item.getObjId(), fn)

