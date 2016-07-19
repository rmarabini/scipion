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
    This protocol wraps the *simple_prime2D* program,
    which is a reference-free 2D alignment/clustering algorithm adopted
    from the prime3D probabilistic ab initio 3D reconstruction algorithm.

    It is assumed that the images are phase-flipped
    (phase flipping can be done inside the protocol using simple_stackops).

    Do not search the origin shifts initially, when the cluster centers are
    of low quality. If your images are far off centre, you can select to
    pre-align them (internally using simple_stackops with option shalgn=yes).
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
                      label='Search origin shifts (px)',
                      help='')

        form.addParam('lowPassFilter', params.IntParam, default=20,
                      label='Low pass filter (A)')

        form.addParam('doPhaseFlip', params.BooleanParam, default=True,
                      label='Phase flip input particles?',
                      help='Phase flip the input particles if they contains '
                           'CTF information and have not already been flipped.')

        form.addParam('doAlign', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Align input particles?',
                      help="")

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.getUniqueId())

        if self.doPhaseFlip:
            self._insertFunctionStep('phaseFlipStep')

        if self.doAlign:
            self._insertFunctionStep('alignStep')

        if self.generateReferences:
            self._insertFunctionStep('runPrime2DInitStep')

        self._insertFunctionStep('runPrime2DStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        stackFn = self.getParticlesStack()
        ctfFn = self._getExtraPath('ctfparams.txt')

        if self.generateReferences:
            writeSetOfParticles(self.inputParticles.get(), stackFn, None, ctfFn)
        else:
            inputClasses = self.inputClasses.get()
            clsStack = self._getExtraPath('startcavgs.mrc')
            docFn = self._getExtraPath('prime2D_startdoc.txt')
            writeSetOfClasses2D(inputClasses, clsStack, stackFn, docFn, ctfFn)
            # SIMPLE_DOC2CAVGS stk=<stack.ext> smpd=<sampling distance(in A)>
            # msk=<mask radius(in pixels)> ncls=<nr of clusters>
            # oritab=<previous clustering doc> which_iter=<iteration nr>
            # [mul=<shift multiplication factor{1}>] [nthr=<nr of OpenMP threads{1}>]
            # args = " stk=particles.mrcs "
            # args += " oritab=prime2D_startdoc.txt "
            # args += " smpd=%f" % inputClasses.getFirstItem().getSamplingRate()
            # args += " ncls=%d" % inputClasses.getSize()
            # args += " msk=%d" % self.getMaskRadius()
            # args += " nthr=%d" % self.numberOfThreads
            # args += " which_iter=1"
            # self.runJob("simple_doc2cavgs", args, cwd=self._getExtraPath())

    def phaseFlipStep(self):
        inputParticles = self.inputParticles.get()

        if not inputParticles.hasCTF():
            self.info('Input particles does not have CTF information. '
                      'NOT phase flipping.')
            return

        if inputParticles.isPhaseFlipped():
            self.info('Input particles are already phase flipped. '
                      'NOT phase flipping again.')
            return

        # simple_stackops stk=ptcls.mrc smpd=2 deftab=ctfparams.txt
        # ctf=flip kv=300 cs=2.7 fraca=0.07 outstk=ptcls_phflip.mrc
        acq = inputParticles.getAcquisition()
        outputName = "particles_phflip.mrcs"

        args = " stk=particles.mrcs deftab=ctfparams.txt ctf=flip"
        args += " smpd=%f" % inputParticles.getSamplingRate()
        args += " kv=%f" % acq.getVoltage()
        args += " cs=%f" % acq.getSphericalAberration()
        args += " fraca=%f" % acq.getAmplitudeContrast()
        args += " outstk=%s" % outputName

        self.info("Phase flipping input particles.")
        self.runJob("simple_stackops", args, cwd=self._getExtraPath())

        # Replace the initial converted stack with the phase flipped one
        pwutils.moveFile(self._getExtraPath(outputName),
                         self.getParticlesStack())

    def alignStep(self):
        inputParticles = self.inputParticles.get()

        # $ simple_stackops stk=particles.spi smpd=1.62 msk=60
        # shalgn=yes trs=3.5 lp=20 nthr=8 outstk=particles_sh.spi
        acq = inputParticles.getAcquisition()
        outputName = "particles_sh.mrcs"

        args = " stk=particles.mrcs shalgn=yes "
        args += " smpd=%f" % inputParticles.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " lp=%d" % self.lowPassFilter
        args += " trs=3.5" #FIXME
        args += " nthr=%d" % self.numberOfThreads
        args += " outstk=%s" % outputName

        self.info("Aligning input particles.")
        self.runJob("simple_stackops", args, cwd=self._getExtraPath())

        # Replace the initial converted stack with the phase flipped one
        pwutils.moveFile(self._getExtraPath(outputName),
                         self.getParticlesStack())

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
        args += " ncls=%d" % self.numberOfClasses
        args += " nthr=%d" % self.numberOfThreads

        return args

    def runPrime2DInitStep(self):
        self.runJob("simple_prime2D_init", self._getCommonArgs(),
                    cwd=self._getExtraPath())

    def runPrime2DStep(self):
        args = self._getCommonArgs()
        if self.originShift > 0:
            args += " trs=%d" % self.originShift
        args += " refs=startcavgs.mrc"
        args += " oritab=prime2D_startdoc.txt"
        self.runJob("simple_prime2D", args, cwd=self._getExtraPath())

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
        return ['Elmlund2013']
    
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
        while os.path.exists(self.getDocFile(lastIter)):
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

    def _updateParticle(self, item, row):
        self._count += 1
        item.setClassId(int(float(row['class'])))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        item.setLocation((self._count, self.getParticlesStack()))

    def _updateClass(self, item):
        item.setAlignment2D()
        fn = self.getClassesStack(self._iter) + ":mrcs"
        item.getRepresentative().setLocation(item.getObjId(), fn)
