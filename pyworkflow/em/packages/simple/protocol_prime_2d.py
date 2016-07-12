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
import pyworkflow.em as em

import simple



class ProtPrime2D(em.ProtClassify2D):
    _label = 'prime 2d'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputParticles', params.PointerParam,
                      #pointerCondition='hasCTF', # or just a warning?
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='')

        form.addParam('numberOfClasses', params.IntParam,
                      label="Number of classes",
                      help="")

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runPrime2D')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        self.inputParticles.get().writeStack(self._getExtraPath("particles.mrcs"))
            
    def runPrime2D(self):
        # simple_prime2D_init stk=<stack.ext> smpd=<sampling distance(in A)>
        # msk=<mask radius(in pixels)> ncls=<nr of clusters>
        # [nthr=<nr of OpenMP threads{1}>] [oritab=<input doc>]

        inputParticles = self.inputParticles.get()
        xdim, _, _ = inputParticles.getDimensions()
        # We will run simple_prime2d in the extra folder, so 'particles.mrcs'
        # should be there
        args = " stk=particles.mrcs"
        args += " smpd=%f" % inputParticles.getSamplingRate()
        args += " msk=30" #FIXME
        args += " ncls=%d" % self.numberOfClasses
        args += " nthr=%d" % self.numberOfThreads

        self.runJob("simple_prime2D_init", args, cwd=self._getExtraPath())

        # simple_prime2D stk=<stack.ext> smpd=<sampling distance(in A)>
        # msk=<mask radius(in pixels)> ncls=<nr of clusters>
        # refs=<initial_references.ext> oritab=<previous clustering doc>
        # [lp=<low-pass limit(in A){20}>] [trs=<origin shift(in pixels){0}>]
        # [nthr=<nr of OpenMP threads{1}>] [startit=<start iteration>]
        # [hp=<high-pass limit(in A)>] [srch_inpl=<yes|no{yes}>]
        # ** less commonly used**
        # [maxits=<max iterations{500}>] [inner=<inner mask radius(in pixels)>]
        # [width=<pixels falloff inner mask(in pixels){10}>]

        # Reusing same arguments from simple_prime2D_init
        args += " refs=startcavgs.mrc"
        args += " oritab=prime2D_startdoc.txt"
        self.runJob("simple_prime2D", args, cwd=self._getExtraPath())


    def createOutputStep(self):
        lastIter = self.getLastIteration()
        docFile = self._getExtraPath('prime2D_doc%d.txt' % lastIter)
        doc = simple.SimpleDocFile(docFile)

        for row in doc:
            print row

        doc.close()

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

    # --------------------------- INFO functions -------------------------------
    def getLastIteration(self):
        lastIter = 1
        pattern = self._getExtraPath("prime2D_doc%d.txt")
        while os.path.exists(pattern % lastIter):
            lastIter += 1
        return lastIter - 1