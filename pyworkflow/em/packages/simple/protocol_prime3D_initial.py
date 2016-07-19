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

        form.addParam('extraParams', params.StringParam,
                      label='Extra parameters',
                      help='')

        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputSet.getUniqueId())

        sym = self.getSym()
        initSym = 'c1' if self.searchSymAxis else sym
        self._insertFunctionStep('prime3DInitStep', initSym)
        # The following files should be output from the prime 3d - init step
        initVol, initOritab = 'startvol_state1.mrc', 'prime3D_startdoc.txt'
        self.finalRoot = self.getRoot(initSym)
        self._insertFunctionStep('prime3DStep', initSym, initVol, initOritab)

        if self.searchSymAxis:
            self._insertFunctionStep('findSymAxisStep')
            # The following files should be produced after finding symmetry axis
            vol, oritab = "recvol_state1.mrc", "sym_%s.txt" % sym
            self.finalRoot = self.getRoot(sym)
            self._insertFunctionStep('prime3DStep', sym, vol, oritab)

        if isinstance(self.inputSet.get(), em.writeSetOfClasses2D):
            self._insertFunctionStep('mapClassesToParticlesStep')

        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        inputSet = self.inputSet.get()
        stkFn = self._getExtraPath(self.getInputStack())

        if isinstance(inputSet, em.SetOfAverages):
            inputSet.writeStack(stkFn)
        elif isinstance(inputSet, em.SetOfClasses2D):
            firstClass = inputSet.getFirstItem()
            writeSetOfClasses2D(inputSet, stkFn,
                                stackFn=self._getExtraPath('particles.mrcs'),
                                docFn=self._getExtraPath('particles.txt'),
                                ctfFn=self.getCtfFile())
        elif isinstance(inputSet, em.SetOfParticles):
            writeSetOfParticles(inputSet, stkFn, docFn=None, ctfFn=None)
        else:
            raise Exception('Unexpected input type: %s' % type(inputSet))

    def prime3DInitStep(self, initSym):
        args = self._getCommonArgs()
        args += " pgrp=%s" % initSym

        self.runJob("simple_prime3D_init", args, cwd=self._getExtraPath())

    def prime3DStep(self, sym, vol, oritab):
        root = self._getExtraPath(self.getRoot(sym))
        pwutils.makePath(root)
        args = self._getCommonArgs(prefix='../')
        args += " vol1=../%s" % vol
        args += " oritab=../%s" % oritab
        args += " pgrp=%s" % sym
        args += " %s " % self.extraParams.get()

        self.runJob("simple_prime3D", args, cwd=root)

    def findSymAxisStep(self):
        root = self.getRoot('c1') # Get root from previous prime 3D
        docFile = self.getDocFile(root)
        volFile = self.getVolFile(root)
        symFile = 'sym_%s.txt' % self.getSym()
        doc = simple.SimpleDocFile(docFile)
        row = doc.getLastRow()
        doc.close()

        args = self._getCommonArgs(stk=False)
        args += " vol1=%s" % volFile.replace(self._getExtraPath(), '.')
        args += " oritab=%s" % docFile.replace(self._getExtraPath(), '.')
        args += " pgrp=%s" % self.getSym()
        args += " lp=%s" % row['lp']
        args += " outfile=%s" % symFile
        # Find symmetry axis
        self.runJob("simple_symsrch", args, cwd=self._getExtraPath())

        # Check some files that need to be produced
        if not os.path.exists(self._getExtraPath(symFile)):
            raise Exception('Expected file %s not produced.' % symFile)

        # Now reconstruct with the defined symmetry
        args = self._getCommonArgs()
        args += " oritab=%s" % symFile
        args += " pgrp=%s" % self.getSym()
        self.runJob("simple_eo_recvol", args, cwd=self._getExtraPath())

    def mapClassesToParticlesStep(self):
        # SIMPLE_MAP2PTCLS stk=<particles.ext> stk2=<selected_cavgs.ext>
        # stk3=<orig_cavgs.ext> oritab=<PRIME 2D doc> [oritab2=<prime3D shc doc>]
        # [comlindoc=<shc_clustering_nclsX.txt>] [doclist=<list of oritabs for the different states>]
        # [deftab=<text file defocus values>] [outfile=<output parameter file{mapped_ptcls_params.txt}>]
        # [nthr=<nr of OpenMP threads{1}>]
        root = self._getExtraPath('map2ptcls')
        pwutils.makePath(root)
        lastDoc = self.getDocFile(self.finalRoot).replace(self._getExtraPath(),
                                                          '..')
        mappedDoc = 'mapped_ptcls_params.txt'
        args = "stk=../particles.mrcs stk2=../averages.mrcs stk3=../averages.mrcs"
        args += " oritab=../particles.txt oritab2=%s" % lastDoc
        args += " outfile=%s" % mappedDoc
        ctfFn = self.getCtfFile()
        if ctfFn is not None:
            args += " deftab=%s" % ctfFn
        self.runJob("simple_map2ptcls", args, cwd=root)

        # Now reconstruct mapped particles
        args = self._getCommonArgs(stk='particles.mrcs', prefix='../')
        args += " oritab=%s" % mappedDoc
        args += " pgrp=%s" % self.getSym()
        self.runJob("simple_eo_recvol", args, cwd=root)

    def _getCommonArgs(self, stk=None, prefix=''):
        inputSet = self.inputSet.get()

        if stk is None:
            stk = self.getInputStack()

        args = " stk=%s%s" % (prefix, stk) if stk else ""
        args += " smpd=%f" % inputSet.getSamplingRate()
        args += " msk=%d" % self.getMaskRadius()
        args += " nthr=%d" % self.numberOfThreads

        return args

    def createOutputStep(self):
        vol = em.Volume(self.getVolFile(self.finalRoot))
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

    def getRoot(self, sym):
        return 'prime3D_%s' % sym

    def getInputStack(self):
        if isinstance(self.inputSet.get(), em.SetOfParticles):
            return "particles.mrcs"
        else:
            return "averages.mrcs"

    def getCtfFile(self):
        firstClass = self.inputSet.get().getFirstItem()
        return 'particles_ctf.txt' if firstClass.hasCTF() else None

    def getSym(self):
        return self.symmetry.get().lower()

    def getDocFile(self, root, iteration=None):
        """ Return the document file with alignment parameters for the given
        iteration, if None passed, return the last iteration.
        """
        if iteration is None:
            iteration = self.getLastIteration(root)

        return self._getExtraPath(root, "prime3D_doc%d.txt" % iteration)

    def getVolFile(self, root, iteration=None):
        """ Return the document file with alignment parameters for the given
        iteration, if None passed, return the last iteration.
        """
        if iteration is None:
            iteration = self.getLastIteration(root)

        return self._getExtraPath(root, "recvol_state1_iter%d.mrc" % iteration)

    def getLastIteration(self, root):
        lastIter = 1
        while os.path.exists(self.getVolFile(root, lastIter)):
            lastIter += 1
        return lastIter - 1