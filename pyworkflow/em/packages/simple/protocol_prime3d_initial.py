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
import convert


class ProtPrime3DInitial(em.ProtInitialVolume):
    """
    This protocol wraps *ini3D_from_cavgs*, to generate an initial
    volume from 2D averages.

    josem@artemis:~/work/development/SIMPLE3.0$ simple_distr_exec prg=ini3D_from_cavgs
USAGE:
bash-3.2$ simple_exec prg=simple_program key1=val1 key2=val2 ...

REQUIRED
stk    = particle stack with all images(ptcls.ext)
smpd   = sampling distance, same as EMANs apix(in A)
msk    = mask radius(in pixels)
pgrp   = point-group symmetry(cn|dn|t|o|i)
nparts = # partitions in distributed exection

OPTIONAL
nthr        = # OpenMP threads{1}
nthr_master = # OpenMP threads on master node{1}
ncunits     = # computing units, can be < nparts {nparts}
hp          = high-pass limit(in A)
lp          = low-pass limit(in A)
frac        = fraction of ptcls(0-1){1}
automsk     = envelope masking(yes|no|cavg){no}
mw          = molecular weight(in kD)
amsklp      = low-pass limit for envelope mask generation(in A)
edge        = edge size for softening molecular envelope(in pixels)
binwidth    = binary layers grown for molecular envelope(in pixels){1}
inner       = inner mask radius(in pixels)
width       = falloff of inner mask(in pixels){10}
nspace      = # projection directions
shbarrier   = use shift search barrier constraint(yes|no){yes}
    """


    _label = 'prime 3d initial'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfClasses,SetOfAverages',
                      label="Input references", important=True,
                      help='You can selected the following type of sets:'
                           'SetOfClasses or SetOfAverages '
                           'in order to produce an initial volume. ')

        form.addParam('maskRadius', params.IntParam,
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

        form.addParallelSection(threads=4, mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',
                                 self.inputSet.getUniqueId())
        self._insertFunctionStep('createInitialVolumeStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, inputId):
        inputSet = self.inputSet.get()
        inputSet.writeStack(self._getExtraPath('input_references.mrcs'))

    def createInitialVolumeStep(self):
        inputSet = self.inputSet.get()
        args = " stk=input_references.mrcs"
        args += " smpd=%0.3f" % inputSet.getSamplingRate()
        args += " msk=%d" % self.maskRadius
        args += " pgrp=%s" % self.symmetry
        args += " %s " % self.extraParams.get('')
        args += " nparts=%d" % self.numberOfMpi
        args += " nthr=%d" % self.numberOfThreads

        self.runJob(simple.getProgram('ini3D_from_cavgs', distr=True),
                    args, cwd=self._getExtraPath())

    def createOutputStep(self, v=1):
        return
        inputSet = self.inputSet.get()

        if isinstance(inputSet, em.SetOfClasses2D):
            root = self._getExtraPath('map2ptcls')
            # We should take the volume from the output of map2ptcls
            volFile = os.path.join(root, 'recvol_state1.mrc')
            docFile = os.path.join(root, 'mapped_ptcls_params.txt')
            # Create an output particles with new angular assignment
            partSet = self._createSetOfParticles()
            partSet.copyInfo(inputSet.getImages())
            partSet.setAlignmentProj()
            convert.particlesFromClasses(inputSet, partSet, docFile)
            self._defineOutputs(outputParticles=partSet)
            self._defineSourceRelation(self.inputSet, partSet)

        else:
            volFile = self.getVolFile(self.finalRoot)

        vol = em.Volume(volFile)
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