# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
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

from os.path import exists

from pyworkflow.object import Pointer
import pyworkflow.utils as pwutils
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.constants import (STEPS_PARALLEL, LEVEL_ADVANCED)
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtParticlePicking

from convert import writeSetOfCoordinates, readSetOfCoordinates
from xmipp3 import XmippProtocol


class XmippProtPickNoise(ProtParticlePicking, XmippProtocol):
    """Protocol to pick noise particles"""
    _label = 'pick noise'
    
    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label="Input coordinates",
                      help='Select the SetOfCoordinates ')

        form.addParam('extractNoiseNumber', params.IntParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label='Number of noise particles',
                      help='Number of noise particles to extract from each micrograph. '
                           'Set to -1 for extracting the same amount of noise particles as the number true particles for that micrograph')

        form.addParallelSection(threads=4, mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        """for each micrograph insert the steps to preprocess it"""
        # Write pos files for each micrograph
        firstStepId = self._insertFunctionStep('writePosFilesStep')
        
        # For each micrograph insert the steps, run in parallel
        deps = []
        
        boxSize = self.inputCoordinates.get().getBoxSize()
        for mic in self.getInputMicrographs():
            localDeps = [firstStepId]
 
            baseMicName = pwutils.removeBaseExt(mic.getFileName())
            deps.append(self._insertFunctionStep('pickNoiseStep', 
                                                 mic.getObjId(), baseMicName, mic.getFileName(), self.extractNoiseNumber.get(), boxSize,
                                                 prerequisites=localDeps))
# 
        self._insertFunctionStep('_createOutput', self._getExtraPath(),
                                 prerequisites=deps)

    #--------------------------- STEPS functions -------------------------------
    def writePosFilesStep(self):
        """ Write the pos file for each micrograph on metadata format. """
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoordinates.get(), scale=self.inputCoordinates.get().getBoxSize())
                        
    def pickNoiseStep(self, micId, baseMicName, micrographToExtract, extractNoiseNumber, boxSize):
        """ Pick noise from one micrograph """
        outputRoot = str(self._getExtraPath(baseMicName))
        fnPosFile = self._getExtraPath(baseMicName + ".pos")

        # If it has coordinates extract the particles      
        particlesMd = 'particles@%s' % fnPosFile

        if exists(fnPosFile):
            args = " -i %s --pos %s" % (micrographToExtract, particlesMd)
            args += " -o %s --Xdim %d" % (outputRoot, boxSize)
            args += " --extractNoise %d"%extractNoiseNumber

            self.runJob("xmipp_micrograph_scissor", args)
            cleanPath("%s.stk"%outputRoot)
            cleanPath("%s.xmd"%outputRoot)
        else:
            self.warning("The micrograph %s hasn't coordinate file! " % baseMicName)
            self.warning("Maybe you picked over a subset of micrographs")
    
    def getInputMicrographs(self):
        return self.inputCoordinates.get().getMicrographs()
    
    def getInputMicrographsPointer(self):
        ptr = Pointer()
        ptr.set(self.getInputMicrographs())
        return ptr
    
    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.getInputMicrographs(), coordSet)

    def _summary(self):
        summary = []
        if hasattr(self, 'outputCoordinates'):
            summary.append('%d noisy particles were picked'%self.outputCoordinates.getSize())
        return summary
