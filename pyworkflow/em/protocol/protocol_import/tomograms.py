# **************************************************************************
# *
# * Author:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol base classes related to EM imports of Tomograms
"""


from pyworkflow.protocol import params
import pyworkflow.utils.path as putils

from images import ProtImportFiles
from micrographs import ProtImportMicBase
from pyworkflow.em import TomoRec

class ProtImportTomograms(ProtImportMicBase):
    """Protocol to import a set of tomograms to the project"""
    _label = 'import tomograms'
    _outputClassName = 'SetOfTomograms'
#     _checkStacks = False    
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportMicBase._defineParams(self, form)    
        
        form.addParam('bfactor', params.FloatParam, default=4.0,
                      label='Provide B-factor:',
                      help= '3D CTF model weighting B-factor per e-/A2')
        form.addParam('totalDose', params.FloatParam, default=40.0,
                      label='Provide acummulated dose:',
                      help= 'Total dose for the whole tomogram.')
    
    def setSamplingRate(self, tomoSet):
        ProtImportMicBase.setSamplingRate(self, tomoSet)
        tomoSet.setBfactor(self.bfactor.get())
        tomoSet.setDose(self.totalDose.get())


class ProtImportTomoRecss(ProtImportFiles):
    """Protocol to import a set of reconstructed tomograms to the project"""
    _label = 'import reconstructed tomograms'
    _outputClassName = 'SetOfTomoRecs'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportFiles._defineParams(self, form)    
        
        form.addParam('inputTomograms', params.PointerParam, pointerClass='SetOfTomograms', 
                          label='Input tomograms',
                          help='Select the subtomograms that you want to import coordinates.')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep', self.filesPath.get())

    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self, importFrom, *args):
        tomoSet = self.inputTomograms.get()
        tomoRecsSet = self._createSetOfTomoRecs()
        tomoRecsSet.copyInfo(tomoSet)
        tomoRecsSet.setTomograms(tomoSet)
        
        for tomoRecFile, _ in self.iterFiles():
            tomo = self.getMatchingTomo(tomoRecFile)
            if tomo is not None:
                tomoRec = TomoRec()
                dest = self._getExtraPath(putils.replaceBaseExt(tomoRecFile, "mrc"))

                copyOrLink = self.getCopyOrLink()
                copyOrLink(tomoRecFile, dest)
                
                tomoRec.setFileName(dest)
                tomoRec.copyInfo(tomo)
                tomoRec.setTomogram(tomo)
                tomoRecsSet.append(tomoRec)
        
        self._defineOutputs(outputTomoRecs=tomoRecsSet)
        self._defineSourceRelation(tomoSet, tomoRecsSet)
        
    def getMatchingTomo(self, tomoRecFile):
        tomoRecBase = putils.removeBaseExt(tomoRecFile)
        for tomo in self.inputTomograms.get():
            tomoBase = putils.removeBaseExt(tomo.getFileName())
            if tomoRecBase in tomoBase or tomoBase in tomoRecBase: #temporal use of in
                return tomo
        return None

    def getCopyOrLink(self):    
        # Set a function to copyFile or createLink
        # depending in the user selected option 
        if self.copyFiles:
            return putils.copyFile
        else:
            return putils.createLink
