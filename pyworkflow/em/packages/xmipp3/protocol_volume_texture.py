# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
"""
Protocol to split a volume in two volumes based on a set of images
"""

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, StringParam
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume
from pyworkflow.em.convert import ImageHandler
from convert import writeSetOfParticles

class XmippProtVolumeTexture(ProtClassify3D):
    """Split volume in two"""
    _label = 'volume texture'
    
    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, label="Volume", important=True, 
                      pointerClass='Volume', help='Select a set of particles with angles. Preferrably the output of a run of directional classes')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      label="Reference", help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('mask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be binary: 0 (remove these voxels) and 1 (let them pass).')
      
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep',self.inputVolume.getObjId(),
                                                self.referenceVolume.getObjId())
        self._insertFunctionStep('findTexture')
        # self._insertFunctionStep('createOutput')

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputVolumeId, referenceVolumeId):
        inputVolume = self.inputVolume.get()
        referenceVolume = self.referenceVolume.get()

        img = ImageHandler()
        img.convert(inputVolume, self._getInputVolFn())

        imgR = ImageHandler()
        img.convert(referenceVolume, self._getRefVolFn())





    # def createOutput(self):
    #     inputParticles = self.directionalClasses.get()
    #     volumesSet = self._createSetOfVolumes()
    #     volumesSet.setSamplingRate(inputParticles.getSamplingRate())
    #     for i in range(2):
    #         vol = Volume()
    #         vol.setLocation(1, self._getExtraPath("split_v%d.vol"%(i+1)))
    #         volumesSet.append(vol)

    #     self._defineOutputs(outputVolumes=volumesSet)
    #     self._defineSourceRelation(inputParticles, volumesSet)
        
    #     vol = Volume()
    #     vol.setLocation(1, self._getExtraPath("split_pc1.vol"))

    #     volumesSet2 = self._createSetOfVolumes()
    #     volumesSet2.append(vol)

    #     self._defineOutputs(outputVolumes=volumesSet2)
    #     self._defineSourceRelation(inputParticles, volumesSet2)

    def findTexture(self):
        # inputParticles = self.directionalClasses.get()
        # Xdim = inputParticles.getDimensions()[0]
        fnMask = ""
        if self.mask.hasValue():
            fnMask = self._getExtraPath("mask.vol")
            img=ImageHandler()
            img.convert(self.mask.get(), fnMask)
            # self.runJob('xmipp_image_resize',"-i %s --dim %d"%(fnMask,Xdim),numberOfMpi=1)
            # self.runJob('xmipp_transform_threshold',"-i %s --select below 0.5 --substitute binarize"%fnMask,numberOfMpi=1)

        args="-i %s -r %s --patchSize 10"%\
             (self._getInputVolFn(),self._getRefVolFn())
        if fnMask!="":
            args+=" --mask %s"%fnMask
        self.runJob("xmipp_volume_texture",args)



    def _getInputVolFn(self):
        return self._getTmpPath('volume.vol')

    def _getRefVolFn(self):
        return self._getTmpPath('reference.vol')

    # def _getRefVolFn(self):
    #     return self._getTmpPath('reference.vol')