# **************************************************************************
# *
# * Authors:     Javier Vargas (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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

from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
import pyworkflow.protocol as pwprot
from pyworkflow.em.packages.xmipp3.convert import (getImageLocation,readSetOfVolumes,Volume)
import xmipp
from pyworkflow.utils import getExt
from numpy import zeros, double
from pyworkflow.object import Float, String


# import xmipp

        
class XmippProtVolumeOccupancy(ProtAnalysis3D):
    """Compare two states of a volume to analyze the local strains and rotations"""
    _label = 'calculate occupancy'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('refMap', PointerParam, label="Reference Map", important=True,
                      pointerClass='Volume',
                      help='Reference map used to calculate the occupancies. Usually comes from a PDB')
        
        form.addParam('inputMasks', pwprot.params.MultiPointerParam, 
                                            label="Input Masks", important=True,
                                            pointerClass='Volume', minNumObjects=2, maxNumObjects=1,
                                            help='Select two or more masks aligned with the reference map to calculate the occupancies.'
                                                 'inside this regions for the input maps')

        form.addParam('inputMaps', pwprot.params.MultiPointerParam, 
                                            label="Input Volumes", important=True,
                                            pointerClass='Volume', minNumObjects=2, maxNumObjects=1,
                                            help='Input maps in which we want to calculate the occupancies')
        
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
             
        #Read the reference map        
        self.refVol = ''
        self.mapVol = ''
        self.mask = ''
        
        self.numOfMaps =  sum(1 for _ in self.inputMaps)
        self.numOfMasks =  sum(1 for _ in self.inputMasks)
        
        self.volDensityInMask = zeros((self.numOfMasks,self.numOfMaps),dtype=double)
        self.refDensityInMask = zeros((self.numOfMasks,1),dtype=double)
        
        self._insertFunctionStep("readRefInMemory", self.refMap.get())
        numV = 0
        numMask = 0
        for itemVol in self.inputMaps:       # Take a look  to protocol_sets line 151 
            self._insertFunctionStep("readMapInMemory", itemVol.get())
            for itemMask in self.inputMasks:
                self._insertFunctionStep("readMaskInMemory", itemMask.get())                                            
                self._insertFunctionStep("compareVoxelDensities",numV,numMask)
                numMask += 1
            numMask = 0
            numV += 1
            
        self._insertFunctionStep('createOutputStep')

        
    #--------------------------- STEPS functions ---------------------------------------------------    
    def changeExtension(self, vol):
        extVol = getExt(vol)
        if (extVol == '.mrc') or (extVol == '.map'):
            vol = vol + ':mrc'
        return vol
    
    def readRefInMemory(self, mapLoc):
        refVolumeLoc = self.changeExtension(getImageLocation(mapLoc))
        self.refVol = xmipp.Image(refVolumeLoc)

    def readMapInMemory(self, mapLoc):
        refVolumeLoc = self.changeExtension(getImageLocation(mapLoc))
        self.mapVol = xmipp.Image(refVolumeLoc)
    
    def readMaskInMemory(self, maskLoc):
        refVolumeLoc = self.changeExtension(getImageLocation(maskLoc))
        self.mask = xmipp.Image(refVolumeLoc)          

    def compareVoxelDensities(self,numV,numMask):        
        dims = self.refVol.getDimensions()
        n = 0
        computeRef = (self.refDensityInMask[numMask,0] == 0)
        for i in range(0,dims[0]):
            for j in range(0,dims[1]):
                for k in range(0,dims[2]):
                    if (self.mask.getPixel(n,i,j,k) > 0.5):
                        self.volDensityInMask[numMask,numV] += self.mapVol.getPixel(n,i,j,k)
                        if (computeRef):
                            self.refDensityInMask[numMask,0] += self.refVol.getPixel(n,i,j,k)
                                
    def createOutputStep(self):
        fnOccupancy = self._getExtraPath('occupancy.xmd')
        md = xmipp.MetaData()        

        
        outputSetOfVolumes = self._createSetOfVolumes()        
        outputSetOfVolumes.setSamplingRate(self.inputMaps[0].get().getSamplingRate())    
            
        #outVol.setSamplingRate(self.refMap.get().getSamplingRate())        
        #outVol.setLocation(self.refMap.get().getFileName())
        #outVol.setObjLabel(self.refMap.get().getObjLabel())
                                
        for maskIdx in range(0,self.numOfMasks):
            outVol = Volume()
            outVol.copyAttributes(self.refMap.get())
            outVol.setFileName(self.refMap.get().getFileName())
            outVol.setObjComment("Reference Map")
            
            outVol.mask = String(self.inputMasks[maskIdx].get().getFileName())
            outVol.weight = Float(self.refDensityInMask[maskIdx,0])
            outputSetOfVolumes.append(outVol)                

            objId = md.addObject()
            md.setValue(xmipp.MDL_IMAGE,self.refMap.get().getFileName(),objId)                       
            md.setValue(xmipp.MDL_MASK,self.inputMasks[maskIdx].get().getFileName(),objId) 
            md.setValue(xmipp.MDL_WEIGHT,self.refDensityInMask[maskIdx,0],objId)
        
        for volIdx in range(0,self.numOfMaps):
            for maskIdx in range(0,self.numOfMasks):            
                
                outVol = Volume()
                outVol.copyAttributes(self.inputMaps[volIdx].get())
                outVol.setFileName(self.inputMaps[volIdx].get().getFileName())
                outVol.setObjComment("Input Map %s " % String(volIdx+1))
                outVol.mask = String(self.inputMasks[maskIdx].get().getFileName())
                outVol.weight = Float(self.volDensityInMask[maskIdx,volIdx])
                outputSetOfVolumes.append(outVol)                
                
                objId = md.addObject()
                md.setValue(xmipp.MDL_IMAGE,self.inputMaps[volIdx].get().getFileName(),objId)
                md.setValue(xmipp.MDL_MASK,self.inputMasks[maskIdx].get().getFileName(),objId) 
                md.setValue(xmipp.MDL_WEIGHT,self.volDensityInMask[maskIdx,volIdx],objId)
        
        md.write(fnOccupancy)
        outputArgs = {'outputVolumes': outputSetOfVolumes}
        self._defineOutputs(**outputArgs)
 

                                             
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        
        '''xdim0 = self.inputVolume0.get().getDim()[0]
        xdimF = self.inputVolumeF.get().getDim()[0]
        if xdim0 != xdimF:
            errors.append("Make sure that the two volumes have the same size")'''

        return errors    
