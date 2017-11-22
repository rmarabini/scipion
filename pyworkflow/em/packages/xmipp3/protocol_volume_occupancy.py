# **************************************************************************
# *
# * Authors:     Javier Vargas (javier.vargasbalbuena@mcgill.ca)
# *
# * Departament of Anatomy and Cell Biology, McGill University
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

from pyworkflow.object import Float, String
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam
from pyworkflow.em.protocol.protocol import EMProtocol
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.utils import cleanPath
from xmipp3 import getMatlabEnviron

import pyworkflow.protocol as pwprot
from pyworkflow.em.packages.xmipp3.convert import (getImageLocation)
import xmipp
from pyworkflow.utils import getExt


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
                                            pointerClass='SetOfVolumes,VolumeMask', minNumObjects=2, maxNumObjects=1,
                                            help='Select two or more masks aligned with the reference map to calculate the occupancies.'
                                                 'inside this regions for the input maps')

        form.addParam('inputMaps', pwprot.params.MultiPointerParam, 
                                            label="Input Volumes", important=True,
                                            pointerClass='SetOfVolumes,Volume', minNumObjects=2, maxNumObjects=1,
                                            help='Input maps in which we want to calculate the occupancies')
        
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
             
        self._insertFunctionStep("readDataInMemory")
        self._insertFunctionStep("compareVoxelDensities")
        
    #--------------------------- STEPS functions ---------------------------------------------------
    
    def changeExtension(self, vol):
        extVol = getExt(vol)
        if (extVol == '.mrc') or (extVol == '.map'):
            vol = vol + ':mrc'
        return vol
    
    def readDataInMemory(self):
        refVolumeLoc = self.changeExtension(getImageLocation(self.refMap.get()))
        print refVolumeLoc
        
        #print self.inputMasks[1].get()
        #print self.inputMasks[0].get()
        #print len(self.inputMasks)


        self.reVol =  xmipp.Image(refVolumeLoc)
        print(self.reVol)

    def compareVoxelDensities(self):
        dims = self.reVol.getDimensions()
        print dims
        n = 0
        for i in range(0,dims[0]):
            for j in range(0,dims[1]):
                for k in range(0,dims[2]):
                    pass
                    #a = self.reVol.getPixel(n,i,j,k)
                         
        
        
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        
        '''xdim0 = self.inputVolume0.get().getDim()[0]
        xdimF = self.inputVolumeF.get().getDim()[0]
        if xdim0 != xdimF:
            errors.append("Make sure that the two volumes have the same size")'''

        return errors    
