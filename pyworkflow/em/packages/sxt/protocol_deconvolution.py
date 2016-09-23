# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Joaquin Oton   (joton@cnb.csic.es)
# *              Marc Rosanes   (mrosanes@cells.es)
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

import pyworkflow.protocol.params as params
from pyworkflow.em import Protocol
#from os.path import basename
#from pyworkflow.utils.path import removeExt
#from h5py import File
#import xmipp
#from pyworkflow.utils import getFloatListFromValues
#import numpy as np
#import pyworkflow.em as em


class ProtDeconvolution(Protocol):

    """    
    This protocol is aimed to deconvolve tiltSeries and 3D PSF.        
    """
    
    _label = 'deconvolving tiltSeries and 3d psf'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTiltSeries', params.PointerParam, 
                      pointerClass='TiltSeries', 
                      label="------",  
                      help="----")
        form.addParam('psf', params.PointerParam, 
                      pointerClass='?????',
                      label="---",
                      help="----")
        form.addParam('kw', params.FloatParam, 
                      label='---',
                      help="----")
        
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):         
        
               
        #self._insertFunctionStep('validateSiemensStar', inputSS)   
        #self._insertFunctionStep('getPsfFromSiemensStar', inputSS, nRef, orders, ringPos, fnOutPsf)
        #self._insertFunctionStep('createOutputStep', fnOutPsf)            
    #--------------------------- STEPS functions --------------------------------------------
    
    #def validateSiemensStar(self, inputSS):
         
    
    #def getPsfFromSiemensStar(self, inputSS, nRef, orders, ringPos, fnOutPsf):        
        
    
    
    #def createOutputStep(self, fnOutPsf): 
        
        
    #--------------------------- INFO functions -------------------------------------------- 
    
#    def _summary(self):
#        summary = []
#        summary.append[]        
#        summary.append[]
#        return summary
#    
#    def _methods(self):
#        messages = []
#        messages.append('Joton')
#        return messages

#    def _validate(self):
#        errors = []
#        hdf5 = self.inputSiemensStar.get().endswith('hdf5')
#        if self.inputSiemensStar.get():
#            if not hdf5:
#                errors.append ("Expected hdf5 files for importing!!!") 
#        else:
#            errors.append("The path can not be empty!!!")      
#        return errors              
    #--------------------------- UTILS functions --------------------------------------------
    
#    def _defineOutputName(self):
#        return self._getExtraPath('psf.mrc')