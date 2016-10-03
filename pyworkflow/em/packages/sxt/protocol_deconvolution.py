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
from pyworkflow.em.packages.xmipp3.convert import writeSetOfImages, imageToRow
#from em.packages.xmipp3.convert import imageToRow

#from os.path import basename
#from pyworkflow.utils.path import removeExt
#from h5py import File
#import xmipp
#from pyworkflow.utils import getFloatListFromValues
import numpy as np
import pyworkflow.em as em


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
                      label="Input tilt series",  
                      help="Select input tilt series from the project")
        form.addParam('inputPsf', params.PointerParam, 
                      pointerClass='Volume',
                      label="Input 3D PSF",
                      help="Select input PSF based on the imaging date")
        form.addParam('kw', params.FloatParam, 
                      label='Inverse SNR',
                      help="----")
        
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):         
        
               
        self._insertFunctionStep('deconvolutionStep')   
        #self._insertFunctionStep('getPsfFromSiemensStar', inputSS, nRef, orders, ringPos, fnOutPsf)
        #self._insertFunctionStep('createOutputStep', fnOutPsf)            
    #--------------------------- STEPS functions --------------------------------------------
    
    def deconvolutionStep(self):
         
        
        inputTiltSeries = self.inputTiltSeries.get()
        inputPSF = self.inputPsf.get()
        
        ih = em.ImageHandler()
        inputTiltSeriesArray = ih.read(inputTiltSeries).getData()
        inputPsfArray = ih.read(inputPSF).getData()
        
        
        print "inputTiltSeriesArray dimension = ", np.shape(inputTiltSeriesArray)
        print "inputPsfArray dimension = ", np.shape(inputPsfArray)
        
        #numrows = infoshape[0][1]
        #numcols = infoshape[0][2]
        #offset = [num_img, 0, 0]
        #img_pixels = [1, numrows, numcols]
        #img_slab = input_nexusfile.getslab(offset, img_pixels)
        #img = np.zeros((numrows, numcols))
        #img[:,:] = img_slab[0,:,:]

        tiltSeriesPixelSize = self.inputTiltSeries.get().getSamplingRate() / 10
        psfPixelSizeX = self.inputPsf.get().getSamplingRate() / 10
        #psfPixelSizeZ = self.inputPsf.get().getZpixelSize()
        psfPixelSizeZ = 150
        
        print "PSF pixelSizeX(nm) = ", psfPixelSizeX
        print "PSF pixelSizeZ(nm) = ", psfPixelSizeZ
        print "TiltSeries pixelSize(nm) = ", tiltSeriesPixelSize
        
        xCenter = np.floor(np.shape(inputPsfArray)[1]/2)
        centerValues = inputPsfArray[:, xCenter, xCenter]
        thr = 0.3
        thrv = max(centerValues) * thr
        
        
        mp = np.where(centerValues > thrv);
        
        zm = mp[0][np.floor((mp[0][0] - mp[0][-1])/2)]
        
        
        
        zN = np.floor(np.shape(inputPsfArray)[0]/psfPixelSizeZ/2);
        
        zPos = []
        for i in range(2*int(zN)+1):
            zPos.append(int(zm-zN+i))
        
        print "PSF image indices to calculate their mean are: ", zPos
        
        
        
        psfArrayToDeconv = np.mean(inputPsfArray[zPos, :, :], axis=0)
        
        

        
        
        
        kw = self.kw.get()
        
        from xpytools.mtf_deconv_wiener import MTFDeconvWiener
        deconvolutionObj = MTFDeconvWiener()
        deconvTiltSeriesArray = deconvolutionObj.mtf_deconv_wiener(inputTiltSeriesArray, tiltSeriesPixelSize, psfArrayToDeconv, psfPixelSizeX, kw, pad=20, fc=-1)
        
        
        print "deconvTiltSeriesArray dimension = ", np.shape(deconvTiltSeriesArray)
        
        
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