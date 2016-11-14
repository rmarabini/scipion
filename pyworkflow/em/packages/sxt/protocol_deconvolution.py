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
from pyworkflow.em.packages.sxt.data import TiltSeries, SetOfTiltSeries
import numpy as np
import pyworkflow.em as em
from pyworkflow.utils.path import cleanPath
from pyworkflow.mapper.sqlite_db import SqliteDb

class ProtDeconvolution(Protocol):

    """    
    This protocol is aimed to deconvolve tiltSeries and 3D PSF.        
    """
    
    _label = 'deconvolving tiltSeries and 3d psf'
    
    TILT_SERIES = 0
    SET_OF_TILT_SERIES = 1
    FOCAL_SERIES = 2
      
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', params.EnumParam, 
                            choices=['Tilt series', 'Set of tilt series', 'Focal series'], 
                            default=self.TILT_SERIES,
                            label='Type of input tilt series',
                            help='Select the type of input tilt series.')
        form.addParam('inputTiltSeries', params.PointerParam,
                      condition = '(inputType == %d)' % self.TILT_SERIES,
                      pointerClass='TiltSeries', 
                      label="Input tilt series",  
                      help="Select input tilt series from the project")
        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      condition = '(inputType == %d)' % self.SET_OF_TILT_SERIES,
                      pointerClass='SetOfTiltSeries', 
                      label="Input set of tilt series",  
                      help="Select input set of tilt series from the project")
        form.addParam('inputFocalSeries', params.PointerParam,
                      condition = '(inputType == %d)' % self.FOCAL_SERIES,
                      pointerClass='SetOfTiltSeries', 
                      label="Input focal series",  
                      help="Select input focal series from the project")        
        form.addParam('inputPsf', params.PointerParam, 
                      pointerClass='PSF3D', ########################################333new obj ... do necessary changes in the code ... PSF3D
                      label="Input 3D PSF",
                      help="Select input PSF based on the imaging date")
        form.addParam('kw', params.FloatParam,default = 0.05, 
                      label='k-Factor',
                      help="----")
        
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        inputType = self.inputType.get()
        if inputType == self.TILT_SERIES:
            self._insertFunctionStep('createOutputStepTiltSeries')
        elif inputType == self.SET_OF_TILT_SERIES:
            self._insertFunctionStep('createOutputStepSetOfTiltSeries')
        else:
            self._insertFunctionStep('createOutputStepFocalSeries')            
    #--------------------------- STEPS functions --------------------------------------------
        
    def createOutputStepTiltSeries(self):
        inputTiltSeries = self.inputTiltSeries.get()
        inputPSF = self.inputPsf.get()
        psfPixelSizeX = self.inputPsf.get().getSamplingRate()/10 ######################################## ba obj jadid baray PSF3D bayad modify shavad
        psfArrayToDeconv = self._preDeconvolution(inputPSF)
        kw = self.kw.get()
                
        ih = em.ImageHandler()
        inputTiltSeriesArray = ih.read(inputTiltSeries).getData()
        tiltSeriesPixelSize = inputTiltSeries.getSamplingRate()/10        
        
        from xpytools.mtf_deconv_wiener import MTFDeconvWiener
        deconvolutionObj = MTFDeconvWiener()
        deconvTiltSeriesArray = deconvolutionObj.mtf_deconv_wiener(
                            inputTiltSeriesArray, tiltSeriesPixelSize, 
                            psfArrayToDeconv, psfPixelSizeX, kw, pad=20, fc=-1)
        
        outputTilt = ih.createImage()        
        fnOutTiltSeries = self._defineOutputName(1)
        k = 0
        for j in range(np.shape(deconvTiltSeriesArray)[0]):
            outputTilt.setData(deconvTiltSeriesArray[j, :, :])
            k += 1
            outputTilt.write((k,fnOutTiltSeries))
        
        tiltSeries = TiltSeries()
        tiltSeries.setFileName(fnOutTiltSeries)
        acquisition = tiltSeries.getXrayAcquisition()        
        acquisition.setLensLabel(inputTiltSeries.getXrayAcquisition().getLensLabel())
        acquisition.setEnergy(inputTiltSeries.getXrayAcquisition().getEnergy())
        acquisition.setDate(inputTiltSeries.getXrayAcquisition().getDate())
        tiltSeries.setXrayAcquisition(acquisition)
        tiltSeries.setSamplingRate(inputTiltSeries.getSamplingRate())
        angles = inputTiltSeries.getAngles()        
        tiltSeries.setAngles(angles)        
        tiltSeries.setSize(inputTiltSeries.getSize())
        self._defineOutputs(outputTiltSeries=tiltSeries)            
            
    def createOutputStepSetOfTiltSeries(self):
    
        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()
        inputPSF = self.inputPsf.get()
        psfPixelSizeX = self.inputPsf.get().getSamplingRate()/10 ######################################## ba obj jadid baray PSF3D bayad modify shavad
        psfArrayToDeconv = self._preDeconvolution(inputPSF)
        kw = self.kw.get()
        
        
        from xpytools.mtf_deconv_wiener import MTFDeconvWiener
        deconvolutionObj = MTFDeconvWiener()
        ih = em.ImageHandler()
        outputSetOfTilt = self._createSetOfTiltSeries()
        
        for i, tilt in enumerate (inputSetOfTiltSeries.iterItems()):
            inputTiltSeriesArray = ih.read(tilt).getData()
            tiltSeriesPixelSize = tilt.getSamplingRate()/10
            deconvTiltSeriesArray = deconvolutionObj.mtf_deconv_wiener(
                                inputTiltSeriesArray, tiltSeriesPixelSize, 
                                psfArrayToDeconv, psfPixelSizeX, kw, pad=20, fc=-1)
        
            outputTiltSeries = ih.createImage()        
            fnOutTiltSeries = self._defineOutputName(i)
            k = 0
            for j in range(np.shape(deconvTiltSeriesArray)[0]):
                outputTiltSeries.setData(deconvTiltSeriesArray[j, :, :])
                k += 1
                outputTiltSeries.write((k,fnOutTiltSeries))
        
            tiltSeries = TiltSeries()
            tiltSeries.setFileName(fnOutTiltSeries)
            acquisition = tiltSeries.getXrayAcquisition()        
            acquisition.setLensLabel(tilt.getXrayAcquisition().getLensLabel())
            acquisition.setEnergy(tilt.getXrayAcquisition().getEnergy())
            acquisition.setDate(tilt.getXrayAcquisition().getDate())
            tiltSeries.setXrayAcquisition(acquisition)
            tiltSeries.setSamplingRate(tilt.getSamplingRate())
            angles = tilt.getAngles()        
            tiltSeries.setAngles(angles)        
            tiltSeries.setSize(tilt.getSize())
            outputSetOfTilt.append(tiltSeries)           
            
        self._defineOutputs(outputSetOfTiltSeries=outputSetOfTilt)
    
    def createOutputStepFocalSeries(self): 
        
        
        x=1     
    #--------------------------- INFO functions -------------------------------------------- 
    
    #def _summary(self):
    #    summary = []
    #    outputSet = self._getOutputSetTiltSeries()
    #    if outputSet is None:
    #        summary.append("Output TiltSeries is not ready yet. Deconvolution "
    #                       "process is in progress.") 
    #    if  outputSet is not None:  
    #        summary.append("*%d* "% outputSet.getSize() +
    #                       "images related to the input TiltSeries "
    #                       "deconvolved with input 3D PSF.")             
    #        summary.append("Sampling rate : *%0.2f* A/px" % (
    #                        outputSet.getSamplingRate()))
    #        summary.append("*Imaging info*:\n" +
    #                       "Lens label: %s \n" % (
    #                        outputSet.getXrayAcquisition().getLensLabel()) +
    #                       "Energy (ev): %f \n" % (
    #                        outputSet.getXrayAcquisition().getEnergy()) +
    #                       "Date of imaging (ddmmyyy): %s" % (
    #                        outputSet.getXrayAcquisition().getDate()))             
    #    return summary
    
    #def _methods(self):
    #    methods = []
    #    outputSet = self._getOutputSetTiltSeries()
    #    if outputSet is not None:
    #        methods.append("*%d* %s were deconvolved with the input 3D PSF. "
    #                       " Output set is %s with a sampling rate of *%0.2f* A/px, "
    #                       "Lens label: %s, Energy (ev): %f, and "
    #                       "Date of imaging (ddmmyyy): %s)."                           
    #                       % (outputSet.getSize(), 
    #                          'images related to the input TiltSeries',
    #                          self.getObjectTag('outputTiltSeries'),
    #                          outputSet.getSamplingRate(),
    #                          outputSet.getXrayAcquisition().getLensLabel(),
    #                          outputSet.getXrayAcquisition().getEnergy(),
    #                          outputSet.getXrayAcquisition().getDate()))          
    #    return methods 
    #--------------------------- UTILS functions --------------------------------------------
    
    def _defineOutputName(self, suffix):
        return self._getExtraPath('deconvolved-tiltSeries_%01d.stk' % suffix)
    
    #def _getOutputSetTiltSeries(self):
    #    return getattr(self, 'outputTiltSeries', None)
    
    def _preDeconvolution(self, inputPSF):
        ih = em.ImageHandler()
        inputPsfArray = ih.read(inputPSF).getData()        
        psfPixelSizeZ = self.inputPsf.get().getZpixelSize() ######################################## ba obj jadid baray PSF3D bayad modify shavad
        #psfPixelSizeZ = 150
        #claculating best psf image to get their mean and use in deconvolution process
        xCenter = np.floor(np.shape(inputPsfArray)[1]/2)
        centerValues = inputPsfArray[:, xCenter, xCenter]
        thr = 0.3
        thrv = max(centerValues) * thr
        #finding indices based on centerValues curne and thrv      
        #mp: indices array above threshold, zm: central index
        #zN: indices range that will use to select psf images to get their average
        mp = np.where(centerValues > thrv)
        zm = mp[0][np.floor((mp[0][0] - mp[0][-1])/2)]
        zN = np.floor(np.shape(inputPsfArray)[0]/psfPixelSizeZ/2);
        zPos = []
        for i in range(2*int(zN)+1):
            zPos.append(int(zm-zN+i))
        print "PSF image indices to calculate their mean are: ", zPos        
        psfArrayToDeconv = np.mean(inputPsfArray[zPos, :, :], axis=0)
        return psfArrayToDeconv
    
    def _createSetOfTiltSeries(self):
        """ Create a set and set the filename. 
        If the file exists, it will be delete. """
        setFn = self._getPath('focalSeries.sqlite')
        # Close the connection to the database if
        # it is open before deleting the file
        cleanPath(setFn)        
        SqliteDb.closeConnection(setFn)        
        setObj = SetOfTiltSeries(filename=setFn)
        return setObj
    
    