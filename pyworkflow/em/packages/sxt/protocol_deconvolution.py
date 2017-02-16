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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys
import pyworkflow.protocol.params as params
from pyworkflow.em.packages.sxt.protocol_import import ProtImportTiltSeries
from pyworkflow.em.packages.sxt.data import TiltSeries, SetOfTiltSeries
import numpy as np
import pyworkflow.em as em
from pyworkflow.utils.path import cleanPath, removeExt
from pyworkflow.mapper.sqlite_db import SqliteDb
import xmipp

class ProtDeconvolution(ProtImportTiltSeries):

    """    
    This protocol is aimed to deconvolve tiltSeries/setOfTiltSeies/focalSeries and 3D PSF.        
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
        form.addParam('firstAngle', params.FloatParam,default = -65.0021,
                      condition = '(inputType == %d)' % self.FOCAL_SERIES, 
                      label='First angle',
                      help="The angle of the first image in degrees")
        form.addParam('incrementStep', params.FloatParam,default = 1.004, 
                      condition = '(inputType == %d)' % self.FOCAL_SERIES,
                      label='Angles increment step',
                      help="Angles increment step in degrees")
        
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        inputType = self.inputType.get()
        inputPSF = self.inputPsf.get()
        if inputType == self.TILT_SERIES:
            self._insertFunctionStep('createOutputStepTiltSeries', inputPSF)
        elif inputType == self.SET_OF_TILT_SERIES:
            self._insertFunctionStep('createOutputStepSetOfTiltSeries', inputPSF)
        else:
            self._insertFunctionStep('focalSeriesStdDeconvolution', inputPSF)
            #self._insertFunctionStep('createOutputStepFocalSeries')            
    #--------------------------- STEPS functions --------------------------------------------
        
    def createOutputStepTiltSeries(self, inputPSF):
        inputTiltSeries = self.inputTiltSeries.get()
        
        psfPixelSizeX = self.inputPsf.get().getSamplingRate()/10 ##if any new object like XpixelSize will be defined, it maybe replaces ## divided by 10 is to go through A and nm
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
        self._createTiltSeriesMd(1, fnOutTiltSeries, angles)
        self._defineOutputs(outputTiltSeries=tiltSeries)            
            
    def createOutputStepSetOfTiltSeries(self, inputPSF):    
        inputSetOfTiltSeries = self.inputSetOfTiltSeries.get()
        
        psfPixelSizeX = self.inputPsf.get().getSamplingRate()/10 ##if any new object like XpixelSize will be defined, it maybe replaces ## divided by 10 is to go through A and nm
        psfArrayToDeconv = self._preDeconvolution(inputPSF)
        kw = self.kw.get()
        
        
        from xpytools.mtf_deconv_wiener import MTFDeconvWiener
        deconvolutionObj = MTFDeconvWiener()
        ih = em.ImageHandler()
        outputSetOfTilt = self._createSetOfTiltSeries('SetOfTilt')
        outputSetOfTilt.setSamplingRate(inputSetOfTiltSeries.getSamplingRate())
        
        mdOut = xmipp.MetaData()
        for i, tilt in enumerate (inputSetOfTiltSeries.iterItems()):            
            inputTiltSeriesArray = ih.read(tilt).getData()
            tiltSeriesPixelSize = tilt.getSamplingRate()/10
            deconvTiltSeriesArray = deconvolutionObj.mtf_deconv_wiener(
                                inputTiltSeriesArray, tiltSeriesPixelSize, 
                                psfArrayToDeconv, psfPixelSizeX, kw, pad=20, fc=-1)
        
            outputTiltSeries = ih.createImage()        
            fnOutTiltSeries = self._defineOutputName(i+1)
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
            fnOutMd = self._createTiltSeriesMd(i+1, fnOutTiltSeries, angles)
            objId = mdOut.addObject()
            mdOut.setValue(xmipp.MDL_TOMOGRAMMD, fnOutMd, objId)
            outputSetOfTilt.append(tiltSeries)
            sys.stdout.write("\rinputTiltSeries number %d Deconvolved \n\n" % (i+1))
            sys.stdout.flush()           
        mdOut.write(self._getExtraPath('deconvolvedSetOfTiltSeries.xmd'))
        self._defineOutputs(outputSetOfTiltSeries=outputSetOfTilt)
    
    def focalSeriesStdDeconvolution(self, inputPSF):
        inputFocalSeries = self.inputFocalSeries.get()
        
        psfPixelSizeX = self.inputPsf.get().getSamplingRate()/10 ##if any new object like XpixelSize will be defined, it maybe replaces ## divided by 10 is to go through A and nm
        psfArrayToDeconv = self._preDeconvolution(inputPSF)
        kw = self.kw.get()
        
        
        from xpytools.mtf_deconv_wiener import MTFDeconvWiener
        deconvolutionObj = MTFDeconvWiener()
        ih = em.ImageHandler()
        
        ### for test
        for i, tilt in enumerate (inputFocalSeries.iterItems()):            
            inputTiltSeriesArray = ih.read(tilt).getData()
            tiltSeriesPixelSize = tilt.getSamplingRate()/10
            deconvTiltSeriesArray = deconvolutionObj.mtf_deconv_wiener(
                                inputTiltSeriesArray, tiltSeriesPixelSize, 
                                psfArrayToDeconv, psfPixelSizeX, kw, pad=20, fc=-1)
        
            outputTiltSeries = ih.createImage()        
            fnOutTiltSeries = self._defineOutputName(i+1)
            k = 0
            for j in range(np.shape(deconvTiltSeriesArray)[0]):
                outputTiltSeries.setData(deconvTiltSeriesArray[j, :, :])
                k += 1
                outputTiltSeries.write((k,fnOutTiltSeries))
            
            self. _sxtFastalign(fnOutTiltSeries, self.firstAngle.get(), self.incrementStep.get())
        
        
    #def createOutputStepFocalSeries(self): 
        
        
      
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        # to ignore validation function of import protocol!!!
        pass
    #def _summary(self):
    #    summary = []
            
    #    return summary
    
    #def _methods(self):
    #    methods = []
    
    #    return methods 
    #--------------------------- UTILS functions --------------------------------------------
    
    def _defineOutputName(self, suffix):
        return self._getExtraPath('deconvolved-tiltSeries_%01d.mrc' % suffix)
    
    def _preDeconvolution(self, inputPSF):
        ih = em.ImageHandler()
        inputPsfArray = ih.read(inputPSF).getData()        
        psfPixelSizeZ = self.inputPsf.get().getZpixelSize() 
        
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
    
    def _createSetOfTiltSeries(self, prefix):
        """ Create a set and set the filename. 
        If the file exists, it will be delete. """
        setFn = self._getPath('%sSeries.sqlite'%prefix)
        # Close the connection to the database if
        # it is open before deleting the file
        cleanPath(setFn)        
        SqliteDb.closeConnection(setFn)        
        setObj = SetOfTiltSeries(filename=setFn)
        return setObj
    
    def _createTiltSeriesMd(self, suffix, tiltImagesStack, tiltAngles):
        outputMd = self._getExtraPath('deconvolvedTiltSeries_%d.xmd'%suffix)
        anglesArray = np.array(tiltAngles.split(','), dtype=np.float)
        mdOut = xmipp.MetaData()            
        for k in range(np.shape(anglesArray)[0]):
            objId = mdOut.addObject()
            mdOut.setValue(xmipp.MDL_IMAGE, "%d@%s" % (k+1,tiltImagesStack), objId)
            mdOut.setValue(xmipp.MDL_ANGLE_TILT, anglesArray[k], objId)
        mdOut.write(outputMd)
        return outputMd
    
    def _sxtFastalign(self, fnStack, firsttiltAngle, incrementStep):
        fnBase = removeExt(fnStack)
        
        #Prealignment file generation               
        tiltxcorrArgs = "-input %s -output %s" % (fnStack, fnBase + '.prexf')
        tiltxcorrArgs += " -first %f" % firsttiltAngle
        tiltxcorrArgs += " -increment %f" % incrementStep
        tiltxcorrArgs += " -rotation 90.0 -sigma1 0.03"
        tiltxcorrArgs += " -radius2 0.25 -sigma2 0.05"        
        self.runJob('tiltxcorr', tiltxcorrArgs) ### iMod needs to be installed
        
        #alignment file generation        
        xftoxgArgs = "-input %s -goutput %s" % (fnBase + '.prexf', fnBase + '.prexg')
        xftoxgArgs += " -nfit 0" 
        self.runJob('xftoxg', xftoxgArgs) 
        
        newstackArgs = "-input %s -output %s" % (fnStack, fnBase + '.preali')
        newstackArgs += " -mode 0 -float 2"
        newstackArgs += " -xform %s" % (fnBase + '.prexg')
        self.runJob('newstack', newstackArgs) 
        
        tiltxcorrArgs = "-input %s -output %s" % (fnBase + '.preali', fnBase + '.fid')
        tiltxcorrArgs += " -first %f" % firsttiltAngle
        tiltxcorrArgs += " -increment %f" % incrementStep
        tiltxcorrArgs += " -prexf %s" % (fnBase + '.prexg')
        tiltxcorrArgs += " -rotation 90.0 -sigma1 0.03"
        tiltxcorrArgs += " -radius2 0.25 -sigma2 0.05" 
        tiltxcorrArgs += " -border 49,49 -size 300,300"
        tiltxcorrArgs += " -LengthAndOverlap 15,4 -overlap 0.33,0.33"       
        self.runJob('tiltxcorr', tiltxcorrArgs)
        
        #Alignment (last step - creating aligned output .mrc)
        tiltalignArgs = "-ModelFile %s" % (fnBase + '.fid')
        tiltalignArgs += " -ImageFile %s" % (fnBase + '.preali')
        tiltalignArgs += " -OutputTransformFile %s" % (fnBase + '.tltxf')
        tiltalignArgs += " -OutputLocalFile %s" % (fnBase + '_local.xf')
        tiltalignArgs += " -OutputTiltFile %s" % (fnBase + '.tlt')        
        tiltalignArgs += " -first %f" % firsttiltAngle
        tiltalignArgs += " -increment %f" % incrementStep
        tiltalignArgs += " -RotationAngle 90.0 -AngleOffset 0.0"
        tiltalignArgs += " -RotOption -1 -RotDefaultGrouping 5"
        tiltalignArgs += " -TiltOption 0 -MagReferenceView 1"
        tiltalignArgs += " -MagOption 0 -MagDefaultGrouping 4"
        tiltalignArgs += " -XStretchOption 0 -XStretchDefaultGrouping 7"
        tiltalignArgs += " -SkewOption 0 -SkewDefaultGrouping 11"
        tiltalignArgs += " -ResidualReportCriterion 3.0 -SurfacesToAnalyze 1"
        tiltalignArgs += " -MetroFactor 0.25 -MaximumCycles 1000"
        tiltalignArgs += " -AxisZShift 0.0 -LocalAlignments 0"
        tiltalignArgs += " -MinFidsTotalAndEachSurface 8,3 -LocalOutputOptions 1,0,1"
        tiltalignArgs += " -LocalRotOption 3 -LocalRotDefaultGrouping 6"
        tiltalignArgs += " -LocalTiltOption 5 -LocalTiltDefaultGrouping 6"
        tiltalignArgs += " -LocalMagReferenceView 1 -LocalMagOption 3"
        tiltalignArgs += " -LocalMagDefaultGrouping 7 -LocalXStretchOption 0"
        tiltalignArgs += " -LocalXStretchDefaultGrouping 7 -LocalSkewOption 0"
        tiltalignArgs += " -LocalSkewDefaultGrouping 11 -BeamTiltOption 0"
        self.runJob('tiltalign', tiltalignArgs)
        
        xfproductArgs = "-in1 %s -in2 %s" % (fnBase + '.prexg', fnBase + '.tltxf')
        xfproductArgs += " -output %s" % (fnBase + '_fid.xf')
        self.runJob('xfproduct', xfproductArgs)
        
        newstackArgs = "-input %s" % fnStack
        newstackArgs += " -output %s" % (fnBase + '_fastAlignedOutput.mrc')
        newstackArgs += " -offset 0,0 -origin -taper 1,0"
        newstackArgs += " -xform" % (fnBase + '_fid.xf')
        self.runJob('newstack', newstackArgs)
        
        
        