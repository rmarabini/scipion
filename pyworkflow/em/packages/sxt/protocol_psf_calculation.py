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
from os.path import basename
from pyworkflow.utils.path import removeExt
from h5py import File
import xmipp
from pyworkflow.utils import getFloatListFromValues
import numpy as np
import pyworkflow.em as em
import pickle
import scipy as sp


class ProtPsfCalculation(Protocol):

    """    
    This protocol is aimed to calculate PSF from input Siemens star pattern.        
    """
    
    _label = 'calculating PSF from SS'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSiemensStar', params.PathParam, 
                      label="Siemens star pattern",  
                      help="Siemens star pattern is the input image or the "
                           "stack of input images with different ZP "
                           "(and same angle) for calculating the MTF.\n"
                           "Note: There is a different MTF for each different "
                           "ZP position.")
        form.addParam('refNumber', params.IntParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Reference number",
                      help="image index of siemens star pattern used to "
                           "calculate de reference image.\n"
                           "By default, nRef = -1  automatically selects central "
                           "image (floor(N/2)+1).")
        form.addParam('fractionOrders', params.NumericListParam,
                      default="1 3",
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Orders",
                      help="Diffraction orders used to calculate MTF profiles.") 
        form.addParam('ringPosition', params.NumericListParam,
                      default="1.5e3 3e3 6e3 12e3",
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Ring positions (nm)",
                      help="Center radius in nm of the void rings of "
                           "Siemens star pattern to be removed from the "
                           "calculation. Default values are set to match "
                           "Xradia manufactured SS pattern")
        form.addParam('pixelSizeX', params.FloatParam, default=5.0,
                      label="PSF pixel size along X axis (nm)")
        form.addParam('pixelSizeZ', params.FloatParam, default=150,
                      label="PSF pixel size along Z axis (nm)") 
                      
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):         
        fnOutPsf = self._defineOutputName()                   
        inputSS = self.inputSiemensStar.get()        
        nRef = self.refNumber.get()
        orders = getFloatListFromValues(self.fractionOrders.get())
        ringPos = getFloatListFromValues(self.ringPosition.get())
               
        self._insertFunctionStep('validateSiemensStar', inputSS)   
        self._insertFunctionStep('getPsfFromSiemensStar', inputSS, nRef, orders, ringPos, fnOutPsf)
        self._insertFunctionStep('createOutputStep', fnOutPsf)            
    #--------------------------- STEPS functions --------------------------------------------
    
    def validateSiemensStar(self, inputSS):
        fnIn = self._getExtraPath(basename(inputSS))                   
        fnStack = removeExt(fnIn) + '.mrc'
        fnOutMd = removeExt(fnIn) + '.xmd'  
        fhHdf5 = File(inputSS, 'r')
        
        imgArray = fhHdf5 ["NXtomo/data/data"][()]                                 
        ih = em.ImageHandler()
        outputImg = ih.createImage()
        i = 0
        for j in range(np.shape(imgArray)[0]):
            outputImg.setData(imgArray[j, :, :])
            i += 1
            outputImg.write((i,fnStack))    
             
        anglesArray = fhHdf5 ["NXtomo/data/rotation_angle"][()]            
        mdOut = xmipp.MetaData()            
        for k in range(np.shape(anglesArray)[0]):
            objId = mdOut.addObject()
            mdOut.setValue(xmipp.MDL_IMAGE, "%d@%s" % (k+1,fnStack), objId)
            mdOut.setValue(xmipp.MDL_ANGLE_TILT, anglesArray[k], objId)
            if k != 0 and anglesArray[k] != anglesArray[k-1]:
                raise Exception ("Selected input file is not a Siemens Star pattern!!!")                 
        mdOut.write(fnOutMd)  
    
    def getPsfFromSiemensStar(self, inputSS, nRef, orders, ringPos, fnOutPsf):        
        fhHdf5 = File(inputSS, 'r')
        imgSS = fhHdf5 ["NXtomo/data/data"][()]
        print "Input Siemense Star pattern dimensions are:\n", np.shape(imgSS)        
        
        imgNumberTotal = np.shape(imgSS)[0]
        imgNumber = np.floor(imgNumberTotal / 2)
        from xpytools.getResolutionfromSiemensStar import MTFgetResolution
        resolutionObj = MTFgetResolution()
        imgSsResolution = resolutionObj.getResolutionfromSiemensStar(imgSS[imgNumber], ringPos)
        print "Resolution of Siemens star image is:\n" , imgSsResolution [0]        
        dx = imgSsResolution [0]
        
        imgSingle = imgSS[imgNumber]
        from xpytools.getMTFfromSiemensStar import MTFfromSiemensStar
        MTFObj = MTFfromSiemensStar()
        mtfOut = MTFObj.getMTFfromSiemensStar(imgSS[77:81], dx, nRef, orders, ringPos)  
        pickle.dump(mtfOut, open(self._defineMtfDicName(), "wb"))
        fx = mtfOut['fx']
        mtf = mtfOut['mtf']
        print "MTF dimension is:\n" , np.shape(mtf)
        
        from xpytools.mtf2psf import MTF2PSFClass
        mtf2psfObj = MTF2PSFClass()
        psfdict = mtf2psfObj.mtf2psf(mtf, fx, self.pixelSizeX.get(), fov=400, fc=-1) 
        pickle.dump(psfdict, open(self._definePsfDicName(), "wb"))        
        psfArray = psfdict['psf']
        psfPixelSize = psfdict['dx']
        print "PSF dimension is:\n" , np.shape(psfArray)
        print "PSF pixelSize(nm) is:\n" , psfPixelSize
        
        ih = em.ImageHandler()
        outputImg = ih.createImage()        
        i = 0
        for j in range(np.shape(psfArray)[0]):
            outputImg.setData(psfArray[j, :, :])
            i += 1
            outputImg.write((i,fnOutPsf))    
        
    def createOutputStep(self, fnOutPsf): #DoF inja hesab shavad .... 
        psfdict = pickle.load(open(self._definePsfDicName(), "rb"))
        psfPixelSizeX = psfdict['dx']
        psfPixelSizeZ = self.pixelSizeZ.get()
        
        
        
        
        #calculating DoF        
        psfArray = psfdict['psf']
        #claculating best psf image to get their mean and use in deconvolution process
        xCenter = np.floor(np.shape(psfArray)[1]/2)
        centerValues = psfArray[:, xCenter, xCenter]
        print "centerValues", centerValues
        thr = 0.3
        thrv = max(centerValues) * thr
        #finding indices based on centerValues curne and thrv      
        #mp: indices array above threshold, 
        mp = np.where(centerValues > thrv)
        print "mp=", mp, "mp[0]=", mp[0]
        zc = np.ceil(np.shape(psfArray)[0]/2)
        print "zc=", zc
        maxValue = max(centerValues)
        print "maxValue=", maxValue
        maxIndex = [i for i, j in enumerate(centerValues) if j == maxValue]
        print "maxIndex", maxIndex
        argth = sp.optimize.fmin(lambda x:abs(np.power(np.sinc(x), 2)-thr),0)
        print "argth=", argth
        #zm: central index
        zm = mp[0][np.floor((mp[0][0] - mp[0][-1])/2)]
        k0 = argth / zm
        print "zm",zm , "k0",k0 
        print "np.shape(psfArray)[0]", np.shape(psfArray)[0] 

        
        #x0 = [0, maxValue, k0, (psfArray[0][zc]-psfArray[0][maxIndex])]
        x0 = [0, maxValue, k0[0], (zc-maxIndex)]
        print "x0", x0
        
        
        #Iapsf = lambda (x,Dz) : np.power(np.multiply(x(1)+x(2),np.sinc(np.matrix(x(3))*np.matrix(Dz-x(4)))),2)
        print "centerValues[mp]", centerValues[mp]
        
        #fminFunc = lambda x : (np.mean(abs(Iapsf(x,mp)-centerValues[mp]), axis=0))
        
        fminFunc = lambda x : (np.mean(abs(self._testFunction(x, mp[0])-centerValues[mp]), axis=0))
        xf = sp.optimize.fmin(fminFunc,x0,maxiter=3000, maxfun=3000)
        print "np.range(np.shape(psfArray)[0])", np.range(np.shape(psfArray)[0])
        
        tapsf = Iapsf(xf,np.range(np.shape(psfArray)[0]))
        
        maxValueTapsf = max(tapsf)
        print "maxValueTapsf=", maxValueTapsf
       
        maxIndexTapsf = [i for i, j in enumerate(tapsf) if j == maxValueTapsf]
        print "maxIndexTapsf", maxIndexTapsf
        
        #dmin = sp.interpolate.interp1d(tapsf[0:maxIndexTapsf], np.shape(psfArray)[0][0:maxIndexTapsf], maxValueTapsf*0.8, kind='cubic')
        #dmax = sp.interpolate.interp1d(tapsf[maxIndexTapsf:-1], np.shape(psfArray)[0][maxIndexTapsf:-1], maxValueTapsf*0.8, kind='cubic')
        print "np.range(np.shape(psfArray)[0])", np.range(np.shape(psfArray)[0])
        dmin = sp.interpolate.interp1d(tapsf[0:maxIndexTapsf], np.range(np.shape(psfArray)[0])[0:maxIndexTapsf], maxValueTapsf*0.8, kind='cubic')
        dmax = sp.interpolate.interp1d(tapsf[maxIndexTapsf:-1], np.range(np.shape(psfArray)[0])[maxIndexTapsf:-1], maxValueTapsf*0.8, kind='cubic')
        print "dmin", dmin, "dmax", dmax
        
        dof = dmax-dmin;
        print "dof", dof

        
        
        
        outPsf = em.Volume()
        outPsf.setLocation(fnOutPsf)
        outPsf.setSamplingRate(psfPixelSizeX * 10)
        #outPsf.setZpixelSize(psfPixelSizeZ)
        #outPsf.setDoF(psfDof)       
        self._defineOutputs(outputPSF=outPsf)
              
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

    def _validate(self):
        errors = []
        hdf5 = self.inputSiemensStar.get().endswith('hdf5')
        if self.inputSiemensStar.get():
            if not hdf5:
                errors.append ("Expected hdf5 files for importing!!!") 
        else:
            errors.append("The path can not be empty!!!")      
        return errors              
    #--------------------------- UTILS functions --------------------------------------------
    
    def _defineOutputName(self):
        return self._getExtraPath('psf.mrc')
    def _defineMtfDicName(self):
        return self._getExtraPath('mtfDic.p')
    def _definePsfDicName(self):
        return self._getExtraPath('psfDic.p')
    def _testFunction(self,x,Dz):
        Iapsf = np.power(np.multiply(x(1)+x(2),np.sinc(np.matrix(x(3))*np.matrix(Dz-x(4)))),2)
        return Iapsf