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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
from pip._vendor.distlib.metadata import Metadata
from xmipp import MetaData
"""
Protocol to perform high-resolution reconstructions
"""

from glob import glob
import math
import numpy as np

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam, EnumParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume
from pyworkflow.em.metadata.utils import getFirstRow, getSize
from os.path import join, exists, split
from pyworkflow.em.packages.xmipp3.convert import createItemMatrix, setXmippAttributes
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
import pyworkflow.em as em
from convert import writeSetOfParticles
import xmipp


class XmippProtReconstructHeterogeneous(ProtClassify3D):
    """3D Reconstruction with heterogeneous datasets"""
    _label = 'significant heterogeneity'

    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Full-size Images", important=True, 
                      pointerClass='SetOfParticles', allowsNull=False,
                      help='Select a set of images at full resolution')
        form.addParam('inputVolumes', PointerParam, label="Initial volumes", important=True,
                      pointerClass='SetOfVolumes',
                      help='Select a set of volumes')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('particleRadius', IntParam, default=-1, 
                     label='Radius of particle (px)',
                     help='This is the radius (in pixels) of the spherical mask covering the particle in the input images')       
        form.addParam('targetResolution', FloatParam, default=8, 
                     label='Target resolution',
                     help='Target resolution to solve for the heterogeneity')    
        form.addParam('computeDiff', BooleanParam, default=False, label="Compute the difference volumes")   
        
        form.addSection(label='Angular assignment')
        form.addParam('numberOfIterations', IntParam, default=3, label='Number of iterations')
        form.addParam('nextMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views', expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0, expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90, expertLevel=LEVEL_ADVANCED)
        form.addParam('numberOfReplicates', IntParam, label="Max. Number of Replicates", default=3, 
                  expertLevel=LEVEL_ADVANCED, help="Significant alignment is allowed to replicate each image up to this number of times")

        form.addParallelSection(threads=1, mpi=8)
    
    def getNumberOfReconstructedVolumes(self):
        return len(self.inputVolumes.get())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        self._insertFunctionStep('convertInputStep', self.inputParticles.getObjId())
        
        self.TsOrig=self.inputParticles.get().getSamplingRate()
        firstIteration = 1
        for self.iteration in range(self.numberOfIterations.get()):
            self.insertIteration(firstIteration+self.iteration)
        self._insertFunctionStep("createOutput")
    
    def insertIteration(self,iteration):
        self._insertFunctionStep('globalAssignment',iteration)
        self._insertFunctionStep('classifyParticles',iteration)
        self._insertFunctionStep('reconstruct',iteration)
        self._insertFunctionStep('postProcessing',iteration)
        if iteration>1:
            self._insertFunctionStep('evaluateConvergence',iteration)

    def readInfoField(self,fnDir,block,label):
        mdInfo = xmipp.MetaData("%s@%s"%(block,join(fnDir,"info.xmd")))
        return mdInfo.getValue(label,mdInfo.firstObject())

    def writeInfoField(self,fnDir,block,label, value):
        mdInfo = xmipp.MetaData()
        objId=mdInfo.addObject()
        mdInfo.setValue(label,value,objId)
        mdInfo.write("%s@%s"%(block,join(fnDir,"info.xmd")),xmipp.MD_APPEND)
    
    def convertInputStep(self, inputParticlesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn)
        self.runJob('xmipp_metadata_utilities','-i %s --fill image1 constant noImage'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "image1=image"'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --fill particleId constant 1'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "particleId=itemId"'%self.imgsFn,numberOfMpi=1)
        imgsFnId=self._getExtraPath('imagesId.xmd')
        self.runJob('xmipp_metadata_utilities','-i %s --operate keep_column particleId -o %s'%(self.imgsFn,imgsFnId),numberOfMpi=1)

        TsCurrent=max(self.TsOrig,self.targetResolution.get()/3)
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        if newXdim<40:
            newXdim=long(40)
            TsCurrent=Xdim*(self.TsOrig/newXdim)
        fnDir = self._getExtraPath()
        self.writeInfoField(fnDir,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
        self.writeInfoField(fnDir,"size",xmipp.MDL_XSIZE,newXdim)

        # Prepare images
        print "Preparing images to sampling rate=",TsCurrent
        fnNewParticles=join(fnDir,"imagesWiener.stk")
        fnNewMetadata = join(fnDir,"imagesWiener.xmd")
        if newXdim!=Xdim:
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim),
                        numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        else:
            self.runJob("xmipp_image_convert","-i %s -o %s --save_metadata_stack %s"%(self.imgsFn,fnNewParticles,fnNewMetadata),numberOfMpi=1)
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())

        row=getFirstRow(fnNewMetadata)
        if row.containsLabel(xmipp.MDL_CTF_MODEL) or row.containsLabel(xmipp.MDL_CTF_DEFOCUSU):
            args="-i %s --sampling_rate %f --correct_envelope"%(fnNewMetadata,TsCurrent)
            if self.inputParticles.get().isPhaseFlipped():
                args+=" --phase_flipped"
            self.runJob("xmipp_ctf_correct_wiener2d",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        
        args="-i %s --sampling %f --fourier low_pass %f"%(fnNewMetadata,TsCurrent,self.targetResolution)
        self.runJob("xmipp_transform_filter",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        
        # Prepare mask
        img = ImageHandler()
        if self.nextMask.hasValue():
            fnMask=join(fnDir,"mask.vol")
            maskObject = self.nextMask.get()
            img.convert(maskObject, fnMask)
            self.runJob('xmipp_image_resize',"-i %s --factor %f"%(fnMask,maskObject.getSamplingRate()/TsCurrent),numberOfMpi=1)
            self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnMask,newXdim),numberOfMpi=1)

        # Prepare volumes
        fnDir=self._getExtraPath('Iter000')
        makePath(fnDir)
        i=1
        for vol in self.inputVolumes.get():
            fnVol=join(fnDir,"volume%02d.mrc"%i)
            img.convert(vol, fnVol)
            TsVol = vol.getSamplingRate()
            if TsVol!=TsCurrent:
                self.runJob('xmipp_image_resize',"-i %s --factor %f"%(fnVol,TsVol/TsCurrent),numberOfMpi=1)
            self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnVol,newXdim),numberOfMpi=1)
            args="-i %s --sampling %f --fourier low_pass %f"%(fnVol,TsCurrent,self.targetResolution)
            self.runJob("xmipp_transform_filter",args,numberOfMpi=1)
            i+=1

    def prepareReferences(self,fnDirPrevious,fnDir,TsCurrent,Xdim):
        fnMask=''
        if self.nextMask.hasValue():
            fnMask=self._getExtraPath("mask.vol")
        for i in range(0,self.getNumberOfReconstructedVolumes()):
            fnPreviousVol=join(fnDirPrevious,"volume%02d.mrc"%(i+1))
            fnReferenceVol=join(fnDir,"volumeRef%02d.mrc"%(i+1))
            copyFile(fnPreviousVol, fnReferenceVol)
            self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%\
                        (fnReferenceVol,self.targetResolution.get(),TsCurrent),numberOfMpi=1)
            R=self.particleRadius.get()
            if R<=0:
                R=self.inputParticles.get().getDimensions()[0]/2*self.TsOrig
            self.runJob('xmipp_transform_mask','-i %s --mask circular -%d'%\
                        (fnReferenceVol,round(R*self.TsOrig/TsCurrent)),numberOfMpi=1)
            self.runJob('xmipp_transform_threshold','-i %s --select below 0 --substitute value 0'%fnReferenceVol,numberOfMpi=1)
            if fnMask!='':
                self.runJob('xmipp_image_operate','-i %s --mult %s'%(fnReferenceVol,fnMask),numberOfMpi=1)

    def calculateAngStep(self,newXdim,TsCurrent,ResolutionAlignment):
        k=newXdim*TsCurrent/ResolutionAlignment # Freq. index
        return math.atan2(1,k)*180.0/math.pi # Corresponding angular step

    def globalAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)

        TsCurrent=self.readInfoField(self._getExtraPath(),"sampling",xmipp.MDL_SAMPLINGRATE)
        newXdim=self.readInfoField(self._getExtraPath(),"size",xmipp.MDL_XSIZE)
        self.prepareReferences(fnDirPrevious,fnDirCurrent,TsCurrent,newXdim)

        # Calculate angular step at this resolution
        angleStep=self.calculateAngStep(newXdim, TsCurrent, self.targetResolution.get())
        angleStep=max(angleStep,5.0)
        self.writeInfoField(fnDirCurrent,"angleStep",xmipp.MDL_ANGLE_DIFF,float(angleStep))
        
        # Generate projections
        fnImgs = self._getExtraPath("imagesWiener.xmd")
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            if not exists(fnAngles):
                fnReferenceVol=join(fnDirCurrent,"volumeRef%02d.mrc"%i)
                fnGallery=join(fnDirCurrent,"gallery%02d.stk"%i)
                fnGalleryXmd=join(fnDirCurrent,"gallery%02d.doc"%i)
                args="-i %s -o %s --sampling_rate %f --perturb %f --sym %s --min_tilt_angle %f --max_tilt_angle %f"%\
                     (fnReferenceVol,fnGallery,angleStep,math.sin(angleStep*math.pi/180.0)/4,self.symmetryGroup,self.angularMinTilt.get(),self.angularMaxTilt.get())
                args+=" --compute_neighbors --angular_distance -1 --experimental_images %s"%fnImgs
                self.runJob("xmipp_angular_project_library",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
                cleanPath(join(fnDirCurrent,"gallery%02d_sampling.xmd"%i))
                cleanPath(join(fnDirCurrent,"gallery%02d_angles.doc"%i))
    
                maxShift=round(self.angularMaxShift.get()*newXdim/100)
                args='-i %s --initgallery %s --maxShift %d --odir %s --dontReconstruct --useForValidation %d'%\
                     (fnImgs,fnGalleryXmd,maxShift,fnDirCurrent,self.numberOfReplicates.get()-1)
                self.runJob('xmipp_reconstruct_significant',args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
                fnAnglesSignificant = join(fnDirCurrent,"angles_iter001_00.xmd")
                if exists(fnAnglesSignificant): 
                    moveFile(fnAnglesSignificant,fnAngles)
                    self.runJob("xmipp_metadata_utilities",'-i %s --operate sort itemId'%fnAngles,numberOfMpi=1)
                    cleanPath(join(fnDirCurrent,"images_iter001_00.xmd"))
                    cleanPath(join(fnDirCurrent,"images_significant_iter001_00.xmd"))
                cleanPath(fnGallery)
                cleanPath(fnGalleryXmd)

    def classifyParticles(self,iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)

        # Gather all angles
        fnAnglesAll = join(fnDirCurrent,"anglesAll.xmd")
        mdVolumes = MetaData()
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnReferenceVol=join(fnDirCurrent,"volumeRef%02d.mrc"%i)
            fnOut = "angles_%02d@%s"%(i,fnAnglesAll)
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            if exists(fnAngles):
                md=MetaData(fnAngles)
            else:
                md=MetaData()
            objId = mdVolumes.addObject()
            mdVolumes.setValue(xmipp.MDL_IMAGE,fnReferenceVol,objId)
            md.write(fnOut,xmipp.MD_APPEND)
        fnVols = join(fnDirCurrent,"referenceVolumes.xmd")
        mdVolumes.write(fnVols)

        # Classify the images
        fnImgsId = self._getExtraPath("imagesId.xmd")
        fnOut = join(fnDirCurrent,"classes.xmd")
        self.runJob("xmipp_classify_significant","--id %s --angles %s --ref %s -o %s"%(fnImgsId,fnAnglesAll,fnVols,fnOut), numberOfMpi=1)
        cleanPath(fnVols)
        cleanPath(fnAnglesAll)

    def reconstruct(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnOut = join(fnDirCurrent,"classes.xmd")
        fnRootVol = join(fnDirCurrent,"class_")
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            self.runJob("xmipp_metadata_split","-i class%06d_images@%s --oroot %s_%06d_"%(i,fnOut,fnRootVol,i),numberOfMpi=1)
            for half in range(1,3):
                args="-i %s_%06d_%06d.xmd -o %s_%02d_half%d.vol --sym %s --weight --thr %d"%(fnRootVol,i,half,fnRootVol,i,half,
                                                                                             self.symmetryGroup,self.numberOfThreads.get())
                self.runJob("xmipp_reconstruct_fourier",args,numberOfMpi=self.numberOfMpi.get())
                cleanPath("%s_%06d_%06d.xmd"%(fnRootVol,i,half))

    def postProcessing(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(self._getExtraPath(),"sampling",xmipp.MDL_SAMPLINGRATE)
        fnRootVol = join(fnDirCurrent,"class_")

        fnMask=''
        if self.nextMask.hasValue():
            fnMask=self._getExtraPath("mask.vol")
        else:
            R=self.particleRadius.get()
            if R<=0:
                R=self.inputParticles.get().getDimensions()[0]/2*self.TsOrig
            fnMask=self._getExtraPath("mask.mrc")
            self.runJob('xmipp_transform_mask','-i %s_000001_half1.vol --mask circular -%d --create_mask %s'%\
                        (fnRootVol,round(R*self.TsOrig/TsCurrent),fnMask),numberOfMpi=1)
        fnCentered = join(fnDirCurrent,"volumeCentered.mrc")
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            # Align the two volumes
            fnVol1="%s_%02d_half1.vol"%(fnRootVol,i)
            fnVol2="%s_%02d_half2.vol"%(fnRootVol,i)
            fnVolAvg=join(fnDirCurrent,"volume%02d.mrc"%i)
            self.runJob('xmipp_image_operate','-i %s --plus %s -o %s'%(fnVol1,fnVol2,fnVolAvg),numberOfMpi=1)
            self.runJob('xmipp_image_operate','-i %s --mult 0.5'%fnVolAvg,numberOfMpi=1)
            self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol1),numberOfMpi=1)
            self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol2),numberOfMpi=1)

            # Remove untrusted background voxels
            fnRootRestored=join(fnDirCurrent,"volumeRestored")
            args='--i1 %s --i2 %s --oroot %s --denoising 1 --mask binary_file %s'%(fnVol1,fnVol2,fnRootRestored,fnMask)
            self.runJob('xmipp_volume_halves_restoration',args,numberOfMpi=1)
            moveFile("%s_restored1.vol"%fnRootRestored,fnVol1)
            moveFile("%s_restored2.vol"%fnRootRestored,fnVol2)

            # Filter bank denoising
            args='--i1 %s --i2 %s --oroot %s --filterBank 0.01 --mask binary_file %s'%(fnVol1,fnVol2,fnRootRestored,fnMask)
            self.runJob('xmipp_volume_halves_restoration',args,numberOfMpi=1)
            moveFile("%s_restored1.vol"%fnRootRestored,fnVol1)
            moveFile("%s_restored2.vol"%fnRootRestored,fnVol2)
            cleanPath("%s_filterBank.vol"%fnRootRestored)
     
            # Laplacian Denoising
            args = "-i %s --retinex 0.95 "+fnMask
            self.runJob('xmipp_transform_filter',args%fnVol1,numberOfMpi=1)
            self.runJob('xmipp_transform_filter',args%fnVol2,numberOfMpi=1)

            # Blind deconvolution
            args='--i1 %s --i2 %s --oroot %s --deconvolution 1 --mask binary_file %s'%(fnVol1,fnVol2,fnRootRestored,fnMask)
            self.runJob('xmipp_volume_halves_restoration',args,numberOfMpi=1)
            moveFile("%s_restored1.vol"%fnRootRestored,fnVol1)
            moveFile("%s_restored2.vol"%fnRootRestored,fnVol2)
            self.runJob("xmipp_image_convert","-i %s_convolved.vol -o %s -t vol"%(fnRootRestored,fnVolAvg),numberOfMpi=1)
            cleanPath("%s_convolved.vol"%fnRootRestored)
            cleanPath("%s_deconvolved.vol"%fnRootRestored)

            # Difference evaluation and production of a consensus average
            args='--i1 %s --i2 %s --oroot %s --difference 2 2 --mask binary_file %s'%(fnVol1,fnVol2,fnRootRestored,fnMask)
            self.runJob('xmipp_volume_halves_restoration',args,numberOfMpi=1)
            self.runJob("xmipp_image_convert","-i %s_avgDiff.vol -o %s -t vol"%(fnRootRestored,fnVolAvg),numberOfMpi=1)
            cleanPath("%s_avgDiff.vol"%fnRootRestored)
            moveFile("%s_restored1.vol"%fnRootRestored,fnVol1)
            moveFile("%s_restored2.vol"%fnRootRestored,fnVol2)
            
            # FSC
            self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol1,fnMask),numberOfMpi=1)
            self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol2,fnMask),numberOfMpi=1)
            self.runJob('xmipp_transform_threshold','-i %s --select below 0 --substitute value 0 '%fnVol1,numberOfMpi=1)
            self.runJob('xmipp_transform_threshold','-i %s --select below 0 --substitute value 0 '%fnVol2,numberOfMpi=1)

            fnFsc=join(fnDirCurrent,"fsc%02d.xmd"%i)
            self.runJob('xmipp_resolution_fsc','--ref %s -i %s -o %s --sampling_rate %f'\
                        %(fnVol1,fnVol2,fnFsc,TsCurrent),numberOfMpi=1)
            self.runJob('xmipp_transform_filter','-i %s --fourier fsc %s --sampling %f'%(fnVolAvg,fnFsc,TsCurrent),numberOfMpi=1)
            cleanPath(fnVol1)
            cleanPath(fnVol2)

            self.runJob('xmipp_image_header','-i %s --sampling_rate %f'%(fnVolAvg,TsCurrent),numberOfMpi=1)
            
            if i==1:
                copyFile(fnVolAvg,fnCentered)
            else:
                self.runJob("xmipp_image_operate","-i %s --plus %s"%(fnCentered,fnVolAvg),numberOfMpi=1)

        # Align all volumes
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnVolAvg=join(fnDirCurrent,"volume%02d.mrc"%i)
            self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnCentered,fnVolAvg),numberOfMpi=1)
        cleanPath(fnCentered)
        
        # Rewrite the first block of the classes.xmd with the representative
        fnClasses = "classes@"+join(fnDirCurrent,"classes.xmd")
        classesmd = md.MetaData(fnClasses)
        for objId in classesmd:
            ref3d = classesmd.getValue(md.MDL_REF3D,objId)
            classesmd.setValue(md.MDL_IMAGE,join(fnDirCurrent,"volume%02d.mrc"%ref3d),objId)
        classesmd.write(fnClasses,md.MD_APPEND)
        
        # Write the images.xmd
        mdAll = md.MetaData()
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            mdi = md.MetaData("class%06d_images@"%i+join(fnDirCurrent,"classes.xmd"))
            mdAll.unionAll(mdi)
        mdAll.write(join(fnDirCurrent,"images.xmd"))

    def evaluateConvergence(self,iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnIntersection = join(fnDirCurrent,"intersection.xmd")
        fnUnion = join(fnDirCurrent,"union.xmd")
        fnAngleDistance = join(fnDirCurrent,"angle_distance")
        N=self.getNumberOfReconstructedVolumes()
        coocurrence = np.zeros((N, N))
        sizeClasses = np.zeros((N, N))
        for i in range(1,N+1):
            for j in range(1,N+1):
                fnCurrent = "class%06d_images@%s"%(i,join(fnDirCurrent,"classes.xmd"))
                fnPrevious = "class%06d_images@%s"%(j,join(fnDirPrevious,"classes.xmd"))
                self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s itemId -o %s"%(fnCurrent,fnPrevious,fnIntersection),numberOfMpi=1)
                self.runJob("xmipp_metadata_utilities","-i %s --set union %s itemId -o %s"%(fnCurrent,fnPrevious,fnUnion),numberOfMpi=1)
                
                sizeIntersection = float(md.getSize(fnIntersection))
                sizeUnion = float(md.getSize(fnUnion))
                sizeClasses[i-1,j-1]= sizeIntersection                
                coocurrence[i-1,j-1]=sizeIntersection/sizeUnion
                
                if i==j:
                    self.runJob("xmipp_angular_distance","--ang1 %s --ang2 %s --oroot %s --sym %s --check_mirrors --compute_weights 1 particleId 0.5"%\
                                (fnPrevious,fnCurrent,fnAngleDistance,self.symmetryGroup),numberOfMpi=1)
                    distances = md.MetaData(fnAngleDistance+"_weights.xmd")
                    angDistance=distances.getColumnValues(md.MDL_ANGLE_DIFF)
                    avgAngDistance = reduce(lambda x, y: x + y, angDistance) / len(angDistance)
                    shiftDistance=distances.getColumnValues(md.MDL_SHIFT_DIFF)
                    avgShiftDistance = reduce(lambda x, y: x + y, shiftDistance) / len(shiftDistance)
                    print("Class %d: average angular diff=%f      average shift diff=%f"%(i,avgAngDistance,avgShiftDistance))
        
        print("Size of the intersections")
        print(sizeClasses)
        print(' ')
        print('Stability of the classes (coocurrence)')
        print(coocurrence)
                    
        cleanPath(fnIntersection)
        cleanPath(fnUnion)
        cleanPath(fnAngleDistance+"_weights.xmd")            

    def createOutput(self):
        # get last iteration
        fnIterDir=glob(self._getExtraPath("Iter*"))
        lastIter=len(fnIterDir)-1
        self.fnLastDir=self._getExtraPath("Iter%03d"%lastIter)

        fnLastImages=join(self.fnLastDir,"images.xmd")
        if not exists(fnLastImages):
            raise Exception("The file %s does not exist"%fnLastImages)
        partSet = self.inputParticles.get()
        self.Ts=self.readInfoField(self._getExtraPath(),"sampling",xmipp.MDL_SAMPLINGRATE)
        self.scaleFactor=self.Ts/partSet.getSamplingRate()

        classes3D = self._createSetOfClasses3D(partSet)
        # Let use an special iterator to skip missing particles
        # discarded by classification (in case of cl2d)
        setMdIter = md.SetMdIterator(fnLastImages,
                                     sortByLabel=md.MDL_PARTICLE_ID,
                                     updateItemCallback=self._updateParticle)
        
        classes3D.classifyItems(updateItemCallback=setMdIter.updateItem,
                             updateClassCallback=self._updateClass)
        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles, classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(self.Ts)
        
        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)
        
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles, volumes)

    def _updateParticle(self, particle, row):
        particle.setClassId(row.getValue(md.MDL_REF3D))
        row.setValue(xmipp.MDL_SHIFT_X,row.getValue(xmipp.MDL_SHIFT_X)*self.scaleFactor)
        row.setValue(xmipp.MDL_SHIFT_Y,row.getValue(xmipp.MDL_SHIFT_Y)*self.scaleFactor)
        setXmippAttributes(particle, row, xmipp.MDL_SHIFT_X, 
                           xmipp.MDL_SHIFT_Y, xmipp.MDL_ANGLE_ROT, 
                           xmipp.MDL_ANGLE_TILT, xmipp.MDL_ANGLE_PSI, 
                           xmipp.MDL_MAXCC, xmipp.MDL_WEIGHT)
        createItemMatrix(particle, row, align=em.ALIGN_PROJ)

    def _updateClass(self, item):
        classId = item.getObjId()
        item.setAlignment3D()
        item.setSamplingRate(self.Ts)
        item.getRepresentative().setFileName(join(self.fnLastDir,"volume%02d.mrc"%classId))
        