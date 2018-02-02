# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), March 2014
# *           Javier Mota Garcia (jmota@cnb.csic.es), January 2018
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

from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import XmippProtConvertToPseudoAtomsBase
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from protocol_nma_base import *
from pyworkflow.utils.path import createLink, cleanPath
from pyworkflow.protocol.params import BooleanParam, MultiPointerParam
from xmipp import MetaData, MDL_NMA, MDL_ENABLED, MDL_NMA_MINRANGE, MDL_NMA_MAXRANGE, MDL_NMA_ABS
import xmipp
from pyworkflow.em.packages.xmipp3 import XmippMdRow
from convert import rowToMode

class XmippProtNMAChoose(XmippProtConvertToPseudoAtomsBase, XmippProtNMABase):
    """ Protocol for choosing a volume to construct an NMA analysis """
    _label = 'NMA selection'
    def __init__(self, **args):
        XmippProtConvertToPseudoAtomsBase.__init__(self, **args)
        XmippProtNMABase.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputRefVolume', PointerParam,
                      pointerClass='Volume',
                      label="Input reference volume", important=True,
                      help='Select the volume you have calculated the normal modes')
        form.addParam('inputVolumes', MultiPointerParam,
                      pointerClass='Volume',
                      label="Input volume(s)", important=True,
                      help='Select one or more volumes (Volume or SetOfVolumes)\n'
                           'to be aligned against the reference volume.')
        form.addParam('inputModes', PointerParam, label="Input normal modes", important=True,
                      pointerClass = 'SetOfNormalModes')
        form.addParam('pseudoAtomRadius', FloatParam, default=1,
                      label='Pseudoatom radius (vox)',
                      help='Pseudoatoms are defined as Gaussians whose \n'
                           'standard deviation is this value in voxels')
        #form.addParam('alignVolumes', BooleanParam, label="Align volumes", default=False,
        #              help="Align deformed PDBs to volume to maximize match")
        #XmippProtConvertToPseudoAtomsBase._defineParams(self,form)
        form.addParallelSection(threads=4, mpi=1)    

        #form.addSection(label='Normal Mode Analysis')
        #XmippProtNMABase._defineParamsCommon(self,form)
             
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        inputModes = self.inputModes.get().getFileName()
        inputModes = inputModes.split(".")
        self.Modes = inputModes[0]+".xmd"
        self.pseudoatoms = self.inputModes.get().getPdb().getFileName()

        print self.inputModes.get().getPdb().getVolume()
        RefVolume = self.inputRefVolume.get().getFileName()
        print self.inputModes.get().getPdb()

        filenames=[]
        self.files = 0
        for inputStructure in self.inputVolumes:
            filenames.append(inputStructure.get().getFileName())
            self.sampling = inputStructure.get().getSamplingRate()
            self.files += 1

        deps = []
        for volCounter in range(1,len(filenames)+1):
            fnIn=filenames[volCounter-1]
            outVolFn = self._getPath('outputRigidAlignment_vol_%s_to_%d.vol' % ("Ref", volCounter))
            self._insertFunctionStep('alignVolumeStep', RefVolume, fnIn, outVolFn, volCounter)
            args="-i %s --pdb %s --modes %s --sampling_rate %f -o %s --fixed_Gaussian %f --opdb %s"%\
                 (filenames[volCounter-1],self.pseudoatoms, \
                  self.Modes,self.sampling,\
                  self._getExtraPath('alignment_%s_%02d.xmd'%("Ref",volCounter)),\
                  self.sampling*self.pseudoAtomRadius.get(),
                  self._getExtraPath('alignment_%s_%02d.pdb'%("Ref",volCounter)))
            #if self.alignVolumes.get():
             #   args+=" --alignVolumes"
            stepId=self._insertRunJobStep("xmipp_nma_alignment_vol",args)
            deps.append(stepId)
            
        self._insertFunctionStep('evaluateDeformationsStep',prerequisites=deps)
        

    #--------------------------- Step functions --------------------------------------------
    def alignVolumeStep(self, refFn, inVolFn, outVolFn, volId):
        args = "--i1 %s --i2 %s --apply %s" % (refFn, inVolFn, outVolFn)
        args += " --local --rot 0 0 1 --tilt 0 0 1 --psi 0 0 1 -x 0 0 1 -y 0 0 1 -z 0 0 1 --dontScale"
        args += " --copyGeo %s" % (
                self._getExtraPath('transformation-matrix_vol%06d.txt'%volId))
        self.runJob("xmipp_volume_align", args)
    
    def evaluateDeformationsStep(self):
        N = self.files
        import numpy
        distances=numpy.zeros([N])
        pdb=open(self.pseudoatoms).readlines()
        for volCounter in range(1,N+1):
            davg=0.
            Navg=0.
            pdb2=open(self._getExtraPath('alignment_%s_%02d.pdb'%("Ref",volCounter))).readlines()
            for i in range(len(pdb)):
                line1=pdb[i]
                if line1.startswith("ATOM"):
                    line2=pdb2[i]
                    x1=float(line1[30:37])
                    y1=float(line1[38:45])
                    z1=float(line1[46:53])
                    x2=float(line2[30:37])
                    y2=float(line2[38:45])
                    z2=float(line2[46:53])
                    dx=x1-x2
                    dy=y1-y2
                    dz=z1-z2
                    d=math.sqrt(dx*dx+dy*dy+dz*dz)
                    davg+=d
                    Navg+=1
            if Navg>0:
                davg/=Navg
            distances[volCounter-1]=davg
        distances=0.5*(distances+numpy.transpose(distances))
        numpy.savetxt(self._getPath('distances.txt'),distances)
        distances1D=numpy.mean(distances,axis=0)
        print("Average distance to rest of volumes=",distances1D)
        imin=numpy.argmin(distances1D)
        print("The volume in the middle is pseudoatoms_%02d.pdb"%(imin+1))
        #createLink(self._getPath("pseudoatoms_%02d.pdb"%(imin+1)),self._getPath("pseudoatoms.pdb"))
        createLink(self._getPath("modes_%02d.xmd"%(imin+1)),self._getPath("modes.xmd"))
        createLink(self._getExtraPath("pseudoatoms_%02d_distance.hist"%(imin+1)),self._getExtraPath("pseudoatoms_distance.hist"))

        # Measure range
        minDisplacement= 1e38*numpy.ones([self.inputModes.get().getSize(),1])
        maxDisplacement=-1e38*numpy.ones([self.inputModes.get().getSize(),1])
        print "Error1"
        absDisplacement = 1e38*numpy.ones([self.inputModes.get().getSize(),1])
        print "Error2"
        mdNMA=MetaData(self.Modes)
        for volCounter in range(1,N+1):
            if volCounter!=imin+1:
                md=MetaData(self._getExtraPath("alignment_%s_%02d.xmd"%("Ref",volCounter)))
                displacements=md.getValue(MDL_NMA, md.firstObject())
                print displacements
                idx1=0
                idx2=0
                for idRow in mdNMA:
                    if mdNMA.getValue(MDL_ENABLED,idRow)==1:
                        minDisplacement[idx2]=min(minDisplacement[idx2],displacements[idx1])
                        maxDisplacement[idx2]=max(maxDisplacement[idx2],displacements[idx1])
                        absDisplacement[idx2]=abs(maxDisplacement[idx2])
                        idx1+=1
                    else:
                        minDisplacement[idx2]=0
                        maxDisplacement[idx2]=0
                        absDisplacement[idx2]=0
                    idx2+=1
        print "Error3"
        idx2=0
        for idRow in mdNMA:
            mdNMA.setValue(MDL_NMA_MINRANGE,float(minDisplacement[idx2]),idRow)
            mdNMA.setValue(MDL_NMA_MAXRANGE,float(maxDisplacement[idx2]),idRow)
            mdNMA.setValue(MDL_NMA_ABS,float(absDisplacement[idx2]),idRow)
            idx2+=1
        mdNMA.write(self._getPath("modes.xmd"))
        print "Error4"
        # Create output
        volCounter=0

        '''for inputStructure in self.inputVolumes:
            if volCounter==imin:
                print("The corresponding volume is %s"%(inputStructure.get().getFileName()))
                finalStructure=inputStructure
                break
            volCounter+=1'''

        '''pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputVolumes.get(), pdb)
        modes = SetOfNormalModes(filename=self._getPath('modes.sqlite'))
        self._defineOutputs(outputModes=modes)
        self._defineSourceRelation(self.inputVolumes,modes)'''
        self._insertFunctionStep('createOutputStep')

        # ToDo: the self.outputPdb should be a Pointer, not an object
#         self._defineSourceRelation(self.outputPdb, self.outputModes)

    #--------------------------- OUTPUT step -----------------------------------------------

    def createOutputStep(self):
        modes = SetOfNormalModes(filename=self._getPath('modes.sqlite'))
        md = xmipp.MetaData(self._getPath('modes.xmd'))
        row = XmippMdRow()

        for objId in md:
            row.readFromMd(md, objId)
            mode = rowToMode(row)
            mode.nmaMin = Float(row.getValue(MDL_NMA_MINRANGE))
            mode.nmaMax = Float(row.getValue(MDL_NMA_MAXRANGE))
            mode.nmaAbs = Float(row.getValue(MDL_NMA_ABS))
            modes.append(mode)
        inputPdb = self.inputModes.get().getPdb()
        modes.setPdb(inputPdb)
        self._defineOutputs(outputModes=modes)
        self._defineSourceRelation(self.inputRefVolume, modes)

        '''inputVol = self.inputRefVolume.get()
        volume = Volume()
        volume.setFileName(self._getExtraPath("pseudoatoms_nma_selection"))
        volume.setSamplingRate(inputVol.getSamplingRate())

        pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
        pdb.setVolume(volume)
        #self.createChimeraScript(inputVol, pdb)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputRefVolume, pdb)'''

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append('Pseudoatom radius (voxels): %f'%self.pseudoAtomRadius.get())
        summary.append('Approximation target error (%%): %f'%self.pseudoAtomTarget.get())
        return summary

    def _methods(self):
        summary = []
#        summary.append('We converted the volume %s into a pseudoatomic representation with Gaussian atoms (sigma=%f A and a target error'\
#                       ' of %f%%) [Nogales2013].'%(self.inputStructure.get().getNameId(),
#                                     self.pseudoAtomRadius.get()*self.inputStructure.get().getSamplingRate(),
#                                     self.pseudoAtomTarget.get()));
#        if self.hasAttribute('outputPdb'):
#            summary.append('We refer to the pseudoatomic model as %s.'%self.outputPdb.getNameId())
        return summary

    def _citations(self):
        return ['Nogales2013','Jin2014']
        