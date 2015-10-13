# **************************************************************************
# *
# * Author:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
from pyworkflow.em.data_tomo import CTF3DModel
"""
This script is a re-implementation of 'prepare_subtomograms.py' script that
was written by Tanmay Bharat to support sub-tomogram averaging in RELION.
"""
import sys
from os.path import basename, join
from pyworkflow.protocol import params
import pyworkflow.utils.path as pwutils
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol, ProtProcessTomograms, STEPS_PARALLEL, LEVEL_ADVANCED
# from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER
from imodpath import CTFFIND_PATH, CTFFIND4_PATH


class ProtPrepareSubtomograms(ProtProcessTomograms):
    """sub-tomogram averaging in RELION
    """
    _label = 'subtomogram averaging'
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addParam('inputTomoCoords', params.PointerParam, important=True,
                      label=Message.LABEL_INPUT_CRD_TOM,
                       pointerClass='SetOfTomoCoordinates')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
              label='CTF Downsampling factor',
              help='Set to 1 for no downsampling. Non-integer downsample factors are possible. '
              'This downsampling is only used for estimating the CTF and it does not affect '
              'any further calculation. Ideally the estimation of the CTF is optimal when '
              'the Thon rings are not too concentrated at the origin (too small to be seen) '
              'and not occupying the whole power spectrum (since this downsampling might '
              'entail aliasing).')
        form.addParam('useCftfind4', params.BooleanParam, default=True,
              label="Use ctffind4 to estimate the CTF?",
              help='If is true, the protocol will use ctffind4 instead of ctffind3')
        form.addParam('astigmatism', params.FloatParam, default=1000.0,
              label='Expected (tolerated) astigmatism (A)')
        form.addParam('findPhaseShift', params.BooleanParam, default=False,
              label="Find additional phase shift?", condition='useCftfind4',
              expertLevel=params.LEVEL_ADVANCED,)
        
        line = form.addLine('Resolution',
                            help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                                 'These cut-offs prevent the typical peak at the center of the PSD and high-resolution'
                                 'terms where only noise exists, to interfere with CTF estimation. The default lowest '
                                 'value is 0.05 but for micrographs with a very fine sampling this may be lowered towards 0.'
                                 'The default highest value is 0.35, but it should '+'be increased for micrographs with '
                                 'signals extending beyond this value. However, if your micrographs extend further than '
                                 '0.35, you should consider sampling them at a finer rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05,
                      label='Lowest' )
        line.addParam('highRes', params.FloatParam, default=0.35,
                      label='Highest')
        # Switched (microns) by 'in microns' by fail in the identifier with jquery
        line = form.addLine('Defocus search range (microns)', expertLevel=LEVEL_ADVANCED,
                            help='Select _minimum_ and _maximum_ values for defocus search range (in microns).'
                                 'Underfocus is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=0.25, 
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=10.,
                      label='Max')
        form.addParam('stepFocus', params.FloatParam, default=2000.0,
              label='Defocus step size for search (A)')
        form.addParam('windowSize', params.IntParam, default=256, expertLevel=LEVEL_ADVANCED,
                      label='Window size',
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier.')
#         form.addSection(label="Reconstruct")
#         form.addParam('filesPath', params.PathParam, 
#               label="Directory of the coordinates",
#               help="Directory with the files you want to import.\n\n"
#                    "The path can also contain wildcards to select\n"
#                    "from several folders.\n\n"
#                    "For example:\n"
#                    "  ~/Coordinates/\n"
#                    "  data/day??_coords/")
#         form.addParam('filesPattern', params.StringParam, default="*.coord",
#                       label='Coordinates pattern', 
#                       help="Pattern of the coordinates to be imported.\n\n"
#                            "The pattern can contain standard wildcards such as\n"
#                            "*, ?, etc.\n"
#                            "The coordinates files must have the same name"
#                            "as the tomograms")
#         form.addParam('boxSize', params.IntParam, default=0,
#                       label='Particle box size', validators=[params.Positive],
#                       help='In pixels. The box size is the size of the boxed subtomograms.')
        
        form.addParallelSection(threads=3, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to preprocess the tomograms
        """
        ctfDepsList = []
        reconsDeps = []
        self.inputCoords = self.inputTomoCoords.get()
        self.tomoSet = self.inputCoords.getTomoRecs().getTomograms()
        
        for tomo in self.tomoSet:
            tomoFn = tomo.getFileName()
            numbOfMics = tomo.getDim()[3]
            
            extractDeps = self._insertFunctionStep("extractTiltAnglesStep", tomoFn, prerequisites=ctfDepsList)
            for i in range(1, numbOfMics+1):
                tiltImgDeps = self._insertFunctionStep("extractTiltImgStep", tomo.getFileName(), i, prerequisites=[extractDeps])
                ctfDeps = self._insertFunctionStep("estimateCtfStep", tomoFn, i, prerequisites=[tiltImgDeps])
                ctfDepsList.append(ctfDeps)
            writeDeps = self._insertFunctionStep("writeStarCtf3DStep", tomo.getFileName(), numbOfMics, prerequisites=ctfDepsList)
            
            for coord in self.inputCoords:
                reconsDep = self._insertFunctionStep("reconstructCtf3DStep", tomoFn, coord.getObjId(), prerequisites=[writeDeps])
                reconsDeps.append(reconsDep)
            self._insertFunctionStep("createOutputStep", prerequisites=reconsDeps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def extractTiltAnglesStep(self,  tomoFn):
        from imodpath import EXTRACTTILTS_PATH
        
        extractProg = EXTRACTTILTS_PATH
        
        args = '-InputFile %(tomogram)s -tilts -OutputFile %(tiltangles)s'
        params = {"tomogram" : tomoFn, 
                  "tiltangles" : self._getFnPath(tomoFn, "tiltAngles.txt")
                  }
        
        self.runJob(extractProg, args % params)
    
    def extractTiltImgStep(self, tomoFn, i):
        from imodpath import NEWSTACK_PATH
        stackProg = NEWSTACK_PATH
        
        args = '-secs %(numStk)d %(tomogram)s %(tomoStck)s'
        params = {"tomogram" : tomoFn,
                  "numStk" : i-1,
                  "tomoStck" : self._getImgName(tomoFn, i)
                  }
        
        self.runJob(stackProg, args % params)
    
    def estimateCtfStep(self, tomoFn, i):
        """ Run ctffind, 3 or 4, with required parameters """
        args, program, params = self._prepareCommand()
        downFactor = self.ctfDownFactor.get()
        micFn = self._getImgName(tomoFn, i)
        if downFactor != 1:
            from pyworkflow.em.packages import xmipp3
            
            micFnMrc = self._getTmpPath(basename(micFn))
            self.runJob("xmipp_transform_downsample","-i %s -o %s --step %f --method fourier" % (micFn, micFnMrc, downFactor), env=xmipp3.getEnviron())
            
        else:
            micFnMrc = micFn
        
        # Update _params dictionary
        params['micFn'] = micFnMrc
        params['micDir'] = self._getTomoPath(tomoFn)
        params['ctffindPSD'] = self._getFnPath(tomoFn, '%s_ctfEstimation.mrc' % pwutils.removeBaseExt(micFn))
        params['ctffindOut'] = self._getFnPath(tomoFn, '%s_ctfEstimation.txt' % pwutils.removeBaseExt(micFn))
        
        try:
            self.runJob(program, args % params)
        except Exception, ex:
            print >> sys.stderr, "ctffind has failed with micrograph %s" % micFnMrc
        pwutils.cleanPath(micFnMrc)
    
    def writeStarCtf3DStep(self, tomoFn, numbOfMics):
        import math
        import pyworkflow.em.metadata as md
        
#         tomoSet = self.inputTomograms.get()
        sizeX, _, _, sizeZ = self.tomoSet.getFirstItem().getDim()
        
#         micKey = pwutils.removeBaseExt(tomoFn)
        
#         coordFn = self._getCoordsFns().get(micKey)
#         numOfCoords = self._getNumberOfCoords(tomoFn)
        
        
#         for i in range(1, numOfCoords+1):
        for coord in self.inputCoords:
            
            ctf3DStar = self._getCtfStar(tomoFn, coord.getObjId())
            ctf3DMd = md.MetaData()
            
#             line = self._getCoords(coordFn, i)
            
            for micNumb in range(1, numbOfMics+1):
                tiltRow = md.Row()
                objId = ctf3DMd.addObject()
                
                aveDefocus = self._getAveDefocus(tomoFn, micNumb)
                tiltAngleGrads = self._getTiltAngles(tomoFn, micNumb)
                tiltAngleRads = float(tiltAngleGrads)*math.pi/180
                
                xTomo = float(coord.getX() - (sizeX/2)) * self.tomoSet.getSamplingRate()
                zTomo = float(coord.getZ() - (sizeZ/2)) * self.tomoSet.getSamplingRate()

                # Calculating the height difference of the particle from the tilt axis
                xImg = (xTomo*(math.cos(tiltAngleRads))) + (zTomo*(math.sin(tiltAngleRads)))
                deltaD = xImg*math.sin(tiltAngleRads)
                partDef = aveDefocus + deltaD
                
                # Weighting the 3D CTF model using the tilt dependent scale factor and the dose dependent B-Factor
                tiltScale = math.cos(abs(tiltAngleRads))
                tiltList = self._getTiltList(tomoFn)
                tiltStep = (max(tiltList) - min(tiltList))/(numbOfMics-1)
                bestTiltDiff = tiltStep + 0.5
                
                for k in range(len(tiltList)):
                    tiltDiff = abs(tiltAngleGrads-tiltList[k])
                    
                    if tiltDiff < (bestTiltDiff+0.25):
                        if tiltDiff < bestTiltDiff:
                            bestTiltDiff = tiltDiff
                            currAcumDose = self._getDose(k+1, numbOfMics)
                    
                doseWeight = currAcumDose * self.tomoSet.getBfactor()
                
                tiltRow.setValue(md.RLN_CTF_DEFOCUSU, partDef)
                tiltRow.setValue(md.RLN_CTF_VOLTAGE, self.tomoSet.getAcquisition().getVoltage())
                tiltRow.setValue(md.RLN_CTF_CS, self.tomoSet.getAcquisition().getSphericalAberration())
                tiltRow.setValue(md.RLN_CTF_Q0, self.tomoSet.getAcquisition().getAmplitudeContrast())
                tiltRow.setValue(md.RLN_ORIENT_ROT, 0.0)
                tiltRow.setValue(md.RLN_ORIENT_TILT, tiltAngleGrads)
                tiltRow.setValue(md.RLN_ORIENT_PSI, 0.0)
                tiltRow.setValue(md.RLN_CTF_BFACTOR, doseWeight)
                tiltRow.setValue(md.RLN_CTF_SCALEFACTOR, tiltScale)
                tiltRow.writeToMd(ctf3DMd, objId)
            ctf3DMd.write("data_images@%s" % ctf3DStar)
    
    def reconstructCtf3DStep(self, tomoFn, coordNum):
        from pyworkflow.em.packages import relion
        program = "relion_reconstruct"
        sampling = self.tomoSet.getSamplingRate()
        params = {"sampling" : sampling,
                  "ctfStar" : self._getCtfStar(tomoFn, coordNum),
                  "ctf3D" : self._getCtf3D(tomoFn, coordNum),
                  "boxSize" : self.boxSize.get()
                  }
        
        args = " --i %(ctfStar)s --o %(ctf3D)s --reconstruct_ctf %(boxSize)d --angpix %(sampling)f"
        self.runJob(program, args % params, env=relion.getEnviron())
    
    def createOutputStep(self):
        ctf3DSet = self._createSetOfCTF3D(self.inputCoords)
        for tomo in self.tomoSet:
            tomoFn = tomo.getFileName()
            for coord in self.inputCoords:
                ctf3D = CTF3DModel()        
                ctf3D.setCtfFile(self._getCtf3D(tomoFn, coord.getObjId()))
                ctf3D.setTomoCoordinate(coord)
                ctf3DSet.append(ctf3D)
        
        self._defineOutputs(ctf3DSet)
        self._defineCtfRelation(self.inputTomoCoords, ctf3DSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors    
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        methods = ''
        return [methods]
        
    def _citations(self):
        pass
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getTomoPath(self, tomoFn):
        tomoBaseDir = pwutils.removeBaseExt(tomoFn)
        pwutils.makePath(self._getExtraPath(tomoBaseDir))
        return self._getExtraPath(tomoBaseDir)
    
    def _getFnPath(self, tomoFn, fn):
        return join(self._getTomoPath(tomoFn), basename(fn))
    
    def _getImgName(self, tomoFn, i):
        baseImgFn = pwutils.removeBaseExt(tomoFn) + "_%03d.mrc" %i
        return join(self._getTomoPath(tomoFn), baseImgFn)
    
    def _getCtfStar(self, tomoFn, i):
        ctfDir = join(self._getTomoPath(tomoFn), "Ctf3D")
        pwutils.makePath(ctfDir)
        baseCtfFn = pwutils.removeBaseExt(tomoFn) + "_ctf_%06d.star" %i
        return join(ctfDir, baseCtfFn)
    
    def _getCtf3D(self, tomoFn, i):
        ctfDir = join(self._getTomoPath(tomoFn), "Particles")
        pwutils.makePath(ctfDir)
        baseCtfFn = pwutils.removeBaseExt(tomoFn) + "_ctf_%06d.mrc" %i
        return join(ctfDir, baseCtfFn)
    
    def _prepareCommand(self):
        tomoSet = self.inputTomograms.get()
        samRate = tomoSet.getSamplingRate()
        acquisition = tomoSet.getAcquisition()
        
        params = {'voltage': acquisition.getVoltage(),
                  'sphericalAberration': acquisition.getSphericalAberration(),
                  'magnification': acquisition.getMagnification(),
                  'ampContrast': acquisition.getAmplitudeContrast(),
                  'samplingRate': samRate * self.ctfDownFactor.get(),
                  'scannedPixelSize': tomoSet.getScannedPixelSize() * self.ctfDownFactor.get(),
                  'windowSize': self.windowSize.get(),
                  'lowRes': self.lowRes.get(),
                  'highRes': self.highRes.get(),
                  # Convert from microns to Amstrongs
                  'minDefocus': self.minDefocus.get() * 1e+4, 
                  'maxDefocus': self.maxDefocus.get() * 1e+4,
                  'astigmatism' : self.astigmatism.get(),
                  'focus' : self.stepFocus.get()
                  }
        
        # Convert digital frequencies to spatial frequencies
        params['lowRes'] = samRate / self.lowRes.get()
        if params['lowRes'] > 50:
            params['lowRes'] = 50
        params['highRes'] = samRate / self.highRes.get()
        if not self.useCftfind4:
            args, program = self._argsCtffind3()
        else:
            if self.findPhaseShift:
                params['phaseShift'] = "yes"
            else:
                params['phaseShift'] = "no"
            args, program = self._argsCtffind4()
        
        return args, program, params
    
    def _argsCtffind3(self):
        program = CTFFIND_PATH
        args = """ --old-school-input << eof > %(ctffindOut)s 
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)0.1f, %(voltage)0.1f, %(ampContrast)0.1f, %(magnification)0.1f, %(scannedPixelSize)0.1f
%(windowSize)d, %(lowRes)0.2f, %(highRes)0.2f, %(minDefocus)0.1f, %(maxDefocus)0.1f, %(focus)0.1f,  %(astigmatism)0.1f
eof
"""
        return args, program
    
    def _argsCtffind4(self):
        program = 'export OMP_NUM_THREADS=1; ' + CTFFIND4_PATH
        args = """ << eof
%(micFn)s
%(ctffindPSD)s
%(samplingRate)f
%(voltage)f
%(sphericalAberration)f
%(ampContrast)f
%(windowSize)d
%(lowRes)f
%(highRes)f
%(minDefocus)f
%(maxDefocus)f
%(focus)f
%(astigmatism)f
%(phaseShift)s
eof
"""
        return args, program
    
#     def _getCoordsFns(self):
#         """ Return a sorted list with the paths of files that matched the pattern"""
#         from glob import glob
#         coorDict = {}
#         filesPath = self.filesPath.get('').strip()
#         filesPattern = self.filesPattern.get('').strip()
#         
#         if filesPattern:
#             fullPattern = join(filesPath, filesPattern)
#         else:
#             fullPattern = filesPath
#             
#         pattern = pwutils.expandPattern(fullPattern)
#         filePaths = glob(pattern)
#         filePaths.sort()
#         for fn in filePaths:
#             coordKey = pwutils.removeBaseExt(fn)
#             coorDict[coordKey] = fn
#         return coorDict
    
#     def _getCoords(self, coordFn, coordNum):
#         crdFile = open(coordFn, 'r')
#         lines = crdFile.readlines()
#         line = lines[coordNum].split()
#         crdFile.close()
#         return line
     
#     def _getNumberOfCoords(self, tomoFn):
#         micKey = pwutils.removeBaseExt(tomoFn)
#         coordFn = self._getCoordsFns().get(micKey)
#         crdFile = open(coordFn, 'r')
#         num = len(crdFile.readlines())
#         crdFile.close()
#         return num
    
    def _getAveDefocus(self, tomoFn, micNum):
        from pyworkflow.em.packages.grigoriefflab.convert import parseCtffindOutput, parseCtffind4Output
        
        micFn = self._getImgName(tomoFn, micNum)
        ctfFn = self._getFnPath(tomoFn, '%s_ctfEstimation.txt' % pwutils.removeBaseExt(micFn))
        
        if not self.useCftfind4:
            result = parseCtffindOutput(ctfFn)
            defocusU, defocusV, _ = result
        else:
            result =  parseCtffind4Output(ctfFn)
            defocusU, defocusV, _, _, _, _ = result
        
        return (defocusU + defocusV) / 2
    
    def _getTiltAngles(self, tomoFn, micNumb):
        tiltList = self._getTiltList(tomoFn)
        return tiltList[micNumb-1]
    
    def _getTiltList(self, tomoFn):
        tiltFn = self._getFnPath(tomoFn, "tiltAngles.txt")
        tiltFile = open(tiltFn, "r")
        tiltList = tiltFile.readlines()
        tiltFile.close()
        tiltList = [float(i) for i in tiltList]
        return tiltList
    
    def _getDose(self, micNumb, numbOfMics):
        totalDose = self.inputTomograms.get().getDose()
        return totalDose * float(micNumb) / float(numbOfMics)
