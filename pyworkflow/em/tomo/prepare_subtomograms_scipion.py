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
        form.addParam('inputTomograms', params.PointerParam, important=True,
                      label=Message.LABEL_INPUT_TOM,
                       pointerClass='SetOfTomograms')
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
        form.addParam('astigmatism', params.FloatParam, default=100.0,
              label='Expected (tolerated) astigmatism', expertLevel=params.LEVEL_ADVANCED,
              condition='useCftfind4', )
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
        line.addParam('maxDefocus', params.FloatParam, default=4.,
                      label='Max')
        
        form.addParam('windowSize', params.IntParam, default=256, expertLevel=LEVEL_ADVANCED,
                      label='Window size',
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier.')
        
        form.addParallelSection(threads=1, mpi=3)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to preprocess the tomograms
        """
        ctfDepsList = []
        tomoSet = self.inputTomograms.get()
        for tomo in tomoSet:
            extractDeps = self._insertFunctionStep("extractTiltAnglesStep", tomo.getFileName(), prerequisites=ctfDepsList)
            for i in range(1, tomo.getDim()+1):
                tiltImgDeps = self._insertFunctionStep("extractTiltImgStep", tomo.getFileName(), i, prerequisites=[extractDeps])
                ctfDeps = self._insertFunctionStep("estimateCtfStep", tomo.getFileName(), i, prerequisites=[tiltImgDeps])
                ctfDepsList.append(ctfDeps)
    
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
    
    
    def _getTomoPath(self, tomoFn):
        tomoBaseDir = pwutils.removeBaseExt(tomoFn)
        pwutils.makePath(self._getExtraPath(tomoBaseDir))
        return self._getExtraPath(tomoBaseDir)
    
    def _getFnPath(self, tomoFn, fn):
        return join(self._getTomoPath(tomoFn), basename(fn))
    
    def _getImgName(self, tomoFn, i):
        baseImgFn = pwutils.removeBaseExt(tomoFn) + "_%03d.mrc" %i
        return join(self._getTomoPath(tomoFn), baseImgFn)

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
                  'maxDefocus': self.maxDefocus.get() * 1e+4
                  }
        
        # Convert digital frequencies to spatial frequencies
        params['lowRes'] = samRate / self.lowRes.get()
        if params['lowRes'] > 50:
            params['lowRes'] = 50
        params['highRes'] = samRate / self.highRes.get()
        params['step_focus'] = 500.0
        if not self.useCftfind4:
            args, program = self._argsCtffind3()
        else:
            params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                params['phaseShift'] = "yes"
            else:
                params['phaseShift'] = "no"
            args, program = self._argsCtffind4()
        
        return args, program, params
    
    def _argsCtffind3(self):
        program = 'export NATIVEMTZ=kk ; ' + CTFFIND_PATH
        args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
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
%(step_focus)f
%(astigmatism)f
%(phaseShift)s
eof
"""            
        return args, program
