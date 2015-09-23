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

#import os, sys, commands, math, stat
from pyworkflow.protocol import params
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol, ProtProcessTomograms, STEPS_PARALLEL, LEVEL_ADVANCED
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER

class ProtPrepareSubtomograms(ProtProcessTomograms):
    """sub-tomogram averaging in RELION
    """
    _label = 'subtomogram averaging'
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Import')
        form.addParam('filesPath', params.PathParam, 
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select\n"
                           "from several folders.\n\n"
                           "For example:\n"
                           "  ~/Particles/\n"
                           "  data/day??_micrographs/")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern', 
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.")
        form.addParam('copyFiles', params.BooleanParam, default=False, 
                      expertLevel=LEVEL_ADVANCED,
                      label="Copy files?",
                      help="By default the files are not copied into the\n"
                           "project to avoid data duplication and to save\n"
                           "disk space. Instead of copying, symbolic links are\n"
                           "created pointing to original files. This approach\n"
                           "has the drawback that if the project is moved to\n"
                           "another computer, the links need to be restored.\n")
        group = form.addGroup('Acquisition info')
#         form.addSection('Acquisition info')
        group.addParam('voltage', params.FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE, 
                   help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2,
                   label=Message.LABEL_SPH_ABERRATION, 
                   help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE, 
                   help=Message.TEXT_MAGNI_RATE)
        group.addParam('samplingRateMode', params.EnumParam, 
                       choices=[Message.LABEL_SAMP_MODE_1, Message.LABEL_SAMP_MODE_2],
                       default=SAMPLING_FROM_IMAGE,
                       label=Message.LABEL_SAMP_MODE,
                       help=Message.TEXT_SAMP_MODE)
        group.addParam('samplingRate', params.FloatParam,  default=4.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE, 
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('scannedPixelSize', params.FloatParam, default=12.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER,
                       label=Message.LABEL_SCANNED,
                       help='')
        group.addParam('bfactor', params.FloatParam, default=4.0,
                      label='Provide B-factor:',
                      help= '3D CTF model weighting B-factor per e-/A2')
        group.addParam('totalDose', params.FloatParam, default=40.0,
                      label='Provide acummulated dose:',
                      help= 'Total dose for the whole tomogram.')
        
        form.addSection(label=Message.LABEL_CTF_ESTI)
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

        form.addParallelSection(threads=2, mpi=3)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to preprocess the tomograms
        """
        self._insertFunctionStep("prepareDirectoriesStep")
        self._insertFunctionStep("imodExecutionStep")
        self._insertFunctionStep("estimateCtfStep")
        
        
        
    
    
    
    
    
#     
#     
#     
#     
#     
#     extracttilts -InputFile BBa.st -tilts -OutputFile tiltangles.txt > extracttilt_output.txt
#     
#     newstack -secs 0 BBa.st ctffind/BBa_image65.0_0.mrc > ctffind/temp_newstack_out.txt
#     
#     
#     
#     
#     ctffindstarname  /home/josue/PROCESAMIENTO/Tomography/tutorialData/Relion/ctffind/BBa_images.star
#     tutorialData/Relion> cat /home/josue/PROCESAMIENTO/Tomography/tutorialData/Relion/ctffind/BBa_images.star
#     data_
# 
#         loop_
#         _rlnMicrographName #1
#         ctffind/BBa_image65.0_0.mrc
#         
#     
#     
#     
#     
#     