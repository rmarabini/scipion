# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import os, sys

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtParticlePicking

import simple
from convert import readSetOfCoordinates

REF_2D = 0
REF_3D = 1


class ProtSimplePick(ProtParticlePicking):
    """
    """
    _label = 'pick'

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        em.ProtParticlePicking._defineParams(self, form)

        form.addParam('referenceType', params.EnumParam, default=REF_2D,
                      choices=['2D averages', '3D volume'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Reference type: ',
                      help='')

        form.addParam('inputReferences', params.PointerParam,
                      pointerClass="SetOfAverages",
                      condition='referenceType==%d' % REF_2D,
                      label='Input references', important=True, allowsNull=True,
                      help="Template 2D averages. ")

        form.addParam('inputVolume', params.PointerParam,
                      pointerClass="Volume",
                      condition='referenceType==%d' % REF_3D,
                      label='Input volume', important=True, allowsNull=True,
                      help="Volume to generate 2D reference templates")

        form.addParam('symmetry', params.StringParam, default='c1',
                      condition='referenceType==%d' % REF_3D,
                      label='Point group symmetry',
                      help='point-group symmetry(cn|dn|t|o|i)')

        form.addParam('invertTemplatesContrast', params.BooleanParam, default=False,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have '
                           'inverted contrast with respect to the particles '
                           'in the micrographs.\n')

        form.addParam('lowpass', params.IntParam,
                      label='Low-pass filter (A)',
                      help='')

        form.addParam('threshold', params.FloatParam, default=0.5,
                      label='Threshold',
                      help='(binarisation: 0-1; distance filer: in pixels)')

        form.addParam('rmOutliers', params.BooleanParam, default=True,
                      label='Remove outliers?',
                      help='Remove outliers')

        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional parameters')

        form.addParallelSection(mpi=1, threads=4)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('makePickReferencesStep')
        self._insertFunctionStep('convertMicrographsSteps')
        self._insertFunctionStep('runSimplePickStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ---------------------------------------------------

    def makePickReferencesStep(self):
        """ Convert input data to mrc format and generate the references.
        """
        # Convert the templates to mrcs stack or create an mrc volume
        args = ''

        if self.referenceType == REF_2D:
            inputRef = 'ref2D.mrcs'
            self.inputReferences.get().writeStack(self._getExtraPath(inputRef))
            args += 'stk=%s ' % inputRef
            args += 'pgrp=c1 '
        else:
            inputRef = 'ref3D.mrc'
            ih = em.ImageHandler()
            ih.convert(self.inputVolume.get(), self._getExtraPath(inputRef))
            args += 'vol1=%s ' % inputRef
            args += 'pgrp=%s ' % self.symmetry

        # The references are expected to be white over black micrographs
        # So, if they are inverted we don't need to do anything, if not
        # we need to pass the special flag 'neg'
        if not self.invertTemplatesContrast:
            args += 'neg=yes '

        self.runJob(simple.getProgram('makepickrefs'), args,
                    cwd=self._getExtraPath())

    def runSimplePickStep(self):
        args = 'smpd=%0.3f ' % self.getInputMicrographs().getSamplingRate()
        args += 'filetab=%s ' % self.getMicTxt()
        args += 'refs=pickrefs.mrc '
        args += 'nthr=%d ' % self.numberOfThreads
        args += 'rm_outliers=%s ' % ('yes' if self.rmOutliers else 'no')
        self.runJob(simple.getProgram('pick'), args, cwd=self._getExtraPath())

    def createOutputStep(self):
        if self.referenceType == REF_2D:
            boxSize = self.inputReferences.get().getXDim()
        else:
            boxSize = self.inputVolume.get().getXDim()

        micSet = self.getInputMicrographs()

        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(boxSize)
        readSetOfCoordinates(self.getMicrographsDir(), micSet, coordSet)
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methodsMsgs = []
        return methodsMsgs

    def _citations(self):
        return []

    # --------------------------- UTILS functions --------------------------------------------------
    def getMicrographsDir(self):
        return self._getExtraPath('micrographs')

    def getMicTxt(self):
        return "input_micrographs.txt"

    def convertMicrographsSteps(self):
        """ Write the required txt file with the list of micrographs and
        convert them to mrc format if required.
        """
        micDir = self.getMicrographsDir()  # put output and mics in extra dir
        pwutils.makePath(micDir)
        micsTxt = self._getExtraPath(self.getMicTxt())

        with open(micsTxt, 'w') as f:
            for mic in self.inputMicrographs.get():
                micPath = os.path.join(micDir, mic.getBaseName())
                micBase = os.path.relpath(micPath, self._getExtraPath())

                if micBase.endswith('.mrc'):
                    pwutils.createLink(mic.getFileName(), micPath)
                else:
                    em.ImageHandler().convert(mic, micPath)
                # Add the name in the list relative to extra path
                f.write('%s\n' % micBase)

