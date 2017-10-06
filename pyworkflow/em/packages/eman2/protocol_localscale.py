# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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

import os

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import IntParam, FloatParam, PointerParam
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em as em

import eman2
from convert import readSetOfCoordinates



class EmanProtLocScale(em.ProtRefine3D):
    """
    Protocol to pick particles automatically in a set of micrographs
    using sparx gaussian picker.
    For more information see http://sparx-em.org/sparxwiki/e2boxer
    """
    _label = 'local scale'
    # _lastUpdateVersion = VERSION_1_1
        
    # def __init__(self, **args):     
    #     ProtParticlePicking.__init__(self, **args)
    #     self.extraParams = 'pixel_input=1:pixel_output=1:invert_contrast=True:use_variance=True'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label='Input volume', important=True,
                      help='Volume to process.')
        form.addParam('modelVolume', PointerParam, pointerClass='Volume',
                      label='Model volume', important=True,
                      help='Model volume to follow.')
        form.addParam('mask', PointerParam, label="Mask", 
                      pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be binary: 0 (remove these voxels)'
                           'and 1 (let them pass).')

        form.addParallelSection(threads=0,mpi=4)

    # def _createIterTemplates(self, currRun):
    #     """ Setup the regex on how to find iterations. """
    #     self._iterTemplate = self._getFileName('mapFull', run=currRun,
    #                                            iter=1).replace('threed_01',
    #                                                            'threed_??')
    #     # Iterations will be identify by threed_XX_ where XX is the iteration
    #     #  number and is restricted to only 2 digits.
    #     self._iterRegex = re.compile('threed_(\d{2,2})')

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):        
        # self._createFilenameTemplates()
        # self._createIterTemplates(self._getRun())
        self._insertFunctionStep('locScaleStep')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions -------------------------------
    def locScaleStep(self):
        args = self._prepareParams()
        print("localscalaBIN.py" + args)
        # self.runJob('e2locScale.py', args, cwd=self.getWorkingDir()) 
        
    def createOutputStep(self):
        print("inside of createOutputStep")
        # coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        # self.readSetOfCoordinates(self.workingDir.get(), coordSet)
        # coordSet.setBoxSize(self.boxSize.get())
        # self._defineOutputs(outputCoordinates=coordSet)
        # self._defineSourceRelation(self.inputMicrographs, coordSet)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        eman2.validateVersion(self, errors)
        return errors

    #--------------------------- UTILS functions -------------------------------
    def _prepareParams(self):
        volume = self.inputVolume.get()
        volumeFn = volume.getFileName()
        print(volume)
        samRate = volume.getSamplingRate()
        Xdim = volume.getDim()[0]

        modelFn = self.modelVolume.get().getFileName()

        args = "--em_map %s --model_map %s --apix %f"%(volumeFn,modelFn,samRate) 

        fnMask = ""
        if self.mask.hasValue():
            fnMask = self._getExtraPath("mask.vol")
            img=ImageHandler()
            img.convert(self.mask.get(), fnMask)
            self.runJob('xmipp_image_resize',"-i %s --dim %d"%(fnMask,Xdim),
                        numberOfMpi=1)
            self.runJob('xmipp_transform_threshold',
                        "-i %s --select below 0.5 --substitute binarize"%fnMask,
                        numberOfMpi=1)
            args+=" --mask binary_file %s" % fnMask

        args += " --outfile %s" % self._getExtraPath()
        args += " --mpi %d" % self.numberOfMpi
        print(args)
        return args

        # cmdl_parser.add_argument('-em', '--em_map', required=True, help='Input filename EM map')
        # cmdl_parser.add_argument('-mm', '--model_map', required=True, help='Input filename PDB map')
        # cmdl_parser.add_argument('-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
        # cmdl_parser.add_argument('-ma', '--mask', help='Input filename mask')
        # cmdl_parser.add_argument('-w', '--window_size', type=int, help='window size in pixel')
        # cmdl_parser.add_argument('-o', '--outfile', required=True, help='Output filename')
        # cmdl_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,
        #                  help='MPI version call by: \"{0}\"'.format(mpi_cmd))