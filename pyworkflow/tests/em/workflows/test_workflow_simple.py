# **************************************************************************
# *
# * Authors:    Jose Luis Vilas (jlvilas@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from test_workflow import TestWorkflow

from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import XmippProtPreprocessMicrographs
from pyworkflow.em.packages.simple.protocol_simple_pick import *
from pyworkflow.em.packages.simple.protocol_prime3d_initial import ProtPrime3DInitial


class TestSimpleBase(BaseTest):
    
    def runInitialModel(self, samplingRate, symmetry, 
                        maskRadius):
        #Import a set of averages
        print "Import Set of averages"
        protImportAvg = self.newProtocol(ProtImportAverages, 
                                         filesPath=self.averages, 
                                         checkStack=True, 
                                         samplingRate=samplingRate)
        self.launchProtocol(protImportAvg)
        self.assertIsNotNone(protImportAvg.getFiles(),
                             "There was a problem with the import")

     
        print "Run Simple"
        protIniModel = self.newProtocol(ProtPrime3DInitial,
                                        symmetry=symmetry,
                                        maskRadius=maskRadius,
                                        numberOfThreads=4)
        protIniModel.inputSet.set(protImportAvg.outputAverages)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVol,
                             "There was a problem with simple initial "
                             "model protocol")


class TestSimpleMDA(TestSimpleBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.averages = cls.dataset.getFile('averages')
        
    def test_mda(self):
        self.runInitialModel(3.5, 'd6', maskRadius=45)


class TestSimpleGroel(TestSimpleBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        
    def test_groel(self):
        self.runInitialModel(2.1, 'd7', maskRadius=45)


class TestWorkflowSimplePick(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def _runPickWorkflow(self):
        # First, import a set of micrographs
        print "Importing a set of micrographs..."
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.ds.getFile('micrographs'),
                                      filesPattern='*.mrc',
                                      samplingRateMode=1,
                                      magnification=79096,
                                      scannedPixelSize=56, voltage=300,
                                      sphericalAberration=2.0)
        protImport.setObjLabel('import 20 mics')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")

        print "Preprocessing the micrographs..."

        protPreprocess = self.newProtocol(XmippProtPreprocessMicrographs,
                                          doCrop=True, cropPixels=25)
        protPreprocess.inputMicrographs.set(protImport.outputMicrographs)
        protPreprocess.setObjLabel('crop 50px')
        self.launchProtocol(protPreprocess)
        self.assertIsNotNone(protPreprocess.outputMicrographs,
                             "There was a problem with the downsampling")

        print "Importing 2D averages (subset of 4)"
        ih = ImageHandler()
        classesFn = self.ds.getFile('import/classify2d/extra/'
                                    'relion_it015_classes.mrcs')

        outputName = 'input_averages.mrcs'
        inputTmp = os.path.abspath(self.proj.getTmpPath())
        outputFn = self.proj.getTmpPath(outputName)

        for i, index in enumerate([5, 16, 17, 18, 24]):
            ih.convert((index, classesFn), (i + 1, outputFn))

        protAvgs = self.newProtocol(ProtImportAverages,
                                    objLabel='avgs - 5',
                                    filesPath=inputTmp,
                                    filesPattern=outputName,
                                    samplingRate=7.08
                                    )

        self.launchProtocol(protAvgs)
        # Select some good averages from the iterations mrcs a

        protPick1 = self.newProtocol(ProtSimplePick,
                                     objLabel='simple - pick refs2d',
                                     referenceType=REF_2D,
                                     invertTemplatesContrast=False,
                                     distance=20,
                                     )

        protPick1.inputMicrographs.set(protPreprocess.outputMicrographs)
        protPick1.inputReferences.set(protAvgs.outputAverages)

        self.launchProtocol(protPick1)

        return protPick1

    def test_ribo(self):
        protPick1 = self._runPickWorkflow()
        self.assertTrue(protPick1.outputCoordinates.getSize() > 0)

        inputVol = self.ds.getFile('volumes/reference.mrc')
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         objLabel='import volume',
                                         filesPath=inputVol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)


        # Launch the same picking run but now using a volume as reference
        protPick2 = self.proj.copyProtocol(protPick1)
        protPick2.inputVolume.set(protImportVol.outputVolume)
        protPick2.setObjLabel('simple - pick refs3d')
        protPick2.referenceType.set(REF_3D)
        self.launchProtocol(protPick2)
        self.assertTrue(protPick2.outputCoordinates.getSize() > 0)
