# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

import sys, unittest

from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.tomo import *


# Some utility functions to import Tomograms that are used
# in several tests.
class TestTomogramBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='Tomo'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.tomogramsFn = cls.dataset.getFile('allTomograms')
        cls.tomoRecsFn = cls.dataset.getFile('allRec')
        cls.coords = cls.dataset.getFile('allCoordinates')
    
    @classmethod
    def runImportTomograms(cls, pattern, vol, spAbrr, samplingRate):
        """ Run an Import tomograms protocol. """
        protImport = cls.newProtocol(ProtImportTomograms,
                                     voltage=vol,
                                     sphericalAberration=spAbrr,
                                     filesPath=pattern,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputTomograms is None:
            raise Exception('Import of images: %s, failed. outputTomograms is None.' % pattern)
        return protImport
    
    @classmethod
    def runImportTomoRecs(cls, pattern, tomograms, sampling):
        """ Run an Import reconstructed tomograms protocol. """
        protImport2 = cls.newProtocol(ProtImportTomoRecs,
                                     filesPath=pattern,
                                     samplingRate=sampling,
                                     inputTomograms=tomograms)
        cls.launchProtocol(protImport2)
        # check that input images have been imported (a better way to do this?)
        if protImport2.outputTomoRecs is None:
            raise Exception('Import of images: %s, failed. outputTomograms is None.' % pattern)
        return protImport2
    
    @classmethod
    def runImportTomoCoords(cls, pattern, tomoRecs, bSize):
        """ Run an Import reconstructed tomograms protocol. """
        protImportCrds = cls.newProtocol(ProtImportTomoCoordinates,
                                         filesPath=pattern,
                                         inputTomoRecs=tomoRecs,
                                         boxSize=bSize)
        cls.launchProtocol(protImportCrds)
        # check that input images have been imported (a better way to do this?)
        if protImportCrds.outputTomoCoordinates is None:
            raise Exception('Import of images: %s, failed. outputTomograms is None.' % pattern)
        return protImportCrds

    


class TestCtf3DEstimation(TestTomogramBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestTomogramBase.setData()
        cls.protImport = cls.runImportTomograms(cls.tomogramsFn, 300, 2.26, 2.6)
        cls.protImportRec = cls.runImportTomoRecs(cls.tomoRecsFn, cls.protImport.outputTomograms, 10.4)
        cls.protImportCrds = cls.runImportTomoCoords(cls.coords, cls.protImportRec.outputTomoRecs, 50)
    
    def testCtf3DEstimation(self):
        print "Run estimate CTF3Ds"
        protCtf3D = self.newProtocol(ProtCtf3DEstimation,
                                       ctfDownFactor=2,
                                       lowRes=0.002,
                                       minDefocus=1.0,
                                       maxDefocus=10,
                                       stepFocus=2000,
                                       numberOfThreads=5)
        protCtf3D.inputTomoCoords.set(self.protImportCrds.outputTomoCoordinates)
        self.launchProtocol(protCtf3D)
#          
#         self.assertIsNotNone(getattr(prot2D, 'outputClasses', None), 
#                              "There was a problem with Relion 2D:\n" + prot2D.getErrorMessage())
#         self.assertAlmostEquals(self.protNormalize.outputParticles.getSamplingRate(), 
#                                 prot2D.outputClasses.getImages().getSamplingRate(),
#                                 "There was a problem with the sampling rate of the particles")

#         print "Run extracting subtomograms"
#         protSub = self.newProtocol(ProtRelionExtractSubtomograms,
#                                    backRadius=20)
#         protSub.inputCoordinates.set(self.protImportCrds.outputTomoCoordinates)
#         protSub.ctfRelations.set(protCtf3D.outpuCft3Ds)
#         self.launchProtocol(protSub)
        
        
        
