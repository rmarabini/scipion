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
        cls.coords = cls.dataset.getFile('allCoordinates')
    
    @classmethod
    def runImportTomograms(cls, pattern, vol, spAbrr, samplingRate):
        """ Run an Import particles protocol. """
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


class TestSubtomogramAveraging(TestTomogramBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestTomogramBase.setData()
        cls.protImport = cls.runImportTomograms(cls.tomogramsFn, 300, 2.26, 2.6)
    
    def testSubtomogramAverage(self):
        print "Run prepare tomograms"
        protTomoAve = self.newProtocol(ProtPrepareSubtomograms,
                                       ctfDownFactor=2,
                                       lowRes=0.002,
                                       minDefocus=1.0,
                                       maxDefocus=10,
                                       stepFocus=2000,
                                       filesPath=self.coords,
                                       boxSize=200,
                                       numberOfThreads=5)
        protTomoAve.inputTomograms.set(self.protImport.outputTomograms)
        self.launchProtocol(protTomoAve)
#          
#         self.assertIsNotNone(getattr(prot2D, 'outputClasses', None), 
#                              "There was a problem with Relion 2D:\n" + prot2D.getErrorMessage())
#         self.assertAlmostEquals(self.protNormalize.outputParticles.getSamplingRate(), 
#                                 prot2D.outputClasses.getImages().getSamplingRate(),
#                                 "There was a problem with the sampling rate of the particles")


