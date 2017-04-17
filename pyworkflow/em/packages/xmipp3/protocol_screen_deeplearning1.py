# **************************************************************************
# *
# * Authors:  Ruben Sanchez (rsanchez@cnb.csic.es), April 2017
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


from os.path import basename
import numpy as np
from math import ceil

from pyworkflow.utils.path import copyFile, cleanPath 
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtProcessParticles

import pyworkflow.em.metadata as md
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, setXmippAttributes
import xmipp

class DeepTFSupervised:
    def read(self, fnNet):
        pass

    def write(self, fnNet):
        pass
    
    def predict(self, I):
        score = np.random.rand()
        return score

class XmippProtScreenDeepLearning1(ProtProcessParticles):
    """ Protocol for screening particles using deep learning. """
    _label = 'screen deep learning 1'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTrueParticles', params.PointerParam, label="Set of true particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of true particles.')  
        form.addParam('inputNegatives', params.PointerParam, label="Set of negative examples", 
                      pointerClass='SetOfParticles',
                      help='Select the set of non-particles.')  
        form.addParam('inputParticles', params.PointerParam, label="Set of putative particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of possible particles.')
        form.addParam('Nepochs', params.FloatParam, label="Number of epochs", default=3, expertLevel=params.LEVEL_ADVANCED,
                      help='Number of epochs for the neural network training')  
        form.addParam('learningRate', params.FloatParam, label="Learning rate", default=0.001, expertLevel=params.LEVEL_ADVANCED,
                      help='Learning rate of the neural network training')  
                      
        form.addParallelSection(threads=8, mpi=0)    
    
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.inputTrueParticles.get(), self.inputNegatives.get(), self.inputParticles.get()) 
        self._insertFunctionStep('train',self.inputTrueParticles.get(), self.inputNegatives.get(), self.learningRate.get()) 
        self._insertFunctionStep('predict',self.inputParticles.get()) 
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------   
    def convertInputStep(self, inputTrueParticles, inputNegatives, inputParticles):
        writeSetOfParticles(inputTrueParticles, self._getExtraPath("trueParticles.xmd"))
        writeSetOfParticles(inputNegatives, self._getExtraPath("negatives.xmd"))
        writeSetOfParticles(inputParticles, self._getPath("particles.xmd")) # Particles to classify
    
    def getBatch(self, mdIn, stack, batchSize):
        mdAux = md.MetaData()
        mdAux.randomize(mdIn)
        n = 0
        I = xmipp.Image()
        for id in mdAux:
            fnImage = mdAux.getValue(md.MDL_IMAGE, id)
            I.read(fnImage)
            stack[...,n]=I.getData()
            n+=1
            if n>=batchSize:
                break
    
    def train(self, inputTrueParticles, inputNegatives, learningRate):
        batchSize=100
        mdTrue = md.MetaData(self._getExtraPath("trueParticles.xmd"))
        mdNegatives = md.MetaData(self._getExtraPath("negatives.xmd"))
        xdim, ydim, _ = inputTrueParticles.getDim()
        ndim = inputTrueParticles.getSize()
        trueStack = np.zeros((xdim,ydim,batchSize))
        negativeStack = np.zeros((xdim,ydim,batchSize))
        batchNumber = int(ceil(ndim*self.Nepochs.get()/batchSize))
        nnet = DeepTFSupervised()
        for batch in range(batchNumber):
            self.getBatch(mdTrue, trueStack, batchSize) 
            self.getBatch(mdNegatives, negativeStack, batchSize)
            nnet.write(self._getExtraPath("nnet.tf"))

    def predict(self, inputParticles):
        nnet = DeepTFSupervised()
        nnet.read(self._getExtraPath("nnet.tf"))
        mdParticles = md.MetaData(self._getPath("particles.xmd"))
        I = xmipp.Image()
        for id in mdParticles:
            fnImage = mdParticles.getValue(md.MDL_IMAGE, id)
            I.read(fnImage)
            score = nnet.predict(I.getData())
            mdParticles.setValue(md.MDL_ZSCORE_DEEPLEARNING1, score, id)
        mdParticles.write(self._getPath("particles.xmd"))

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(imgSet)
        partSet.copyItems(imgSet,
                            updateItemCallback=self._updateParticle,
                            itemDataIterator=md.iterRows(self._getPath("particles.xmd"), sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        pass
    
    #--------------------------- UTILS functions --------------------------------------------
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_ZSCORE_DEEPLEARNING1)
        if row.getValue(md.MDL_ENABLED) <= 0:
            item._appendItem = False
        else:
            item._appendItem = True