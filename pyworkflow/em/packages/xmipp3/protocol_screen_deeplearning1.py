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


import os

from pyworkflow.utils.path import copyFile, cleanPath 
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtProcessParticles
import pyworkflow.em.metadata as md
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, setXmippAttributes


class XmippProtScreenDeepLearning1(ProtProcessParticles):
    """ Protocol for screening particles using deep learning. """
    _label = 'screen deep learning 1'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label='Continue training from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                      'run of type *%s* class and some of the input parameters'
                      'will be taken from it.' % self.getClassName())
        form.addParam('continueRun', params.PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')

        form.addParam('inPosSetOfParticles', params.PointerParam, label="Consensus particles (mostly true particles)", 
                      pointerClass='SetOfParticles', 
                      help='Select the intersection set of particles (mostly true particles).')  
        form.addParam('inNegSetOfParticles', params.PointerParam, label="Set of negative particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of non-particles.')

        form.addParam('predictSetOfParticles', params.PointerParam, label="Set of putative particles to predict", 
                      pointerClass='SetOfParticles',
                      help='Select the set of putative particles particles to classify.')

        form.addParam('auto_stopping',params.BooleanParam, default=True,
                      label='Auto stop training when convergency is detected?',
                      help='If you set to *Yes*, the program will automatically stop training'
                      'if there is no improvement for consecutive 2 epochs, learning rate will be decreased'
					  'by a factor 10. If learningRate_t < 0.01*learningrate_0 training will stop. Warning: '
                      'Sometimes convergency seems to be reached, but after time, improvement can still happen. '
                      'Not recommended for very small data sets (<100 true particles)')

        if 'CUDA' in os.environ and not os.environ['CUDA']=="False":

            form.addParam('gpuToUse', params.IntParam, default='0',
                          label='Which GPU to use:',
                          help='Currently just one GPU will be use, by '
			                   'default GPU number 0 You can override the default '
                               'allocation by providing other GPU number, p.e. 2')
        else:
            form.addParallelSection(threads=8, mpi=0)

        form.addParam('nEpochs', params.FloatParam, label="Number of epochs", default=10.0, expertLevel=params.LEVEL_ADVANCED,
                      help='Number of epochs for neural network training.')  
        form.addParam('learningRate', params.FloatParam, label="Learning rate", default=1e-3, expertLevel=params.LEVEL_ADVANCED,
                      help='Learning rate for neural network training')
        form.addParam('nModels', params.IntParam, label="Number of models for ensemble", default=3, expertLevel=params.LEVEL_ADVANCED,
                      help='Number of models to fit in order to build an ensamble. Tipical values are 3 to 20. The more the better'
                      'until a point where no gain is obtained')  

        form.addSection(label='testingData')
        form.addParam('doTesting', params.BooleanParam, default=False,
                      label='Perform testing during training?',
                      help='If you set to *Yes*, you should select a testing positive set '
                      'and a testing negative set')
        form.addParam('testPosSetOfParticles', params.PointerParam, condition='doTesting',
                      label="Set of positive test particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of ground true positive particles.')
        form.addParam('testNegSetOfParticles', params.PointerParam, condition='doTesting',
                      label="Set of negative test particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of ground false positive particles.')

    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.inPosSetOfParticles.get(), self.inNegSetOfParticles.get(),
				                     self.predictSetOfParticles.get(), 
				                     self.testPosSetOfParticles.get(), self.testNegSetOfParticles.get()) 
        self._insertFunctionStep('train',self.inPosSetOfParticles.get(), self.inNegSetOfParticles.get(), self.testPosSetOfParticles.get(),
				                     self.testNegSetOfParticles.get(), self.learningRate.get()) 
        self._insertFunctionStep('predict',self.testPosSetOfParticles.get(),self.testNegSetOfParticles.get(),
				                     self.predictSetOfParticles.get(), self.inPosSetOfParticles.get(), self.inNegSetOfParticles.get())
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------   
    def convertInputStep(self, inPosSetOfParticles, inNegSetOfParticles, predictSetOfParticles, testPosSetOfParticles,testNegSetOfParticles):
        writeSetOfParticles(inPosSetOfParticles, self._getExtraPath("inputTrueParticlesSet.xmd"))
        writeSetOfParticles(inNegSetOfParticles, self._getExtraPath("inputFalseParticlesSet.xmd"))
        writeSetOfParticles(predictSetOfParticles, self._getExtraPath("predictSetOfParticles.xmd"))
        if not testPosSetOfParticles is None and not testNegSetOfParticles is None:
          writeSetOfParticles(testPosSetOfParticles, self._getExtraPath("testTrueParticlesSet.xmd"))
          writeSetOfParticles(testNegSetOfParticles, self._getExtraPath("testFalseParticlesSet.xmd"))
        
              
    def train(self, inPosSetOfParticles, inNegSetOfParticles, testPosSetOfParticles, testNegSetOfParticles, learningRate):
        '''
        inPosSetOfParticles, inNegSetOfParticles, testPosSetOfParticles, testNegSetOfParticles: SetOfParticles
        learningRate: float
        '''
        if self.doContinue.get():
            previousRun=self.continueRun.get()
            print( previousRun._getExtraPath('nnetData') )
            copyDir(previousRun._getExtraPath('nnetData'),self._getExtraPath('nnetData'))

      
        from pyworkflow.em.packages.xmipp3.deepLearning1 import  DeepTFSupervised, DataManager, updateEnviron

        if self.gpuToUse:
            updateEnviron( self.gpuToUse.get() )
            numberOfThreads=None
        else:
            numberOfThreads=self.numberOfThreads.get()

        trainDataManager= DataManager(posImagesXMDFname= self._getExtraPath("inputTrueParticlesSet.xmd"),
                                       posImagesSetOfParticles= inPosSetOfParticles,
                                       negImagesXMDFname= self._getExtraPath("inputFalseParticlesSet.xmd"),
                                       negImagesSetOfParticles= inNegSetOfParticles)

        if not testPosSetOfParticles is None and not testNegSetOfParticles is None:
			testDataManager= DataManager(posImagesXMDFname=  self._getExtraPath("testTrueParticlesSet.xmd"),
							    posImagesSetOfParticles= testPosSetOfParticles,
							    negImagesXMDFname= self._getExtraPath("testFalseParticlesSet.xmd"),
							    negImagesSetOfParticles= testNegSetOfParticles)
        else:
			testDataManager= None
        nEpochs= self.nEpochs.get()
        numberOfBatches = trainDataManager.getNBatches(nEpochs)
        numModels= self.nModels.get()
        assert numModels>=1, "Error, nModels>=1"
        for i in range(numModels):
            print("Training model %d/%d"%((i+1), numModels))
            nnet = DeepTFSupervised(rootPath=self._getExtraPath("nnetData"), modelNum=i)        
            nnet.createNet( trainDataManager.shape[0], trainDataManager.shape[1], trainDataManager.shape[2], trainDataManager.nTrue)
            nnet.startSessionAndInitialize(numberOfThreads)
            nnet.trainNet(numberOfBatches, trainDataManager, learningRate, testDataManager, self.auto_stopping.get())
            nnet.close(saveModel= False) #Models will be automatically saved during training, so True no needed
    #        self.predict( testPosSetOfParticles, testNegSetOfParticles, inPosSetOfParticles, inNegSetOfParticles)
    #        raise ValueError("Debug mode")
            del nnet


    def predict(self, testPosSetOfParticles, testNegSetOfParticles, predictSetOfParticles, inPosSetOfParticles, inNegSetOfParticles):
        from pyworkflow.em.packages.xmipp3.deepLearning1 import  DeepTFSupervised, DataManager, updateEnviron
        import numpy as np
        from sklearn.metrics import accuracy_score, roc_auc_score

        if self.gpuToUse:
            updateEnviron( self.gpuToUse.get() )
            numberOfThreads=None
        else:
            numberOfThreads=self.numberOfThreads.get()

        trainDataManager= DataManager(posImagesXMDFname= self._getExtraPath("inputTrueParticlesSet.xmd"),
                                       posImagesSetOfParticles= inPosSetOfParticles,
                                       negImagesXMDFname= self._getExtraPath("inputFalseParticlesSet.xmd"),
                                       negImagesSetOfParticles= inNegSetOfParticles)

        predictDataManager= DataManager(posImagesXMDFname= self._getExtraPath("predictSetOfParticles.xmd"),
                                       posImagesSetOfParticles= predictSetOfParticles,
                                       negImagesXMDFname= None,
                                       negImagesSetOfParticles= None)
        numModels= self.nModels.get()

        resultsDictPos={}
        resultsDictNeg={}
        for i in range(numModels):
            print("Predicting with model %d/%d"%((i+1), numModels))

            nnet = DeepTFSupervised(rootPath=self._getExtraPath("nnetData"), modelNum=i)
            nnet.createNet( trainDataManager.shape[0], trainDataManager.shape[1], trainDataManager.shape[2], trainDataManager.nTrue)
            nnet.startSessionAndInitialize(numberOfThreads)
            y_pred, labels, typeAndIdList = nnet.predictNet(predictDataManager)
            nnet.close(saveModel= False)

            for score, label, (mdIsPosType, mdId) in zip(y_pred , labels, typeAndIdList):
              if mdIsPosType==True:
                 try:
                     resultsDictPos[mdId]+= float(score)/float(numModels)
                 except KeyError:
                     resultsDictPos[mdId]=0.0
              else:
                 try:
                     resultsDictNeg[mdId]+= float(score)/float(numModels)
                 except KeyError:
                     resultsDictNeg[mdId]=0.0
        metadataPos, metadataNeg= predictDataManager.getMetadata()
        for mdId in resultsDictPos:
            metadataPos.setValue(md.MDL_ZSCORE_DEEPLEARNING1, resultsDictPos[mdId], mdId)

        for mdId in resultsDictNeg:
            metadataNeg.setValue(md.MDL_ZSCORE_DEEPLEARNING1, resultsDictNeg[mdId], mdId)

        metadataPos.write(self._getPath("particles.xmd"))

        if not testPosSetOfParticles is None and not testNegSetOfParticles is None:
          testDataManager= DataManager(posImagesXMDFname=  self._getExtraPath("testTrueParticlesSet.xmd"), 
                                posImagesSetOfParticles= testPosSetOfParticles,
                                negImagesXMDFname= self._getExtraPath("testFalseParticlesSet.xmd"),
                                negImagesSetOfParticles= testNegSetOfParticles)

          nnet.close(saveModel= False)
          scores_list=[]
          labels_list=[]
          cum_acc_list=[]
          cum_auc_list=[]
          for i in range(numModels):
            print("Predicting test data with model %d/%d"%((i+1), numModels))

            nnet = DeepTFSupervised(rootPath=self._getExtraPath("nnetData"), modelNum=i)
            nnet.createNet( trainDataManager.shape[0], trainDataManager.shape[1], trainDataManager.shape[2], trainDataManager.nTrue)
            nnet.startSessionAndInitialize(numberOfThreads)
            y_pred, labels, typeAndIdList = nnet.predictNet(testDataManager)
            scores_list.append(y_pred)
            labels= [ 0 if label[0]==1.0 else 1 for label in labels]
            labels_list.append(labels)
            curr_auc= roc_auc_score(labels, y_pred)
            curr_acc= accuracy_score(labels, [1 if y>0.5 else 0 for y in y_pred])
            cum_acc_list.append(curr_acc)
            cum_auc_list.append(curr_auc)
            print("Model %d test accuracy: %f  auc: %f"%(i, curr_acc, curr_auc))
            nnet.close(saveModel= False)

          labels= np.concatenate( labels_list)
          scores= np.concatenate( scores_list)
          auc_val= roc_auc_score(labels, scores)
          scores[ scores>0.5]=1
          scores[ scores<=0.5]=0
          accuracy_score(labels, scores)
          print(">>>>>>>>>>>>\nEnsemble test accuracy: %f  auc: %f"%(accuracy_score(labels, scores) , auc_val))
          print("Mean single model test accuracy: %f  auc: %f"%(np.mean(cum_acc_list) , np.mean(cum_auc_list)))

    def createOutputStep(self):
        imgSet = self.predictSetOfParticles.get()
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
