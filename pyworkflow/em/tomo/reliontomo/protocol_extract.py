# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
#  
import pyworkflow.utils.path as putils
from pyworkflow.em import RELATION_CTF, Subtomogram
from pyworkflow.protocol import params
from pyworkflow.em.protocol import ProtExtractSubtomograms

import convert as conv

class ProtRelionExtractSubtomograms(ProtExtractSubtomograms):
    """ Wrapper to Relion preprocess program.
    This protocol provides an easy way to execute *relion_preprocess* program
    to extract subtomograms.
    """
    _label = 'extract subtomograms'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Coordinates", important=True,
                      pointerClass='SetOfTomoCoordinates',
                      help='Select the SetOfCoordinates ')
        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
                      relationName=RELATION_CTF, attributeName='getCoordinates',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n'
                           'CTF estimation is needed if you want to do phase flipping or \n'
                           'you want to associate CTF information to the particles.')
        
        form.addSection(label='Preprocess')
        form.addParam('doInvert', params.BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize', important=True,
                      help='If set to True, particles will be normalized in the way RELION prefers it.\n'
                           'It is recommended to *always normalize your particles*, and use a reasonable\n'
                           'radius for the circle around your particles outside of which the standard\n'
                           'deviation and average values for the noise are calculated.\n\n'
                           '*Note*: if the particles are re-scaled, the radius for normalize will be\n'
                           'taken over the new dimensions.')
        form.addParam('backRadius', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius (px)',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pixel).')
        line = form.addLine('Dust removal',
                            help='Remove very white/dark pixels from the extracted particles.'
                                 ' Pixels values higher/lower than this many times the image stddev'
                                 ' will be replaced with values from a Gaussian distribution.'
                                 ' Use negative value to switch off dust removal.')
        line.addParam('whiteDust', params.FloatParam, default=-1., label='White')
        line.addParam('blackDust', params.FloatParam, default=-1., label='Black')
        
        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.coordSet = self.inputCoordinates.get()
        self.starTomoFn = self._getPath('all_tomograms.star')
        self._insertFunctionStep("convertInputStep", self.coordSet.getObjId())
        self._insertFunctionStep('extractStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, coordsId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        recSet = self.coordSet.getTomoRecs()
        conv.writeSetOfTomograms(recSet, self.starTomoFn, self._getExtraPath())
        self._writeSetOfTomoCoords()
    
    def extractStep(self):
        args = " --o subtomo --mic_star %(star)s --coord_suffix .coord --extract --extract_size %(boxSize)d "
        paramsDict = {"star" : os.path.basename(self.starTomoFn),
                      "boxSize" : self.coordSet.getBoxSize()
                      }
        params = args % paramsDict
        
        if self.doNormalize:
            radius = self.backRadius.get()
            if radius <= 0:
                radius = 0
            params += ' --norm --bg_radius %d' % radius
             
        if self.doInvert:
            params += ' --invert_contrast'            
        
        wDust = self.whiteDust.get()
        bDust = self.blackDust.get()
        
        if wDust > 0:
            params += ' --white_dust %f' % wDust
        if bDust > 0:
            params += ' --black_dust %f' % bDust
 
        self.runJob(self._getProgram(), params, env=conv.getEnviron(), cwd=self._getPath())
     
    def createOutputStep(self):
        subtomoSet = self._createSetOfSubtomograms()
        subtomoSet.copyInfo(self.coordSet.getTomoRecs())
        subtomoSet.setHasCTF(self.ctfRelations.hasValue())
        self.coordDict = {}
        
        for coord in self.coordSet:
            subtomoFn = self._getPath(self._getSubtomoFn(coord))
            if subtomoFn is not None:
                subtomo = Subtomogram()
                subtomo.setFileName(subtomoFn)
                subtomo.copyInfo(coord.getTomoRec())
                subtomo.setCoordinate(coord)
                subtomo.setTomoId(coord.getTomoId())
                subtomo.setTomoName(coord.getTomoName())
                if self.ctfRelations:
                    subtomo.setCTF(self.ctfRelations.get()[coord.getObjId()])
                subtomoSet.append(subtomo)
        self._defineOutputs(outputSubtomograms=subtomoSet)
        self._defineSourceRelation(self.coordSet, subtomoSet)
    
#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        validateMsgs = []
        return validateMsgs
     
    def _summary(self):
        summary = []
 
        return summary
     
    def _methods(self):
        return []

    #--------------------------- UTILS functions ---------------------------------------------------
    def getCoordinates(self):
        return self.inputCoordinates.get()
    
    def _getProgram(self, program='relion_preprocess'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program
    
    def _writeSetOfTomoCoords(self):
        # Create a dictionary with the coords filenames for each tomogram
        posDict = {}
        for tomo in self.coordSet.iterTomoRecs():
            tomoName = tomo.getFileName()
            posFn = os.path.join(self._getExtraPath(), putils.replaceBaseExt(tomoName, "coord"))
            posDict[tomo.getObjId()] = posFn
        
        f = None
        lastTomoId = None
        c = 0
        
        for coord in self.coordSet.iterItems(orderBy='_tomoId'):
            tomoId = coord.getTomoId()
        
            if tomoId != lastTomoId:
                # we need to close previous opened file
                if f:
                    f.close()
                    print "Micrograph %s (%d)" % (lastTomoId, c)
                    c = 0
                f = open(posDict[tomoId], 'w+')
                lastTomoId = tomoId
            c += 1
            f.write(" %d   %d   %d\n" % (coord.getX(), coord.getY(), coord.getZ()))
        
        if f:
            f.close()
            print "Micrograph %s (%d)" % (lastTomoId, c)


    def _getSubtomoFn(self, coord):
        import pyworkflow.em.metadata as md
        if self.coordDict == {}:
            subTomoMd = md.MetaData(self._getPath("subtomo.star"))
            for objId in subTomoMd:
                tomoBsName = putils.removeBaseExt(subTomoMd.getValue(md.RLN_MICROGRAPH_NAME, objId))
                x = str(int(subTomoMd.getValue(md.RLN_IMAGE_COORD_X, objId)))
                y = str(int(subTomoMd.getValue(md.RLN_IMAGE_COORD_Y, objId)))
                z = str(int(subTomoMd.getValue(md.RLN_IMAGE_COORD_Z, objId)))
                subtomoFn = subTomoMd.getValue(md.RLN_IMAGE_NAME, objId)
                key = tomoBsName + x + y + z
                self.coordDict[key] = subtomoFn
        coordKey = putils.removeBaseExt(coord.getTomoName()) + str(coord.getX()) + str(coord.getY()) + str(coord.getZ())
        return self.coordDict[coordKey] if self.coordDict[coordKey] else None
        