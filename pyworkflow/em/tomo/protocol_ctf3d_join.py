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

from pyworkflow.protocol import params
from pyworkflow.em.protocol import ProtProcessTomograms, STEPS_PARALLEL


class ProtCtf3DJoin(ProtProcessTomograms):
    """ Join Ctf3D and Coordinates coming from independent tomographs.

    NOTE:
    This is a workaround to do tomo processing in Frank's lab where we can use
    threads and process tomographs separately.
    """
    _label = 'ctf3D join'
    
    def __init__(self, **kwargs):
        ProtProcessTomograms.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCtfs', params.MultiPointerParam,
                      pointerClass='SetOfCTF3D',
                      label="Input CTFs 3D", important=True)
        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep("createOutputStep")
    
    #--------------------------- STEPS functions ---------------------------------------------------

    def createOutputStep(self):
        outTomoRecSet = self._createSetOfTomoRecs()
        outCoordSet = self._createSetOfTomoCoordinates(outTomoRecSet)
        outCtf3DSet = self._createSetOfCTF3D(outCoordSet)

        firstCtf3DSet = self.inputCtfs[0].get()
        outCtf3DSet.copyInfo(firstCtf3DSet)
        firstCoordSet = firstCtf3DSet.getTomoCoordinates()
        outCoordSet.copyInfo(firstCoordSet)
        firstTomoRecSet = firstCoordSet.getTomoRecs()
        outTomoRecSet.copyInfo(firstTomoRecSet)
        
        for ctfPointer in self.inputCtfs:
            ctf3DSet = ctfPointer.get()
            
            firstCtf3D = ctf3DSet.getFirstItem()
            firstCoord = firstCtf3D.getTomoCoordinate()
            tomoRec = firstCoord.getTomoRec()
            tomoRec.setObjId(None)
            
            for ctf3D in ctf3DSet:
                ctf3D.setObjId(None)
                coord = ctf3D.getTomoCoordinate()
                coord.copyObjId(ctf3D)
                coord.setTomoRec(tomoRec)
                ctf3D.setTomoCoordinate(coord)
                outCtf3DSet.append(ctf3D)
                outCoordSet.append(coord)
            
            outTomoRecSet.append(tomoRec)
        
        self._defineOutputs(outputTomoRecs=outTomoRecSet)
        self._defineOutputs(outputTomoCoordinates=outCoordSet)
        self._defineOutputs(outputCft3Ds=outCtf3DSet)
        self._defineCtfRelation(outCoordSet, outCtf3DSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors    
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        methods = ''
        return [methods]
        
    def _citations(self):
        pass
    
    #--------------------------- UTILS functions --------------------------------------------
