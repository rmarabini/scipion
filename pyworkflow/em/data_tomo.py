# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
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
This modules contains data classes related to Tomography workflow.
"""

#NOTE: Some of this importS are needed by the mapper,
# not directly in the code
from os.path import exists
from pyworkflow.em.data import (EMObject, EMSet, Micrograph, SetOfMicrographsBase,
                                Acquisition, Particle, Coordinate)
from pyworkflow.object import Float, Pointer, Integer, String, Object 
from pyworkflow.em.convert import ImageHandler


class Tomogram(Micrograph):
    """ Represent a set of frames of micrographs.
    """
    def __init__(self, filename=None, **kwargs):
        Micrograph.__init__(self, filename=filename, **kwargs)
    
    def isCompressed(self):
        return self.getFileName().endswith('bz2') 
        
    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)
        Consider compressed Movie files"""
        if not self.isCompressed():
            fn = self.getFileName()
            if fn is not None and exists(fn.replace(':mrc', '')):
                return ImageHandler().getDimensions(self)[3]
        return None


class SetOfTomograms(SetOfMicrographsBase):
    """ Represents a set of Tomograms. """
    ITEM_TYPE = Tomogram
    
    def __init__(self, **kwargs):
        SetOfMicrographsBase.__init__(self, **kwargs)
        self._totalDose = Float(0)
        self._bfactor = Float(1)
    
    def getDose(self):
        return self._totalDose.get()
    
    def setDose(self, dose):
        self._totalDose.set(dose)
    
    def getBfactorr(self):
        return self._bfactor.get()
    
    def setBfactor(self, bfactor):
        self._bfactor.set(bfactor)
    
    def __str__(self):
        """ String representation of a set of movies. """
        sampling = self.getSamplingRate()

        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0
        ####self._firstFrameNum.set(self.getDimensions()[3])
        if self._firstDim.isEmpty():
            try:
                self._firstDim.set(self.getFirstItem().getDim())
            except Exception, ex:
                print "Error reading dimension: ", ex
                import traceback
                traceback.print_exc()
        dimStr = str(self._firstDim)
        s = "%s (%d items, %d tilt series, %s, %0.2f A/px)" % (self.getClassName(), self.getSize(), self.getFirstItem().getDim(), dimStr, sampling)
        return s