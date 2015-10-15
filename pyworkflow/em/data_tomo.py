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
from pyworkflow.em.data import (EMObject, EMSet, SetOfImages,
                                SetOfMicrographsBase, Image)
from pyworkflow.object import Float, Pointer, Integer, String, Object 
from pyworkflow.em.convert import ImageHandler


class TomogramBase(Image):
    """ Create a base class for both Tomogram and TomoRec,
    but avoid to select TomoRec when Tomogram are required. 
    """
    def __init__(self, filename=None, **kwargs):
        Image.__init__(self, filename=filename, **kwargs)
        self._tomoName = String()
    
    def setTomoName(self, micName):
        self._tomoName.set(micName)
    
    def getTomoName(self):
        if self._tomoName.get():
            return self._tomoName.get()
        else:
            self.getFileName()
    
    def copyInfo(self, other):
        """ Copy basic information """
        Image.copyInfo(self, other)
        self.setTomoName(other.getTomoName())
    
    def isCompressed(self):
        return self.getFileName().endswith('bz2') 
    
    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)
        Consider compressed Movie files"""
        if not self.isCompressed():
            fn = self.getFileName()
            if fn is not None:
                return ImageHandler().getDimensions(self)
        return None


class Tomogram(TomogramBase):
    """ Represent a tomogram tilt series.
    """
    pass


class TomoRec(TomogramBase):
    """ Represent a reconstructed tomogram.
    """
    def __init__(self, filename=None, **kwargs):
        TomogramBase.__init__(self, filename=filename, **kwargs)
        self._tomogramPointer = Pointer(objDoStore=False)
    
    def getTomogram(self):
        """ Return the micrograph object to which
        this coordinate is associated.
        """
        return self._tomogramPointer.get()
    
    def setTomogram(self, tomogram):
        """ Set the tomogram to which this reconstruction belongs. """
        self._tomogramPointer.set(tomogram)


class Subtomogram(Image):
    """ Represents an EM Particle object """
    def __init__(self, **args):
        Image.__init__(self, **args)
        self._ctfModel = None
        self._tomoCoordinate = None
        self._tomoId = Integer()
        self._tomoName = String()
        
    def hasCoordinate(self):
        return self._tomoCoordinate is not None
    
    def setCoordinate(self, coordinate):
        self._tomoCoordinate = coordinate
        
    def getCoordinate(self):
        return self._tomoCoordinate
    
    def getTomoId(self):
        """ Return the tomogram id if the coordinate is not None.
        or have set the _tomoId property.
        """
        if self._tomoId.hasValue():
            return self._tomoId.get()
        if self.hasCoordinate():
            return self.getCoordinate().getTomoId()
        
        return None
    
    def setTomoId(self, micId):
        self._tomoId.set(micId)
        
    def hasTomoId(self):
        return self.getTomoId() is not None
    
    def setTomoName(self, tomoName):
        self._micName.set(tomoName)
    
    def getTomoName(self):
        return self._tomoName.get()
    
    def hasCTF(self):
        return self._ctfModel is not None
    
    def getCTF(self):
        """ Return the CTF model """
        return self._ctfModel
    
    def setCTF(self, newCTF):
        self._ctfModel = newCTF


class TomoCoordinate(EMObject):
    """This class holds the (x,y, z) position and other information
    associated with a tomogram coordinate"""
    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._x = Integer(kwargs.get('x', None))
        self._y = Integer(kwargs.get('y', None))
        self._z = Integer(kwargs.get('z', None))
        self._tomoRecPointer = Pointer()
        self._tomoId = Integer()
        self._tomoName = String()
    
    def getX(self):
        return self._x.get()
    
    def setX(self, x):
        self._x.set(x)
        
    def shiftX(self, shiftX):
        self._x.sum(shiftX)
    
    def getY(self):
        return self._y.get()
    
    def setY(self, y):
        self._y.set(y)
        
    def getZ(self):
        return self._z.get()
    
    def setZ(self, z):
        self._z.set(z)
    
    def getPosition(self):
        """ Return the position of the coordinate as a (x, y) tuple.
        mode: select if the position is the center of the box
        or in the top left corner.
        """
        return (self.getX(), self.getY(), self.getZ())

    def setPosition(self, x, y, z):
        self.setX(x)
        self.setY(y)
        self.setZ(z)
    
    def getTomoRec(self):
        """ Return the tomoRec object to which
        this coordinate is associated.
        """
        return self._tomoRecPointer.get()
    
    def setTomoRec(self, tomoRec):
        """ Set the tomogram to which this reconstruction belongs. """
        self._tomoRecPointer.set(tomoRec)
        self._tomoId.set(tomoRec.getObjId())
        self._tomoName.set(tomoRec.getTomoName())
    
    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setObjId(coord.getObjId())
        
    def getTomoId(self):
        return self._tomoId.get()
    
    def setTomoId(self, tomoId):
        self._tomoId.set(tomoId)
    
    def setTomoName(self, tomoName):
        self._micName.set(tomoName)
    
    def getTomoName(self):
        return self._tomoName.get()


class CTF3DModel(EMObject):
    """ Represents a generic CTF model. """
    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._ctfFile = String()
#         self._tomogramPointer  = Pointer()
        self._tomoCoordinatePointer  = Pointer()
    
    def getCtfFile(self):
        self._ctfFile.get()
    
    def setCtfFile(self, ctfFile):
        self._ctfFile.set(ctfFile)
    
#     def getTomogram(self):
#         """ Return the tomoRec object to which
#         this coordinate is associated.
#         """
#         return self._tomogramPointer.get()
#     
#     def setTomogram(self, tomo):
#         """ Set the tomogram to which this reconstruction belongs. """
#         self._tomogramPointer.set(tomo)
    
    def getTomoCoordinate(self):
        """ Return the tomoRec object to which
        this coordinate is associated.
        """
        return self._tomoCoordinatePointer.get()
    
    def setTomoCoordinate(self, tomoCoordinate):
        """ Set the tomogram to which this reconstruction belongs. """
        self._tomoCoordinatePointer.set(tomoCoordinate)


class SetOfTomogramsBase(SetOfMicrographsBase):
    """ Create a base class to set of Tomograms and a set of TomoRecs. """
    
    def __init__(self, **kwargs):
        SetOfMicrographsBase.__init__(self, **kwargs)
        self._totalDose = Float(0)
        self._bfactor = Float(1)
    
    def getDose(self):
        return self._totalDose.get()
    
    def setDose(self, dose):
        self._totalDose.set(dose)
    
    def getBfactor(self):
        return self._bfactor.get()
    
    def setBfactor(self, bfactor):
        self._bfactor.set(bfactor)
    
    def copyInfo(self, other):
        SetOfMicrographsBase.copyInfo(self, other)
        self.setDose(other.getDose())
        self.setBfactor(other.getBfactor())
    
    def __str__(self):
        """ String representation of a set of movies. """
        sampling = self.getSamplingRate()
        
        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0
        if self._firstDim.isEmpty():
            try:
                self._firstDim.set(self.getFirstItem().getDim())
            except Exception, ex:
                print "Error reading dimension: ", ex
                import traceback
                traceback.print_exc()
        dimStr = str(self._firstDim)
        s = "%s (%d items, %d tilt series, %s, %0.2f A/px)" % (self.getClassName(), self.getSize(), self.getFirstItem().getDim()[3], dimStr, sampling)
        return s


class SetOfTomograms(SetOfTomogramsBase):
    """Represet a set of Tomograms"""
    ITEM_TYPE = Tomogram
    
    def __init__(self, **kwargs):
        SetOfTomogramsBase.__init__(self, **kwargs)


class SetOfTomoRecs(SetOfTomogramsBase):
    """Represet a set of TomogRecs"""
    ITEM_TYPE = TomoRec
    
    def __init__(self, **kwargs):
        SetOfTomogramsBase.__init__(self, **kwargs)
        self._tomogramsPointer  = Pointer()
    
    def getTomograms(self):
        """ Return the tomoRec object to which
        this coordinate is associated.
        """
        return self._tomogramsPointer.get()
    
    def setTomograms(self, tomoSet):
        self._tomogramsPointer.set(tomoSet)
    
    def __str__(self):
        """ String representation of a set of movies. """
        sampling = self.getSamplingRate()
         
        if not sampling:
            print "FATAL ERROR: Object %s has no sampling rate!!!" % self.getName()
            sampling = -999.0
        if self._firstDim.isEmpty():
            try:
                self._firstDim.set(self.getFirstItem().getDim())
            except Exception, ex:
                print "Error reading dimension: ", ex
                import traceback
                traceback.print_exc()
        dimStr = str(self.getFirstItem().getDim()[0]) + " x " + str(self.getFirstItem().getDim()[1]) + " x " + str(self.getFirstItem().getDim()[3])
        s = "%s (%d items, %s, %0.2f A/px)" % (self.getClassName(), self.getSize(), dimStr, sampling)
        return s


class SetOfSubtomograms(SetOfImages):
    """ Represents a set of Subtomograms.
    The purpose of this class is to separate the
    concepts of reconstructed tomograms and subtomograms, even if
    both are considered Images
    """
    ITEM_TYPE = Subtomogram
    REP_TYPE = Subtomogram
    
    def __init__(self, **args):
        SetOfImages.__init__(self, **args)
        self._tomoCoordsPointer = Pointer()
        
        
    def hasCoordinates(self):
        return self._tomoCoordsPointer.hasValue()
    
    def getCoordinates(self):
        """ Returns the SetOfTomoCoordinates associated with 
        this SetOfSubtomograms"""
        return self._tomoCoordsPointer.get()
    
    def setCoordinates(self, coordinates):
        """ Set the SetOfCoordinates associates with
        this set of particles.
         """
        self._tomoCoordsPointer.set(coordinates)    
    
    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and sampling rate)
        from other set of subtomograms to current one.
        """
        SetOfImages.copyInfo(self, other)
        self.setHasCTF(other.hasCTF())


class SetOfTomoCoordinates(EMSet):
    """ Encapsulate the logic of a subtomograms coordinates.
    Each coordinate has a (x,y,z) position and is related to a reconstructed tomogram.
    """
    ITEM_TYPE = TomoCoordinate
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._setOfTomoRecsPointer = Pointer()
        self._boxSize = Integer()

    def getBoxSize(self):
        """ Return the box size of the subtomograms.
        """
        return self._boxSize.get()
    
    def setBoxSize(self, boxSize):
        """ Set the box size of the subtomograms. """
        self._boxSize.set(boxSize)
    
    def iterTomoRecs(self):
        """ Iterate over the micrographs set associated with this
        set of coordinates.
        """
        return self.getTomoRecs()
    
    def iterCoordinates(self, tomoRec=None):
        """ Iterate over the coordinates associated with a tomoRec.
        If tomoRec=None, the iteration is performed over the whole set of coordinates.
        """
        if tomoRec is None:
            tomoId = None
        elif isinstance(tomoRec, int):
            tomoId = tomoRec
        elif isinstance(tomoRec, TomoRec):
            tomoId = tomoRec.getObjId()
        else:
            raise Exception('Invalid input reconstructed tomogram of type %s' % type(tomoRec))
        
        #Iterate over all coordinates if micId is None,
        #otherwise use micId to filter the where selection
        coordWhere = '1' if  tomoId is None else '_tomoId=%d' % tomoId
        for coord in self.iterItems(where=coordWhere):
            yield coord
    
    def getTomoRecs(self):
        """ Returns the SetOfTomoRecs associated with 
        this SetOfTomoCoordinates"""
        return self._setOfTomoRecsPointer.get()
    
    def setTomoRecs(self, tomoRecs):
        """ Set the SetOfTomoRecs associates with 
        this set of coordinates.
         """
        self._setOfTomoRecsPointer.set(tomoRecs)
    
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths
    
    def __str__(self):
        """ String representation of a set of coordinates. """
        if self._boxSize.hasValue():
            boxSize = self._boxSize.get()
            boxStr = ' %d x %d x %d' % (boxSize, boxSize, boxSize)
        else:
            boxStr = 'No-Box'
        s = "%s (%d items, %s)" % (self.getClassName(), self.getSize(), boxStr)
        
        return s


class SetOfCTF3D(EMSet):
    """ Contains a set of CTF models estimated."""
    ITEM_TYPE = CTF3DModel
    
    def __init__(self, **args):
        EMSet.__init__(self, **args)
        self._tomoCoordinatesPointer = Pointer()
    
    def getTomoCoordinates(self):
        return self._tomoCoordinatesPointer.get()
    
    def setTomoCoordinates(self, tomoCoordinates):
        self._tomoCoordinatesPointer.set(tomoCoordinates)
