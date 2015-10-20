# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
# from itertools import izip
"""
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from os.path import join
import pyworkflow.utils as putils
from pyworkflow.em import metadata as md
from pyworkflow.em import constants as const


def getEnviron():
    """ Setup the environment variables needed to launch Relion. """
    environ = putils.Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['RELION_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib') + ":" + join(os.environ['RELION_HOME'], 'lib64'),
            'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
            }, position=putils.Environ.BEGIN)
    return environ


def writeSetOfTomograms(tomoRecSet, starFile, dirPath):
    tomoMd = md.MetaData()
    for tomoRec in tomoRecSet:
        dest = os.path.join(dirPath, os.path.basename(tomoRec.getFileName()))
        putils.createLink(tomoRec.getFileName(), dest)
        objId = tomoMd.addObject()
        imgRow = md.Row()
        imgRow.setValue(md.RLN_MICROGRAPH_NAME, join("extra", os.path.basename(tomoRec.getFileName())))
        imgRow.writeToMd(tomoMd, objId)
    tomoMd.write(starFile)


def writeSetOfSubtomograms(subtomoSet, subtomoStar, path):
    subtomoMd = md.MetaData()
    setOfImagesToMd(subtomoSet, subtomoMd)
    subtomoMd.write(subtomoStar)


def setOfImagesToMd(imgSet, imgMd):
    """ This function will fill Relion metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        md: metadata to be filled
        rowFunc: this function can be used to setup the row before 
            adding to meta
    """
    for img in imgSet:
        objId = imgMd.addObject()
        imgRow = md.Row()
        subtomoToRow(img, imgRow)
        imgRow.writeToMd(imgMd, objId)


def subtomoToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    coord = part.getCoordinate()
    partRow.setValue(md.RLN_MICROGRAPH_NAME, coord.getTomoRec().getFileName())
    partRow.setValue(md.RLN_IMAGE_COORD_X, float(coord.getX()))
    partRow.setValue(md.RLN_IMAGE_COORD_Y, float(coord.getY()))
    partRow.setValue(md.RLN_IMAGE_COORD_Z, float(coord.getZ()))
    partRow.setValue(md.RLN_IMAGE_NAME, part.getFileName())
    if part.hasCTF():
        partRow.setValue(md.RLN_CTF_IMAGE, part.getCTF().getCtfFile())


def relionToLocation(filename):
    """ Return a location (index, filename) given
    a Relion filename with the index@filename structure. """
    if '@' in filename:
        indexStr, fn = filename.split('@')
        return int(indexStr), str(fn)
    else:
        return const.NO_INDEX, str(filename)


def convertBinaryVol(vol, outputDir):
    from pyworkflow.em import ImageHandler
    """ Convert binary volume to a format read by Relion.
    Params:
        vol: input volume object to be converted.
        outputDir: where to put the converted file(s)
    Return:
        new file name of the volume (convrted or not).
    """
    
    ih = ImageHandler()
    # This approach can be extended when
    # converting from a binary file format that
    # is not read from Relion
    def convertToMrc(fn):
        """ Convert from a format that is not read by Relion
        to mrc format.
        """
        newFn = join(outputDir, putils.replaceBaseExt(fn, 'mrc'))
        ih.convert(fn, newFn)
        return newFn
        
    ext = vol.getFileName()
    
    if not ext.endswith('.mrc'):
        fn = convertToMrc(vol.getFileName())
    else:
        fn = vol.getFileName()
    return fn
