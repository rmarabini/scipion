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
from itertools import izip
"""
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from os.path import join
import pyworkflow.utils as putils


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
    import pyworkflow.em.metadata as md
    
    tomoMd = md.MetaData()
    for tomoRec in tomoRecSet:
        dest = os.path.join(dirPath, os.path.basename(tomoRec.getFileName()))
        putils.createLink(tomoRec.getFileName(), dest)
        objId = tomoMd.addObject()
        imgRow = md.Row()
        imgRow.setValue(md.RLN_MICROGRAPH_NAME, join("extra", os.path.basename(tomoRec.getFileName())))
        imgRow.writeToMd(tomoMd, objId)
    tomoMd.write(starFile)


