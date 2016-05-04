# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This sub-package contains data and protocol classes
wrapping imod programs
"""
import os
from os.path import join
from pyworkflow.utils import Environ

EXTRACTTILTS = 'extracttilts'
NEWSTACK = 'newstack'

EXTRACTTILTS_PATH = join(os.environ['IMOD_DIR'], 'bin', EXTRACTTILTS)
NEWSTACK_PATH = join(os.environ['IMOD_DIR'], 'bin', NEWSTACK)

def getEnviron():
    """ Setup the environment variables needed to launch imod. """
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['IMOD_DIR'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['IMOD_DIR'], 'lib'),
            }, position=Environ.BEGIN)
    return environ
