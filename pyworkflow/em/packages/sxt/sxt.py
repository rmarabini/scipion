# **************************************************************************
# *
# * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
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
from os.path import join
from pyworkflow.utils import Environ, runJob


def getEnviron(xmippFirst=True):
    """ Create the needed environment for Xmipp programs. """
    environ = Environ(os.environ)
    pos = Environ.BEGIN if xmippFirst else Environ.END
    environ.update({
            'PATH': join(os.environ['XMIPP_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['XMIPP_HOME'], 'lib')
            }, position=pos)
    return environ

def runXmippProgram(program, args=""):
    """ Internal shortcut function to launch a Xmipp program. """
    runJob(None, program, args, env=getEnviron())

