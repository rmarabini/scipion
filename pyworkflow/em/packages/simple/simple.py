# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
from collections import OrderedDict

from pyworkflow.utils import Environ



def getEnviron():
    """ Create the needed environment for SIMPLE programs. """
    environ = Environ(os.environ)
    SIMPLEBIN = os.path.join(os.environ['SIMPLE_HOME'], 'bin')
    SIMPLEAPPS = os.path.join(os.environ['SIMPLE_HOME'], 'apps')
    environ.update({
                    'SIMPLEBIN': SIMPLEBIN,
                    'SIMPLEPATH': os.environ['SIMPLE_HOME'],
                    'SIMPLESYS': os.environ['SIMPLE_HOME'],
                    'PATH': os.pathsep.join([SIMPLEBIN, SIMPLEAPPS])
                    }, 
                   position=Environ.BEGIN)
    return environ


class SimpleDocFile(object):
    """ Handler class to read/write SIMPLE docfile. """

    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def writeRow(self, row):
        """ Write a row to the data file.
        Each row should contains values in the form of key=value.
        """
        self._count += 1

        for value in row.iteritems(): # value is a tuple key, value
            self._file.write("%s=%s " % value)

        self._file.write("\n")

    def iterValues(self):
        for line in self._file:
            line = line.strip()
            if not line.startswith(';'):
                row = OrderedDict()
                for token in line.split():
                    key, value = token.split("=")
                    row[key] = value
                yield row

    def __iter__(self):
        return self.iterValues()

    def getLastRow(self):
        lastRow = None
        for row in self:
            lastRow = row
        return lastRow

    def close(self):
        self._file.close()