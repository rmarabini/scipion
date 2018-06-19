# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

# 1 copy 8 sections
# param line; boolean and conditional

import os
import sys
import time
import select
import paramiko
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import (BooleanParam,
                                        FloatParam,
                                        IntParam,
                                        MultiPointerParam,
                                        PathParam,
                                        StringParam
                                        )



from pyworkflow.utils.remote import


class ProtBackup(EMProtocol):
    """ copy project, or part of a project, to another disk or computer
    """
    _label = 'Backup project'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doAll', BooleanParam, default=True, label="Backup Whole Project?",
                       help='Check to backup the whole project')
        form.addParam('inputProtocolsExclude', MultiPointerParam,
                      label="Exclude these protocols",
                      pointerClass='EMProtocol',
                      default=None, condition='not doAll',
                      help="Protocols NOT to be excluded")
        form.addParam('doLocal', BooleanParam, default=True, label="Local "
                                                                   "Backup",
                       help='Perform local backup (same computer)')
        form.addParam('targetLocal', PathParam, default=None, condition='doLocal',
                      label='Target Directory', help="Path to local directory where "
                                                     "the backup will be stored")
        form.addParam('targetRemote', StringParam, default=None,
                      condition='not doLocal',
                      label='Target Directory', help="Path to remote directory where "
                                                     "the backup will be "
                                                     "stored. Example: "
                                                     "user@host:path")
        #TODO
        form.addSection('Streaming')

        form.addParam('dataStreaming', BooleanParam, default=False,
              label="Process data in streaming?",
              help="Select this option if you want import data as it is "
                   "generated and process on the fly by next protocols. "
                   "In this case the protocol will keep running to check "
                   "new files and will update the output Set, which can "
                   "be used right away by next steps.")

        form.addParam('timeout', IntParam, default=43200,
              condition='dataStreaming',
              label="Timeout (secs)",
              help="Interval of time (in seconds) after which, if no new file "
                   "is detected, the protocol will end. When finished, "
                   "the output Set will be closed and no more data will be "
                   "added to it. \n"
                    "Note 1:  The default value is  high (12 hours) to avoid "
                   "the protocol finishes during the aqcuisition of the "
                   "microscpe. You can also stop it from right click and press "
                   "STOP_STREAMING.\n"
                   "Note 2: If you're using individual frames when importing "
                   "movies, the timeout won't be refreshed until a whole "
                   "movie is stacked.")

        form.addParam('fileTimeout', IntParam, default=30,
              condition='dataStreaming',
              label="File timeout (secs)",
              help="Interval of time (in seconds) after which, if a file has "
                   "not changed, we consider it as a new file. \n")


    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._copyData()
    #--------------- STEPS functions -----------------------

    def _copyData(self, params):
        if self.doLocal:
            print "Local"
            target =self.targetLocal.get()
        else:
            print "Remote"
            target =self.targetRemote.get()

        projName = self.protocol.getProject().getShortName()
        print "projName=", projName

        #create directory
        exitcode = CopyFiles(projName, target, self.timeout)

    #--------------- INFO functions -------------------------

    def _validate(self):
        return []

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    #--------------- UTILS functions -------------------------


# ===========================

"""
        #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('exportVolume', params.PointerParam, label="Volume to export", important=True,
                      pointerClass='Volume',
                      help='This volume will be exported using mrc format')
        form.addParam('exportFSC', params.PointerParam, label="FSC to export", important=True,
                      pointerClass='FSC',
                      help='This FSC will be exported using XML emdb format')
        form.addParam('filesPath', params.PathParam, important=True,
                      label="Export to directory",
                      help="Directory where the files will be generate.")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._createFileNamesTemplates()
        self._insertFunctionStep('exportVolumeStep')
        self._insertFunctionStep('exportFSCStep')

    #--------------------------- STEPS functions --------------------------------------------

    def exportVolumeStep(self):

        ih = ImageHandler()
        ih.convert(self.exportVolume.get().getLocation(), self.getFnPath())

    def exportFSCStep(self):

        x,y = self.exportFSC.get().getData()
        fo = open(self.getFnPath("fsc"), "w")
        fo.write('<fsc title="FSC(%s)" xaxis="Resolution(A-1)" '
                 'yaxis="Correlation Coefficient">\n' % self._getFileName('volume'))
        for i in range(len(x)):
            fo.write("<coordinate>\n")
            fo.write("<x>%f</x>\n"%x[i])
            fo.write("<y>%f</y>\n" % y[i])
            fo.write("</coordinate>\n")

        fo.write("</fsc>\n")
        fo.close()

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        message = []
        fnPath = self.filesPath.get()
        if fnPath == "" or fnPath is None:
            message.append("You must set a path to export.")
        return message

    def _summary(self):
        message = "Data Available at : *%s*"% self.filesPath.get()
        return [message]

    def _methods(self):
        return []

#--------------------------- UTILS functions ---------------------------------------------------

    def getFnPath(self, label='volume'):
        return os.path.join(self.filesPath.get(), self._getFileName(label))

"""