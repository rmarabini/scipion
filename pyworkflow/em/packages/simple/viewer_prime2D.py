# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
from os.path import exists

import pyworkflow.em as em
import pyworkflow.em.showj as showj
from pyworkflow.em.plotter import EmPlotter
import pyworkflow.protocol.params as params
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO

from protocol_prime_2d import ProtPrime2D
from protocol_prime3D_initial import ProtPrime3DInitial

ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CHIMERADATAVIEW = 0

CLASSES_ALL = 0
CLASSES_SEL = 1

FSC_CORRECTED = 0
FSC_UNMASKEDMAPS = 1
FSC_MASKEDMAPS = 2
FSC_RANDOMIZED = 3
FSC_ALL = 4


class SimpleViewerPrime2D(ProtocolViewer):
    """ This protocol serve to analyze the results of Relion runs.
    (for protocols classify 2d/3d and 3d auto-refine)
    The visualization tools follow the recommendations of Relion 1.3 tutorial:
    http://www2.mrc-lmb.cam.ac.uk/groups/scheres/relion13_tutorial.pdf
    """
    _targets = [ProtPrime2D, ProtPrime3DInitial]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer simple'
    
    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize", 
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', params.NumericRangeParam, 
                      condition='viewIter==%d' % ITER_SELECTION, 
                      label="Iterations list", 
                      help="Write the iteration list to visualize.")

        form.addParam('showImagesInClasses', params.LabelParam,
                      label='Show classification in Scipion', important=True,
                      help='Display each class with the number of particles assigned. \n'
                           '*Note1*: The images of one class can be shown by \n'
                           'right-click on the class and select "Open images".')

        form.addParam('showChanges', params.LabelParam, default=True,
                      label='Plot ',
                      help='')
                                              
    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {
                'showImagesInClasses': self._showImagesInClasses,
                'showChanges': self._showChanges,
                }

        # If the is some error during the load, just show that instead
        # of any viewer
        if self._errors:
            for k in visualizeDict.keys():
                visualizeDict[k] = self._showErrors

        return visualizeDict

    def _showErrors(self, param=None):
        views = []
        self.errorList(self._errors, views)
        return views
        
    def _viewAll(self, *args):
        pass
    
#===============================================================================
# showImagesInClasses     
#===============================================================================
    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to 
        # properly see the rlnClassDistribution label. 
        dim = self.protocol.inputParticles.get().getDim()[0]
        if dim < 128:
            zoom = 128*100/dim
        else:
            zoom = 100
        return zoom
        
    def _showImagesInClasses(self, paramName=None):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        views = []
        
        for it in self._iterations:
            fn = self.protocol.getIterClasses(it)
            v = self.createScipionView(fn)
            views.append(v)
        
        return views
    
    def _showChanges(self, paramName=None):
        logFile, _, _ = self.protocol.getLogPaths()
        f = open(logFile)
        distribution = []
        percentage = []
        correlation = []

        for line in f:
            l = line.strip()
            if l.startswith('>>> DISTRIBUTION OVERLAP:'):
                distribution.append(float(line.split()[-1]))
            elif l.startswith('>>> PERCENTAGE OF SEARCH SPACE SCANNED:'):
                percentage.append(float(line.split()[-1]) / 100)
            elif l.startswith('>>> CORRELATION:'):
                correlation.append(float(line.split()[-1]))

        f.close()

        plotter = EmPlotter()
        plotter.createSubPlot("Values per Iterations", "Iterations", "Value")
        xx = xrange(1, len(distribution)+1)
        plotter.plot(xx, distribution, label='Distribution overlap')
        plotter.plot(xx, percentage, label='Search space scanned')
        plotter.plot(xx, correlation, label='Correlation')
        plotter.legend()

        return [plotter]
    
#===============================================================================
# Utils Functions
#===============================================================================
    def _validate(self):
        if self.lastIter is None:
            return ['There are not iterations completed.'] 
    
    def createScipionView(self, filename):
        labels =  'enabled id _size _representative._filename '
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, 
                      showj.RENDER:'_representative._filename',
                      showj.SORT_BY: '_size desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        inputParticlesId = self.protocol.inputParticles.get().strId()
        view = em.ClassesView(self._project,
                          self.protocol.strId(), filename, other=inputParticlesId,
                          env=self._env,
                          viewParams=viewParams)

        return view

    def _getRange(self, var, label):
        """ Check if the range is not empty.
        :param var: The variable to retrieve the value
        :param label: the labe used for the message string
        :return: the list with the range of values, empty
        """
        value = var.get()
        if value is None or not value.strip():
            self._errors.append('Provide %s selection.' % label)
            result = []
        else:
            result = self._getListFromRangeString(value)

        return result

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1]
        self._errors = []

        self.firstIter = 1
        self.lastIter = self.protocol.getLastIteration()
        
        if self.viewIter == ITER_LAST:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')

