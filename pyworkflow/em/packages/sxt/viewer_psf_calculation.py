# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
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

import os
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from protocol_psf_calculation import ProtPsfCalculation
from pyworkflow.gui.plotter import plt
import numpy as np
import pyworkflow.protocol.params as params
import pickle
from pyworkflow.em.viewer import ChimeraClientView, ObjectView


VOLUME_SLICES = 0
VOLUME_CHIMERA = 1


class ProtPsfCalculationViewer(ProtocolViewer):
    """ Wrapper to visualize different MTF plots, also 3D PSF
    """
    
    _label = 'viewer calculating PSF from SS'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtPsfCalculation]     
   
        
    def _defineParams(self, form):
        form.addSection(label='Show MTF and 3D PSF')
        form.addParam('mtfNumber', params.IntParam, default=80,
                      label="MTF curve number to plot")
        form.addParam('doShowMtf', params.LabelParam,
                      label="Display the MTF curve")
        form.addParam('displayPsfWith', params.EnumParam, 
                      choices=['slices', 'chimera'],
                      display=params.EnumParam.DISPLAY_HLIST, 
                      default=VOLUME_SLICES,
                      label='Display 3D PSF with',
                      help='*slices*: display 3D PSF as 2D slices along z axis.\n'
                           '*chimera*: display 3D PSF as surface with Chimera.')
        form.addParam('doShowPsf', params.LabelParam,
                      label="Display the 3D PSF")
    
    def _getVisualizeDict(self):
        return {'doShowMtf': self._visualizeMTF,
                'doShowPsf': self._visualizePSF}
    
    def _visualizeMTF(self, e=None):        
        fnMtfDic = self.protocol._defineMtfDicName()
        if not os.path.exists(fnMtfDic):
            return [self.errorMessage('The necessary MTF dic was not produced.\n'
                                      'Execute again the protocol.\n',
                                      title='Missing MTF Dic file')]
                    
        mtfDict = pickle.load(open(fnMtfDic, "rb"))
        fx = mtfDict['fx']
        mtf = mtfDict['mtf']        
        mtfNumber = self.mtfNumber.get()
        mtfToPlot = mtf[:, mtfNumber]        
        mtfCurve = plt.plot(fx, mtfToPlot,'-')
        plt.xlabel('Frequency (1/nm)')
        plt.ylabel('MTF (AU)')
        plt.title('MTF based on a single image')
        plt.grid()
        plt.show()
        return mtfCurve
    
    def _visualizePSF(self, e=None):
        fnPsf = self.protocol._defineOutputName()
        if self.displayPsfWith == VOLUME_CHIMERA:
            view = ChimeraClientView(fnPsf, showProjection=True)
            return [view]
        
        elif self.displayPsfWith == VOLUME_SLICES:
            return [ObjectView(self._project, self.protocol.strId(), fnPsf)]
        
    def _validate(self):
        errors = []        
        mtfNumber = self.mtfNumber.get()
        fnMtfDic = self.protocol._defineMtfDicName()
        mtfDict = pickle.load(open(fnMtfDic, "rb"))
        mtf = mtfDict['mtf']
        if mtfNumber > (np.shape(mtf)[1]-1):
            errors.append("The MTF curve number is exceeded "
                          "from input SS number of images!!!")        
        return errors
    
