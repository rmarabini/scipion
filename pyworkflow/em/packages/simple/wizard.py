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

from pyworkflow.em.wizard import *
from protocol_prime2d import ProtPrime2D
from protocol_prime3d_initial import ProtPrime3DInitial
from protocol_prime3d_refine import ProtPrime3DRefine


class SimpleParticleMaskRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtPrime2D, ['maskRadius']),
                (ProtPrime3DInitial, ['maskRadius']),
                (ProtPrime3DRefine, ['maskRadius'])]
    
    def _getParameters(self, protocol):
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = {}
        if isinstance(protocol, ProtPrime3DInitial):
            inputSet = protocol.inputSet
        else:
            inputSet = protocol.inputParticles
        protParams['input']= inputSet
        protParams['label']= label
        protParams['value']= value
        return protParams
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, UNIT_PIXEL)


class SimpleFilterParticlesWizard(FilterParticlesWizard):
    _targets = [(ProtPrime2D, ['lowPassFilter'])]

    def _getParameters(self, protocol):
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = {'unit': UNIT_ANGSTROM}

        protParams['input'] = protocol.inputParticles
        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = 0 # Low pass filter mode
        return protParams

    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']
        return FilterParticlesWizard._getListProvider(self, _objs)

    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        _mode = params['mode']
        _unit = params['unit']
        FilterParticlesWizard.show(self, form, _value, _label, _mode, _unit)