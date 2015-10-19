# **************************************************************************
# *
# * Author:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This sub-package will contains Relion protocols
"""

from bibtex import _bibtex # Load bibtex dict with references

_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Chen2012']

from convert import getEnviron
from pyworkflow.em.tomo.reliontomo.protocol_classify3d_tomo import ProtRelionSubtomoClassify3D
from pyworkflow.em.tomo.reliontomo.protocol_refine3d_tomo import ProtRelionSubtomoRefine3D
from pyworkflow.em.tomo.reliontomo.protocol_extract_tomo import ProtRelionExtractSubtomograms

# # Wizards
# from wizard import *
# 
# from viewer import *

_environ = getEnviron()