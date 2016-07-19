# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from collections import OrderedDict
import numpy as np

from pyworkflow.em import ImageHandler
from pyworkflow.em.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ
from pyworkflow.em.transformations import (translation_from_matrix,
                                           euler_matrix, euler_from_matrix)
from pyworkflow.em.data import Transform

from simple import SimpleDocFile

    
SHIFTX = 'x'
SHIFTY = 'y'

ANGLE_PSI = 'e1' # in-plane, xmipp psi
ANGLE_THE = 'e2' # tilt in xmipp
ANGLE_PHI = 'e3' # rot in xmipp

CLASS = 'class'

FLIP = 'flip'


def writeSetOfParticles(imgSet, stackFn, docFn, ctfFn):
    """ This function will write a SetOfParticles as a Spider stack and selfile.
    Params:
        imgSet: the SetOfParticles instance.
        stackFn: the filename where to write the stack.
        docFn: the filename to write the information about the particles.
        ctfFn: the filename to write the ctf information
    """
    imgSet.writeStack(stackFn, applyTransform=True)

    writeCTF = imgSet.hasCTF() and ctfFn is not None
    writeDoc = docFn is not None

    if writeCTF:
        ctfDoc = SimpleDocFile(ctfFn, 'w+')
        ctfRow = OrderedDict()

    if writeDoc:
        doc = SimpleDocFile(docFn, 'w+')
        row = OrderedDict()

    if writeCTF or writeDoc:
        for particle in imgSet:

            if writeCTF:
                ctfModelToRow(particle.getCTF(), ctfRow)
                ctfDoc.writeRow(ctfRow)

            if writeDoc:
                pass

    if writeCTF:
        ctfDoc.close()

    if writeDoc:
        doc.close()


def writeSetOfClasses2D(clsSet, clsStack, stackFn, docFn, ctfFn):
    """ This function will write a SetOfClasses2D as a MRC stack and docfile.
    Params:
        clsSet: the SetOfClasses2D instance.
        stackFn: the filename where to write the stack.
        docFn: the filename to write the information about the particles.
        ctfFn: the filename to write the ctf information
    """
    ih = ImageHandler()
    i = 0
    writeCTF = clsSet.getFirstItem().hasCTF() and ctfFn is not None
    writeParticles = stackFn is not None

    if writeParticles:
        doc = SimpleDocFile(docFn, 'w+')
        row = OrderedDict()
    if writeCTF:
        ctfDoc = SimpleDocFile(ctfFn, 'w+')
        ctfRow = OrderedDict()

    for c, cls2D in enumerate(clsSet):
        ih.convert(cls2D.getRepresentative(), (c+1, clsStack))
        if writeCTF or writeParticles:
            for particle in cls2D:
                i += 1
                if writeParticles:
                    row[CLASS] = c + 1
                    ih.convert(particle, (i, stackFn))
                    alignmentToRow(particle.getTransform(), row, ALIGN_2D)
                    doc.writeRow(row)
                if writeCTF:
                    ctfModelToRow(particle.getCTF(), ctfRow)
                    ctfDoc.writeRow(ctfRow)

    if writeParticles:
        doc.close()
    if writeCTF:
        ctfDoc.close()

def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    # Convert defocus values to microns
    ctfRow['dfx'] = ctfModel.getDefocusU() / 10000
    ctfRow['dfy'] = ctfModel.getDefocusV() / 10000
    ctfRow['angast'] = ctfModel.getDefocusAngle()


#-------------- Geometry conversions -----------------

def geometryFromMatrix(matrix, inverseTransform):

    if inverseTransform:
        from numpy.linalg import inv
        matrix = np.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    rawAngles = -np.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    
    # Try to have always positives angles
    # angles = [a + 360 if a < 0 else a for a in rawAngles]
    angles = rawAngles
    
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    radAngles = -np.deg2rad(angles)

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')

    if inverseTransform:
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def rowToAlignment(alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
        """
    is2D = alignType == ALIGN_2D
    inverseTransform = alignType != ALIGN_PROJ

    alignment = Transform()
    angles = np.zeros(3)
    shifts = np.zeros(3)
    #flip = alignmentRow.getValue(xmipp.MDL_FLIP)
    flip = None

    shifts[0] = float(alignmentRow[SHIFTX])
    shifts[1] = float(alignmentRow[SHIFTY])

    if not is2D:
        raise Exception('rowToAlignment: 3D is not implemented yet.')
    else:
        psi = float(alignmentRow[ANGLE_PSI])
        rot = float(alignmentRow[ANGLE_PHI])

        if rot != 0. and psi != 0:
            raise Exception("rot and psi are different from zero in 2D case")

        angles[0] = psi + rot # One of two angles should be zero in 2D

    matrix = matrixFromGeometry(shifts, angles, inverseTransform)

    if flip:
        if alignType == ALIGN_2D:
            matrix[0, :2] *= -1.  # invert only the first two columns keep x
            matrix[2, 2] = -1.  # set 3D rot
        elif alignType == ALIGN_3D:
            matrix[0, :3] *= -1.  # now, invert first line excluding x
            matrix[3, 3] *= -1.
        elif alignType == ALIGN_PROJ:
            pass
            # matrix[0,:4] *= -1.#now, invert first line including x
            ##
    alignment.setMatrix(matrix)

    return alignment


def alignmentToRow(alignment, alignmentRow, alignType):
    """
    alignType == ALIGN_2D-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    """
    is2D = alignType == ALIGN_2D
    inverseTransform = alignType != ALIGN_PROJ
    #only flip is meaninfull if 2D case
    #in that case the 2x2 determinant is negative
    flip = False
    matrix = alignment.getMatrix()
    
    if alignType == ALIGN_2D:
        #get 2x2 matrix and check if negative
        flip = bool(np.linalg.det(matrix[0:2,0:2]) < 0)
        if flip:
            matrix[0,:2] *= -1.#invert only the first two columns keep x
            matrix[2,2]   =  1.#set 3D rot
        else:
            pass

    elif alignType==ALIGN_3D:
        flip = bool(np.linalg.det(matrix[0:3,0:3]) < 0)
        if flip:
            matrix[0,:4] *= -1.#now, invert first line including x
            matrix[3,3]   =  1.#set 3D rot
        else:
            pass

    else:
        flip = bool(np.linalg.det(matrix[0:3,0:3]) < 0)
        if flip:
            matrix[0,:4] *= -1.#now, invert first line including x
    shifts, angles = geometryFromMatrix(matrix, inverseTransform)
    alignmentRow[SHIFTX] = -shifts[0]
    alignmentRow[SHIFTY] = -shifts[1]
    
    if is2D:
        angle = angles[0] + angles[2]
        alignmentRow[ANGLE_PHI] = angle
    else:
        #if alignType == ALIGN_3D:
        alignmentRow[ANGLE_PHI] = angles[0]
        alignmentRow[ANGLE_THE] = angles[1]
        alignmentRow[ANGLE_PSI] = -angles[2]
        
    alignmentRow[FLIP] = -1 if flip else 1
    
    