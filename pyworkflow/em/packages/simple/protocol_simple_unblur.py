# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Hans Elmlund (hans.elmlund@monash.edu) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Dept. of Biochemistry and Molecular Biology, Monash University
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

import os, sys

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.protocol import STEPS_PARALLEL, STEPS_SERIAL

import simple


class ProtSimpleUnblur(ProtAlignMovies):
    """
    """
    _label = 'unblur'
    CONVERT_TO_MRC = 'mrc'

    def __init__(self, **args):
        ProtAlignMovies.__init__(self, **args)
        self.stepsExecutionMode = STEPS_SERIAL

    #--------------------------- DEFINE param functions ------------------------
    def _defineAlignmentParams(self, form):
        form.addParam('binFactor', params.FloatParam, default=1.,
                      label='Binning factor',
                      help='Pass a value greater than 1 to down-scale the '
                           'movies. The parameter for the underlying program '
                           'will be 1/binFactor as expected by simple. '
                           'Downscaling is particularly useful for super'
                           'resolution movies. ')

        line = form.addLine('Subset of frames',
                            expertLevel=cons.LEVEL_ADVANCED,
                            help='')
        line.addParam('frame0', params.IntParam, default=1,
                      label='from')
        line.addParam('frameN', params.IntParam, default=0,
                      label='to')

        form.addParam('frameavg', params.IntParam, default=1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Group number of frames',
                      help=''
                      )

        line = form.addLine('Low-pass limit (A)',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='')
        line.addParam('lpstart', params.IntParam, default=15,
                      label='start')
        line.addParam('lpstop', params.IntParam, default=8,
                      label='stop')

        form.addParam('maxShift', params.IntParam, default=5,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Max. halfwidth shift (px)',
                      help="")

        form.addParam('nsig', params.FloatParam, default=6,
                      pertLevel=cons.LEVEL_ADVANCED,
                      label='Number of sigmas',
                      help='')

        form.addParam('extraParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""  """)

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        """
nthr     = # OpenMP threads{1} #TODO: check performance with simple_distr_exec
startit  = start iterating from here # relevant if doing many
        """
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
        outMicDwFn = self._getExtraPath(self._getOutputMicWtName(movie))
        outMicFn = self._getExtraPath(self._getOutputMicName(movie))
        outMicThumb = self._getOutputMicThumbnail(movie)

        # Create a .txt file with the movie name
        inputTxt = 'input_micrograph.txt'
        with open(os.path.join(movieFolder, inputTxt), 'w') as f:
            f.write('%s' % movie.getBaseName())

        args = "filetab=%s smpd=%0.3f " % (inputTxt, movie.getSamplingRate())
        #args += "nthr=1 " # We take care of parallelization for now
        args += "nthr=%d " % self.numberOfThreads
        args += "trs=%d " % self.maxShift
        args += "frameavg=%d " % self.frameavg
        args += "nsig=%0.3f " % self.nsig

        if self.frame0 > 1 or self.frameN > 0:
            args += "fromf=%d tof=%d " % (self.frame0, self.frameN)

        if abs(self.binFactor.get() - 1.) > 0.0001:
            args += "scale=%0.3f " % (1 / self.binFactor.get())

        args += "%s " % self.extraParams.get('')

        try:
            self.runJob(simple.getProgram('unblur'), args, cwd=movieFolder)

            def _getFile(suffix):
                return os.path.join(movieFolder, '1_%s.mrc' % suffix)

            def _moveFile(suffix, dest):
                pwutils.moveFile(_getFile(suffix), dest)

            _moveFile('forctf', outMicFn)
            _moveFile('intg',  outMicDwFn)
            _moveFile('thumb', outMicThumb)
            self.__runEman2Program('e2proc2d.py', "%s %s"
                                   % (_getFile('pspec'),
                                      self._getPsdCorr(movie)))

            #for suffix in ['forctf', 'intg', 'pspec', 'thumb']
        except:
            print("ERROR: Movie %s failed\n" % movie.getName())

    def _preprocessOutputMicrograph(self, mic, movie):
        mic.thumbnail = em.Image(location=self._getOutputMicThumbnail(movie))
        mic.psdCorr = em.Image(location=self._getPsdCorr(movie))

    def _getOutputMicThumbnail(self, movie):
        baseThumb = ProtAlignMovies._getOutputMicThumbnail(self, movie)
        return pwutils.replaceExt(baseThumb, 'mrc')

    def _getPsdCorr(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_psd.png')

    def _createOutputWeightedMicrographs(self):
        return True

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        return []

    #--------------------------- UTILS functions ------------------------------

