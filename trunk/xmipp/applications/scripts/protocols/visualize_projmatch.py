#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_projmatch.py
#
# Example use:
# python visualize_projmatch.py
#
# This script requires that protocol_projmatch.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
#show results for iteration
DisplayIterNo=1
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=False
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=True
# {expert} Width of Matrix-views (multiple of three value!):
MatrixWidth=10
#------------------------------------------------------------------------------------------------
# {section}Reference volume 
#------------------------------------------------------------------------------------------------
#show reference volume 
DisplayReference=False
# Plot the angular distribution of the reference(s)?
VisualizeAngDistribution=False
#Show projection maching library and alignes classes
DisplayProjectionMatching=False
#display angular distribution
DisplayAngularDistribution=False
#display reconstructed volume
DisplayReconstruction=False



#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_projmatch_class:

    #init variables
    def __init__(self,
                _VisualizeVolZ,
                _VisualizeVolX,
                _VisualizeVolY,
                _VisualizeVolChimera,
                _DisplayIterNo,
                _DisplayReference,
                _DisplayProjectionMatching,
                _DisplayReconstruction,
                _MatrixWidth,
                _VisualizeAngDistribution
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log
        import visualization

        # import corresponding protocol
        import protocol_projmatch
        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol_projmatch.WorkDirectory
        self._LogDir=protocol_projmatch.LogDir
        self._ProjectDir=protocol_projmatch.ProjectDir
        self._ReferenceVolume=protocol_projmatch.ReferenceVolume
        self._multi_align2d_sel=protocol_projmatch.multi_align2d_sel
        self._SelFileName=self._ProjectDir+'/'+str(protocol_projmatch.SelFileName)
        self._ReferenceVolume=protocol_projmatch.ReferenceVolume
        self._Proj_Maching_Output_Root_Name=protocol_projmatch.Proj_Maching_Output_Root_Name
        self._Proj_Maching_Output_Root_Name + '.doc'
        self.mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._iteration_number=_DisplayIterNo
        self._Iteration_Working_Directory=self._WorkDirectory+'/Iter_'+\
                                   str(self._iteration_number)
        #os.chdir(Iteration_Working_Directory)

        if (_DisplayReference):
           self.ShowVolumes=[] 
           self.ShowVolumes.append(os.getcwd()+'/'+\
                                   self._Iteration_Working_Directory+'/'+\
                                   self._ReferenceVolume)
           visualization.visualize_volumes(self.ShowVolumes,
                                                _VisualizeVolZ,
                                                _VisualizeVolX,
                                                _VisualizeVolY,
                                                _VisualizeVolChimera)
        if (_VisualizeAngDistribution):
           self._ShowPlots=[] 
           self._ShowPlots.append(os.getcwd()+'/'+\
                                  self._Iteration_Working_Directory+'/'+\
                                  self._Proj_Maching_Output_Root_Name+\
                                  ".doc")
           show_ang_distribution(self._ShowPlots,self._iteration_number)
           
        if (_DisplayProjectionMatching):
            self.ShowSelfiles=[] 
            self.ShowSelfiles.append(self._Iteration_Working_Directory+'/'+\
                                     self._multi_align2d_sel)
            visualization.visualize_images(self.ShowSelfiles,
                             True,
                             self._MatrixWidth)
        #if (_DisplayReconstruction):
        #    self.visualize_Reconstruction(self._SomName,self._SpectraName)
            
        # Return to parent dir
        # os.chdir(os.pardir)

def show_ang_distribution(_ShowPlots,_iteration_number):
        import os
        import docfiles
        import visualization
        for plots in _ShowPlots: 
            doc=docfiles.docfile(plots)
            doc.check_angle_range()
            doc.write_several(plots,
                              10,
                              7,
                              doc.minimum_of_column(7),
                              doc.maximum_of_column(7)
                              )
            plot=visualization.gnuplot()
            title='Angular distribution for projection matching for iteration '+\
                    str(_iteration_number)
            plot.plot_xy1y2_several_angular_doc_files(plots,
                                                      title,
                                                      'degrees',
                                                      'degrees')

def close():
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create projmatch_class object
    visualize_projmatch=visualize_projmatch_class(VisualizeVolZ,
                                                  VisualizeVolX,
                                                  VisualizeVolY,
                                                  VisualizeVolChimera,
                                                  DisplayIterNo,
                                                  DisplayReference,
                                                  DisplayProjectionMatching,
                                                  DisplayReconstruction,
                                                  MatrixWidth,
                                                  VisualizeAngDistribution)
    # close 
    visualize_projmatch.close()

