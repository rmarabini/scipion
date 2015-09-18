#!/usr/bin/env python
print ':: RELION sub-tomogram averaging ::' 
print 'This python script was written by Tanmay Bharat to support sub-tomogram averaging in RELION.' 
print 'Please ensure that you have IMOD and RELION commands in your path and that you have CTFFIND installed.' 
print 'Please report bugs and comments to tbharat@mrc-lmb.cam.ac.uk or scheres@mrc-lmb.cam.ac.uk'
print 'Please read the documentation on the RELION wiki, several questions are answered there.' + '\n' 

import os, sys, commands, math, time, stat

######### INPUT #########################################

## Input STAR file with all tomograms
TomogramStarFileName = 'all_tomograms.star'

## Parameters for ctffind
# Microscope voltage in kV
Voltage = 300              
# Spherical aberration coefficient in mm
Cs = 2.7                   
# Magnification of the image
Magnification = 53000       
# Pixel size of the detector (in micron)
DPixSize = 11.5             
# Path to CTFFIND (version 3 or 4)
PathToCtffind = '/public/EM/ctffind/ctffind.exe'     
# If CTFFIND crashed in the middle, you can turn this to True to resume CTF estimations only for unfinished images
OnlyDoUnfinishedCTFs = False                               

## Rootname for the sub-tomograms. Please give the same rootname when you extract sub-tomograms using the RELION GUI later. 
SubtomoName = 'subtomo'


## CTFFIND inputs
BoxSize = 256
LowResLimit = 100
HighResLimit = 10
LowDefocusLimit = 10000
HighDefocusLimit = 90000
DefocusStep = 5000
AmpContrast = 0.1
Astigmatism = 5000

## 3D CTF model weighting B-factor per e-/A2
Bfactor = 4.0


######## FUNCTIONS #######################################

#
def ensure_dir(f):
  d = os.path.dirname(f)
  if not os.path.exists(d):
    #print 'Making directory'
    os.makedirs(d)
#

# To read the STAR files. Please note that this STAR file parser is only meant for setting up the sub-tomogram averaging scripts.
# RELION has a more comprehensive parser in the main code.
def read_relion_star(filename):
  starfile=open(filename, 'r')
  j=-1

  micnames=[]
  defociu=[]
  defociv=[]
  for line in starfile:

    #print line

    emptycheck = line.isspace()
    if(emptycheck):
      #print 'empty line found'
      continue
    
    fields = line.split()
    firstfield = fields[0]
    if firstfield[0] == 'd':
      #print 'data_ line found'
      continue
    if firstfield[0] == 'l':
      #print 'loop_ line found'
      continue
    j=j+1 

    if firstfield == '_rlnMicrographName':
      imgnamecolumn = j
      continue
    if firstfield == '_rlnDefocusU':
      defocusucolumn = j
      continue
    if firstfield == '_rlnDefocusV':
      defocusvcolumn = j
      continue
    #if firstfield == '_rlnCtfFigureOfMerit':
    #  ctffigureofmeritcolumn = j
    #  continue
    if firstfield[0] == '_':
      continue
    

    micnames.append(fields[imgnamecolumn])
    if 'defocusucolumn' in locals():
      defociu.append(fields[defocusucolumn])
      defociv.append(fields[defocusvcolumn])

  starfile.close()
  if len(defociu) > 0:
    return micnames,defociu,defociv
  if len(defociu) == 0:
    return micnames
#

######## FUNCTIONS #######################################




######## RUNNING THE SCRIPT #################


## Looping through the micrographs
ScriptDir = os.getcwd() + '/'
#print ScriptDir

micnames = read_relion_star(TomogramStarFileName)
#print micnames

# Shell script to do 3D CTF model reconstruction
ctfreconstmastername = ScriptDir + 'do_all_reconstruct_ctfs.sh'
ctfreconstmasterfile = open(ctfreconstmastername, 'w')
os.chmod(ctfreconstmastername, stat.S_IRWXU)

#
# This is the master STAR file for refinement later on
subtomostarname = ScriptDir + 'particles_' + SubtomoName + '.star'
subtomostarfile = open(subtomostarname, 'w')
# writing out the header of the list star file
subtomostarfile.write('data_' + '\n' + '\n')
subtomostarfile.write('loop_' + '\n')
subtomostarfile.write('_rlnMicrographName #1' + '\n')
subtomostarfile.write('_rlnCoordinateX #2' + '\n')
subtomostarfile.write('_rlnCoordinateY #3'+ '\n')
subtomostarfile.write('_rlnCoordinateZ #4' + '\n')
subtomostarfile.write('_rlnImageName #5' + '\n')
subtomostarfile.write('_rlnCtfImage #6' +'\n')
#

for mic in micnames:

  #
  # Parsing the micrograph names 
  micsplit = mic.split('.')
  microot = micsplit[0]
  dirsplit = microot.split('/')
  MicDirName = ""
  for dircount in range(0,(len(dirsplit)-1)):
    MicDirName = MicDirName + dirsplit[dircount]
    MicDirName = MicDirName + '/'
  MicRootName = dirsplit[len(dirsplit)-1]
  
  #print MicDirName
  #print MicRootName
  
  micname = MicDirName + MicRootName + '.mrc'
  stackname = MicDirName + MicRootName + '.st'
  ordername = MicDirName + MicRootName + '.order'
  coordsname = MicDirName + MicRootName + '.coords'
  #print micname, stackname, ordername, coordsname
  # Parsing the micrograph names 
  #

  #sys.exit()

  ##### Running CTFFIND on all images of the tilt series  ##########

  CtffindDirName = 'ctffind/'
  OutputDir = ScriptDir + MicDirName + CtffindDirName
  newstackroot = MicDirName + CtffindDirName + MicRootName +  '_image'
  print OutputDir

  ## Making a new directory to output the results of CTFFIND
  ensure_dir(OutputDir)

  ## Extracting the tilt information with the IMOD command extracttilts
  extracttile_scratchname = OutputDir + 'extracttilt_output.txt'
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Using IMOD extracttilts to get tilt angles' + '\n' 
  exttltline = 'extracttilts -InputFile ' + stackname + ' -tilts -OutputFile ' + OutputDir +  'tiltangles.txt > ' + extracttile_scratchname +  '\n'
  print(exttltline)
  os.system(exttltline)
  os.remove(extracttile_scratchname)
  ##
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Tilt values extracted ' + '\n'

  ##
  tiltanglesfilename = OutputDir + 'tiltangles.txt'
  tiltfile = open(tiltanglesfilename, 'r')
  
  ctffindstarname = OutputDir + MicRootName + '_images.star'
  ctffindstarfile = open(ctffindstarname, 'w')
  ctffindstarfile.write('data_' + '\n' + '\n')
  ctffindstarfile.write('loop_' + '\n')
  ctffindstarfile.write('_rlnMicrographName #1' + '\n')
  
  exttilts=[]
  i=-1
  for line in tiltfile:
    
    pair = line.split()
    #print pair
    i=i+1
    
    # Tilt of the stage for the current image
    tilt = float(pair[0])
    #roundtilt = round(tilt)
    exttilts.append(tilt)
    #print(str(int(roundtilt)))
    #

    # extracting each image using the IMOD command newstack
    print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Extracting tilt series image ' + '\n' 
    newstack_scratchname = OutputDir + 'temp_newstack_out.txt'
    extracted_image_name = newstackroot + str(tilt) + '_' + str(i) + '.mrc'
    newstackline = 'newstack -secs ' + str(i) + ' ' + stackname + ' ' +  extracted_image_name + ' > ' + newstack_scratchname +'\n'
    print(newstackline)
    ctffindstarfile.write(extracted_image_name + '\n')

    os.system(newstackline)
    os.remove(newstack_scratchname)

  ctffindstarfile.close()
  #sys.exit()

  # running CTFFIND using the RELION command relion_run_ctffind
  outputstarname =  OutputDir + MicRootName +  '_ctffind.star'
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Running relion_run_ctffind ' + '\n' 
  relion_ctffindline = 'relion_run_ctffind --i ' + ctffindstarname + ' --o ' + outputstarname + ' --CS ' + str(Cs) + ' --HT ' + str(Voltage) +  ' --ctfWin -1 --AmpCnst ' + str(AmpContrast) +  ' --DStep ' + str(DPixSize) +  ' --XMAG ' + str(Magnification) + ' --Box ' + str(BoxSize) +  ' --dFMin ' + str(LowDefocusLimit) + ' --dFMax ' + str(HighDefocusLimit) + ' --FStep ' + str(DefocusStep) + ' --dAst ' + str(Astigmatism) + ' --ResMin ' + str(LowResLimit) + ' --ResMax ' + str(HighResLimit) + ' --ctffind_exe \"' + PathToCtffind + '  --omp-num-threads 1 --old-school-input\"' 
  
  # If some are unfinished 
  if OnlyDoUnfinishedCTFs == True:
    relion_ctffindline = relion_ctffindline + ' --only_do_unfinished'
  print(relion_ctffindline)
  os.system(relion_ctffindline)
          #
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'CTF Parameters of all tilt series images were estimated using RELION\'s  relion_run_ctffind ' + '\n' 
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Parameters have been saved in ' + outputstarname + '\n' 

  tiltfile.close()

  ##### Running CTFFIND on all images of the tilt series  ##########

  #sys.exit()

  ##### Making .star files for each 3D CTF Volume #################
  
  RelionPartName = 'Particles/'
  RelionPartDir = ScriptDir + RelionPartName
  RelionRecDir = RelionPartDir + MicDirName
  RelionRecFileName =  RelionPartName + MicDirName + MicRootName + '_rec_CTF_volumes.sh'
  RelionRecFileName_for_script =  MicRootName + '_rec_CTF_volumes.sh'

  ## Making a new directory to output the results of CTFFIND
  ensure_dir(RelionRecDir)

  coordfile = open(coordsname, 'r')
  relionfile = open(RelionRecFileName, 'w')

  # Getting the tilt order
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Reading tilt series order file for dose dependent B-Factor weighting ' + '\n' 
  tiltorderfile = open(ordername, 'r')
  tiltorder=[]
  accumulated_dose=[]

  for line in tiltorderfile:
    
    emptycheck = line.isspace()
    if(emptycheck):
      #print 'empty line found'
      continue
    
    pair=line.split()
    tiltorder.append(float(pair[0]))
    accumulated_dose.append(float(pair[1]))

  #print tiltorder, accumulated_dose
  tiltorderfile.close()
  #

  # Reading the output of CTFFIND
  micnames, avgdefoci, defocusv = read_relion_star(outputstarname)

  #print exttilts
  #print avgdefoci
  #print len(exttilts), len(tiltorder)
  

  #sys.exit()

  if len(tiltorder) != len(exttilts):
    print ':: RELION sub-tomogram averaging :: ' + '\n' + 'The number of images in the CTFFIND output file and the tilt order file are different. Exiting'
    sys.exit()

  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'The number of images in the CTFFIND output file and the tilt order file are the same. Continuing.'
  
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Writing out .star files to make 3D CTF volumes ' + '\n' 

  # Pixelsize calculation
  PixelSize = DPixSize/Magnification*10000
  #print PixelSize

  # getting tomogram size using the IMOD program header
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Using IMOD header program to get the size of the tomogram ' + '\n' 
  headerline = 'header -brief -size -input ' + micname 
  print headerline
  status, sizevals = commands.getstatusoutput(headerline)
  tomosize=sizevals.split()
  #print tomosize
  xlimit = float(tomosize[0])
  zlimit = float(tomosize[2])

  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'Writing out .star files to make 3D CTF volumes ' + '\n' 
  subtomonum=0
  for line in coordfile:

    subtomonum = subtomonum+1

    cols = line.split()
    
    # Coordinates of the sub-tomogram in the tomogram
    X = float(cols[0])
    Y = float(cols[1])
    Z = float(cols[2])

    # Output 3D CTF volume and .star file
    outstarname = RelionPartName + MicDirName + MicRootName + '_ctf' + str("%06d" % subtomonum) + '.star' 
    outstarname_for_rec_script = MicRootName + '_ctf' + str("%06d" % subtomonum) + '.star'
    outfile = open(outstarname, 'w')
    outctfname = RelionPartName + MicDirName + MicRootName + '_ctf' + str("%06d" % subtomonum) + '.mrc' 
    outctfname_for_rec_script = MicRootName + '_ctf' + str("%06d" % subtomonum) + '.mrc'
    

    # Writing out the header of the ctf star file    
    outfile.write('data_images' + '\n')
    outfile.write('loop_' + '\n')
    outfile.write('_rlnDefocusU #1 ' + '\n')
    outfile.write('_rlnVoltage #2 ' + '\n') 
    outfile.write('_rlnSphericalAberration #3 ' + '\n')
    outfile.write('_rlnAmplitudeContrast #4 ' + '\n')
    outfile.write('_rlnAngleRot #5 ' + '\n')
    outfile.write('_rlnAngleTilt #6' + '\n')
    outfile.write('_rlnAnglePsi #7 ' + '\n')
    outfile.write('_rlnBfactor #8 ' + '\n')
    outfile.write('_rlnCtfScalefactor #9 ' + '\n')

    for j in range(0,len(exttilts)):
      
      avgdefocus = float(avgdefoci[j])
      tilt_radians = (exttilts[j]*math.pi/180)
      tilt_degrees = exttilts[j]
      #print tilt_radians, tilt_degrees

      xtomo = float(X - (xlimit/2) )*PixelSize
      ztomo = float(Z - (zlimit/2) )*PixelSize
      #print xtomo, ztomo

      # Calculating the height difference of the particle from the tilt axis
      ximg = (xtomo*(math.cos(tilt_radians))) + (ztomo*(math.sin(tilt_radians)))
      deltaD = ximg*math.sin(tilt_radians)
      ptcldefocus = avgdefocus + deltaD
      #print ptcldefocus
      #

      # Now weighting the 3D CTF model using the tilt dependent scale factor and the dose dependent B-Factor
      tiltscale = math.cos(abs(tilt_radians))
      #print tiltscale

      tiltstep = (max(exttilts) - min(exttilts))/(len(exttilts)-1)
      besttiltdiff = tiltstep + 0.5


      for k in range(0,len(tiltorder)):
          
          tiltdiff = abs(tilt_degrees-tiltorder[k])
          
          if tiltdiff < (tiltstep+0.25):
            if tiltdiff < besttiltdiff:
              besttiltdiff = tiltdiff
              accumulated_dose_current = accumulated_dose[k]

      doseweight = accumulated_dose_current * Bfactor
      #print exttilts, tiltorder, accumulated_dose, besttiltdiff, accumulated_dose_current
      #print doseweight
      #

      # Writing parameters in the .star file for each 2D slice of the 3D CTF model volume
      ang_rot = '0.0'
      ang_psi = '0.0'
      ctfline =  str("%.2f" % ptcldefocus) + '\t' + str(Voltage) + '\t' + str(Cs) + '\t' + str(AmpContrast) + '\t' + ang_rot + '\t' + str(tilt_degrees) + '\t' + ang_psi + '\t' + str(doseweight) + '\t' + str("%.2f" % tiltscale) + '\n'
      outfile.write(ctfline)

    # RELION 3D CTF model reconstruction file
    reconstructline = 'relion_reconstruct --i ' + outstarname_for_rec_script + ' --o ' + outctfname_for_rec_script + ' --reconstruct_ctf ' + '$1' + ' --angpix ' + str("%.2f" % PixelSize) + '\n'
    relionfile.write(reconstructline)
    
    # writing the .star file for refinement
    currentsubtomoname = RelionPartName+ MicDirName +  MicRootName + '_' + SubtomoName + str("%06d" % subtomonum) + '.mrc'
    subtomostarline = micname + '\t' + str(X) + '\t' + str(Y) + '\t' + str(Z) + '\t' + currentsubtomoname + '\t' + outctfname + '\n'
    subtomostarfile.write(subtomostarline)

    outfile.close()
  
  relionfile.close()
  ctfreconstmasterfile.write('cd ' + RelionPartName + MicDirName + '\n')
  ctfreconstmasterfile.write( RelionRecFileName_for_script + ' $1\n')
  ctfreconstmasterfile.write('cd ' + ScriptDir + '\n')
  os.chmod(RelionRecFileName, stat.S_IRWXU)
  print ':: RELION sub-tomogram averaging :: ' + '\n' + '.star files to make 3D CTF volumes were written out in ' + RelionRecDir + '\n' 
  print ':: RELION sub-tomogram averaging :: ' + '\n' + 'shell script to reconstruct the 3D CTF volumes is' + RelionRecFileName + '\n'

  ##### Making .star files for each 3D CTF Volume #################

subtomostarfile.close()
ctfreconstmasterfile.close()
print ':: RELION sub-tomogram averaging :: ' 
print 'Please extract sub-tomograms using the RELION GUI. Remember to use the same subtomoname as you gave in this script'
print 'Please run the 3D CTF model volume reconstructions using the .sh scripts written in the working directory' 
print 'run this script from the command line with the command '
print 'do_all_reconstruct_ctfs.sh SubtomogramSize '
print 'STAR file to use for refinement (after sub-tomogram extraction and 3D CTF volume reconstruction) was written in ' + subtomostarname