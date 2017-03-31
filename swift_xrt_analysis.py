import re
import sys
import os
from astropy.io import fits
import math
import random
import subprocess
import random

#Swift XRAY Analysis Pipe
#       Uses xselect and xspec to extract a cts/s
#       Samuel Wyatt 2017
#       samuel.wyatt@ttu.edu

#Usage:
#       python swift_xrt_analysis.py [eventfile]
#
#       Arguments: Eventfile - the swift (level 2) screened event file:
#                  example: sw00034976001xpcw3po_cl.evt

#Requirments:
#       Have xspec and xselect in your pathway
#       If there is already a bg.pi or a source.pi in the directory, it will complain
#       You need to already have the Background region in the directory named: bg.reg

#Outputs:
#       Exposure time
#       Source region cts/s
#       Background region cts/s
#       Extracted Source from Background cts/s

#Future implementations
#       loop over all evt files and find average cts/s
#       HEARSAC API to generate flux
#       develop as an object oriented process because fun

global aper
aper = "18\""
global camera_aper_size
camera_aper_size = 0.2

def TextFileToDictionary(fileName):
        fo = open(fileName, 'r')
        itera = 0
        dict = {str(itera): fileName}
        for line in fo:
                itera = itera + 1
                dict[str(itera)] = line
        fo.close()
        return dict

def create_rand_radec(source_RA, source_DEC, camera_aper_size):
        raDiff = camera_aper_size*random.uniform(0.2, 1)
        decDiff = camera_aper_size*random.uniform(0.2, 1)

        randRA = source_RA + raDiff if random.uniform(0,1) <= 0.5 else source_RA - raDiff
        randDEC = source_DEC + decDiff if random.uniform(0,1) <= 0.5 else source_DEC - decDiff
        return [randRA, randDEC]

def generate_background_files(source_RA, source_DEC, aper, camera_aper_size, iterations=1):
        #DO NOT USE
        #Will generate random background regions around the source based on the number of iterations specified
        for itera in range(0,iterations,1):
                randCoords = create_rand_radec(source_RA, source_DEC, camera_aper_size)
                filename = "bg" if iterations == 1 else "bg"+str(itera)
                Region(str(randCoords[0]), str(randCoords[1]), aper, filename)

def Region(RA, DEC, APER, FILE):
        #Creates a ds9 circle region
        #inputs RA and DEC of source or background center in either hex or dec format.
        #aperature in arcsec around the source
        #FILE is a string, either "source" or "bg"

        color = "red" if FILE == "source" else "green"
        defaultParamLine = "#Region file format: DS9 version 4.1\nglobal color="+color+" dashlist=8 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\ncircle("+RA+","+DEC+","+APER+")"
        fo = open(FILE+".reg", 'w')
        fo.write(defaultParamLine)
        fo.close()

def make_xselectfile(sessionName, eventFile, fileName):
        #Creates the xselect script that extracts a spectrum
        #inputs either the source or background and outputs either accordingly
        # read event                                    // Prompts user for action
        # .                                             // Specifies directory
        # eventfile_cl.tz                               // Specifies the event file     
        # yes                                           // resets the session
        # extract image                                 
        # filter region (either) source.reg (or) bg.reg
        # extract spectrum
        # save spectrum (either) source.pi (or) bg.pi
        # exit                                  
        # no
        region = fileName+".reg"
        saveSpectrum = fileName+".pi"

        paramLine = sessionName+"\nread event\n.\n"+eventFile+"\nyes\nextract image\nfilter region "+region+"\nextract spectrum\nsave spectrum "+saveSpectrum+"\nexit\nno"
        fo = open("xselectSpectrum.xco", 'w')
        fo.write(paramLine)
        fo.close()

def make_xspecfile(itera=0):
        #Creates the the xspec script that extracts the background from the source
        # data source.pi
        # back bg.pi

        bgfile = "bg.pi"
        if itera is not 0:
                bgfile = "bg"+str(itera)+".pi"
        paramLine = "data source.pi\nback "+bgfile
        fo = open("xspecFile.xco", 'w')
        fo.write(paramLine)
        fo.close();

def main():
        #############
        #STARTS HERE#
        #############

        #The only argument should be the event fits file
        args = sys.argv
        eventFile = args[1]

        #extract the ra and dec from the fits header
        eventFile_fits = fits.open(eventFile)
        source_RA = eventFile_fits[0].header['RA_OBJ']
        source_DEC = eventFile_fits[0].header['DEC_OBJ']
        print("RA: " + str(source_RA) + ", DEC: " + str(source_DEC))

        #if you think the wcs is inaccurate with the fits file, you can define them here
        #source_RA = 63.031750
        #source_DEC = -32.853036

        #creates a source region in ds9.reg format.
        #saves as source.reg

        ######################################
        #I would manually create the source region if you don't trust the wcs in the fits image#
        #####################################

        Region(str(source_RA), str(source_DEC), aper, "source")

        #Leave uncommented. Generates a random background outside the source not scientifically significant
        #generate_background_files(source_RA, source_DEC, aper, camera_aper_size)        
	#creates the xselect script to extract the source spectrum
        make_xselectfile("TESTxspec_script", eventFile, "source")
        cmd = "xselect @xselectSpectrum.xco > xsel_source.dat"
        subprocess.call(cmd, shell=True)

        #creates the xselect script to extract the background spectrum
        make_xselectfile("TESTxspec_script", eventFile, "bg")
        cmd = "xselect @xselectSpectrum.xco > xsel_bg.dat"
        subprocess.call(cmd, shell=True)

        #crates the xspec script to remove the background spectrum from the source spectrum
        make_xspecfile()
        cmd = "xspec xspecFile.xco > xspecoutput.dat"
        subprocess.call(cmd, shell=True)

        #Each step outputs a .dat file composed what is output from the xselect commands and xspec command
        #The ex.TextFileToDictionary("Filename") creates a dict[line number][line content]
        xselDict_source = TextFileToDictionary('xsel_source.dat')
        xselDict_bg = TextFileToDictionary('xsel_bg.dat')
        xspecDict = TextFileToDictionary('xspecoutput.dat')

        #output the exposure time
        print("Exposure sec: "+xselDict_source['52'])
        print(xselDict_source['53'])
        #output the source counts 
        print("Source cts/s (unextracted)")
        print(xselDict_source['54'])
        print(xselDict_source['55'])
        print("")
        #outpt the bg counts
        print("BG cts/s")
        print(xselDict_bg['54'])
        print(xselDict_bg['55'])
        print("")
        #output the extracted spectrum cts/s
        print("After extracting the Background from the Source")
        print(xspecDict['20'])

        cmd = 'mv source.pi source_'+eventFile.split('.')[0]+'.pi'
        subprocess.call(cmd, shell=True)
        cmd = 'mv bg.pi bg_'+eventFile.split('.')[0]+'.pi'
        subprocess.call(cmd, shell=True)

if __name__ == "__main__":
        main()
