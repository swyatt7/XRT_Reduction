import re
import sys
import os
from astropy.io import fits
import math
import random
import subprocess
import random
import numpy
from optparse import OptionParser, OptionGroup

#Swift XRAY Analysis Pipe
#       Uses xselect and xspec to extract a cts/s
#       Uses Ftools nH and pimms to extract a flux and luminosity
#       Samuel Wyatt 2017
#       samuel.wyatt@ttu.edu

#Usage:
#       python swift_xrt_analysis.py -a [bool RunsAllNights] --ra [RA] --dec [DEC]
#
#       Arguments: 
#            
#       It will unpack that data if it is a .tar file
#       Sets up a working director /clfiles/
#       loops over each night and extracts the spectrum in cts/sec
#       if there is a background region (bg.reg) to compare to
#               it will extract the background
#               and perform poisson statistics on the file to determine if there is a detection
#       It will convert the cts/sec into Flux
#       if there is a distance (mpc) then it will calculate a luminosity
#       Outputs data in csv and html, along with on the screen
#

#       python swift_xrt_analysis.py -a -n [supernova name] --infopath [/info/path/SN_info.dat]
#       You can have a file called SN_info.dat that has the structure: name,ra,dec,distance
#       It will match the supernova name (can be portions of the name)
#       So you don't have to manually enter the ra, dec, distance everytime.

#Requirments:
#       Have xspec and xselect in your pathway
#       You need to already have the Background region in the directory named: bg.reg
#       In order to calculate flux and luminosity:
#               You will need to have the program pimms
#               you can input the pathway to the program:
#               Download and instructions located here.
#                       https://heasarc.gsfc.nasa.gov/docs/software/tools/pimms.html#latest
#       If you don't have access to the program, there is a web tool:
#               https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
#               converts from (swift/xrt/pc) to flux (ergs/cec)


#Outputs:
#       csv and html formats
#               Date (yyyy/mm/dd/T) and mjd
#               Exposure time
#               Source region cts/s
#               Background region cts/s
#               Extracted Source from Background cts/s
#               Flux and Luminosity
#               Poisson Probability
#               Detection Confidence
#       Flux_light curve graph.

global CanGraph 
CanGraph = True
try:
	import matplotlib.pyplot as plt
except ImportError:
	CanGraph = False
	pass

global aper
aper = "25\""
global camera_aper_size
camera_aper_size = 0.2

def setGlobalVerbose(v):
	global verbose
	verbose = v

class Spectrum:
	def __init__(self, sC, sR, bC, bR, bS, bSError, pT, e, mjd, dater, prob, confidence, flux, lum):
		self.sourceCount = sC
		self.sourceRate = sR
		self.bgCount = bC
		self.bgRate = bR
		self.backedSpectrum = bS
		self.backedSpectrumError = bSError
		self.percentTotal = pT
		self.exposureTime = e
		self.mjd = mjd
		self.date_readable = dater
		self.detection_prob = prob
		self.confidence = confidence
		self.flux = flux
		self.luminosity = lum

class SNInfo:
        def __init__(self, name, ra, dec, dist):
                self.name = name
                self.ra = ra
                self.dec = dec
                self.dist = dist

def printv(text, v):
	if v:
		print(text)

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

def MakePimmsCommand_CountsToFlux(rate, nH, energy_range, powerlaw):
        parma = "model pl "+powerlaw+" "+nH+"\nfrom swift xrt pc "+energy_range+"\ninst flux ergs "+energy_range+"\ngo "+str(rate)+"\nexit\n"
        fo = open("pimms_countstoflux.xco", 'w')
        fo.write(parma)
        fo.close();

def CLSetup(cwd):
	#after the data has been untarred
	#it will create the clfiles directory
	#and grab all the necessary swift cleaned files from the data
	#it also checks if there are any pre-made region files in the directory and moves them
	dirFiles = os.listdir(cwd)

	if Check_Untar(cwd):
		Untar(cwd)
	
	dirFiles = os.listdir(cwd)
	if not os.path.exists(cwd+"/clfiles"):
		printv("clfiles Directory Does Not Exist", verbose)
		printv("Making clfiles Directory", verbose)
		os.makedirs(cwd+"/clfiles")
	printv("", verbose)
	for file in dirFiles:
		if "000" in file and ".tar" not in file:
                	cmd = "cp "+file+"/xrt/event/sw"+file+"xpcw3po_cl.evt.gz clfiles/."
                	printv(cmd, verbose)
                	subprocess.call(cmd, shell=True)
		if file == "source.reg" or file == "bg.reg":
			cmd = "cp "+file+" clfiles/."
			printv(cmd, verbose)
			subprocess.call(cmd, shell=True)
	printv("", verbose)

def Untar(cwd):
	#Untars the data
        dirlis = os.listdir(cwd)
        tarboys = [x for x in dirlis if ".tar" in x]
        for tar in tarboys:
                cmd = "tar -xvf "+tar
                subprocess.call(cmd, shell=True)
	printv("", verbose)

def Check_Untar(cwd):
	#Checks the current directory and determines whether or not the data needs to be untarred
	#tests if there are more tar files than other files (expected to be in the directory)
        dirFiles = os.listdir(cwd)
        tarOrReg = [x for x in dirFiles if (".tar" in x or ".reg" in x) and not "clfiles" in x]
        if len(tarOrReg) == len(dirFiles) - 1:
                return True
        return False

def Check_Clsetup(cwd):
	#Checks whether or not to create the clfiles directory
        dirFiles = os.listdir(cwd)
        if "clfiles" in dirFiles:
                return False
        return True
	
def AverageRADEC(clfiles):
	#last ditch effort to get the RA and DEC
	#grabs the pointing of the camera. Will not get the source
	#mainly used for testing
	RA, DEC = [], []
	for file in clfiles:
		eventHdr = fits.getheader(file)
		RA.append(float(eventHdr['RA_OBJ']))
		DEC.append(float(eventHdr['DEC_OBJ']))

	return str(numpy.mean(RA))[0:10], str(numpy.mean(DEC))[0:10] 
		
def OutputHtmlTable(data, name):
	#outputs the data into an htlml file
	#type elinks SpectrumData.html to view it.
	filename = name + "_SpectrumData.html" if name is not "" else "SpectrumData.html"
	param = "<!DOCTYPE html>\n<html>\n<head>\n<style>\ntable, th, td {\n   border: 1px solid black;\n}\n</style>\n</head>\n<body>"
	param += "<table>"
	
	param += "<tr>\n<th>Date</th>\n<th>MJD</th>\n<th>Exposure</th>\n<th>Source Count</th>\n<th>Source Rate</th>\n<th>BG Count</th>\n<th>BG Rate</th>\n<th>Backed Rate</th>\n<th>Backed Rate Error</th>\n<th>Flux</th>\n<th>Luminosity</th>\n<th>Probability</th>\n<th>Confidence</th>\n"
	for d in data:
		param += "<tr>\n<td>"+str(d.date_readable)+"</td>\n<td>"+str(d.mjd)+"</td>\n<td>"+str(d.exposureTime)+"</td>\n<td>"+str(d.sourceCount)+"</td>\n<td>"+str(d.sourceRate)+"</td>\n<td>"+str(d.bgCount)+"</td>\n<td>"+str(d.bgRate)+"</td>\n<td>"+str(d.backedSpectrum)+"</td>\n<td>"+str(d.backedSpectrumError)+"</td>\n<td>"+str(d.flux)+"</td>\n<td>"+str(d.luminosity)+"</td>\n<td>"+str(d.detection_prob)+"</td>\n<td>"+d.confidence+"</td>\n</tr>\n"
		
	param += "</table>\n"
	param += "</body>\n</html>"
	fo = open(filename,"w")
	fo.write(param)
	fo.close()	
        cmd = "mv " + filename + " ../."
        subprocess.call(cmd, shell=True)

def PoissonCalc(srcCnts, bckAvg):
	#calculates the poisson statistcs of the detection
	#p(x) = (bckAvg^(srcCnts)/(srcCnts!))*exp(-bckAvg)
	#returns the probability and the sigma detection
	if srcCnts == 0:
		return 0.0, "Non Detection"

	prob = (math.pow(bckAvg, srcCnts) / math.factorial(srcCnts))*math.exp(-bckAvg)
	conf_sub = 1-prob
	confidence = ""
	if conf_sub >= 0.999999998027:
		confidence = "6 sigma"
	elif conf_sub >= 0.999999426697 and conf_sub < 0.999999998027:
		confidence = "5 sigma"
	elif conf_sub >= 0.999936657516 and conf_sub < 0.999999426697:
		confidence = "4 sigma" 
	elif conf_sub >= 0.997300203937 and conf_sub < 0.999936657516:
		confidence = "3 sigma"
	elif conf_sub >= 0.954499736104 and conf_sub < 0.997300203937:
		confidence = "2 sigma"
	elif conf_sub >= 0.682689492137 and conf_sub < 0.954499736104:
		confidence = "1 sigma"
	else: 
		confidence = "Non detection"
	return prob, confidence

def GetBackAvg(bckCnts):
	#returns the expected background average
	#tests whether the radius is in arcmin or arcsec
	#returns bckCnts * (source_aper_area/background_aper_area)
	fi = open("bg.reg", 'r')
	data = []
        for line in fi:
                data.append(line)
        radius = data[3].split(',')[2].split(')')[0]
        if len(radius.split("'")) > 1:
                radius = str(float(radius.split("'")[0])*60.0)
        if len(radius.split('"')) > 1:
                radius = radius.split('"')[0]
	return bckCnts*((math.pow(float(aper.split('"')[0]), 2))/ (math.pow(float(radius), 2)))

def OutputData(data, name):
	#outputs the data in to a comma seperated variable file
	filename = name+"_SpectrumData.csv" if name is not "" else "SpectrumData.csv"
	fo = open(filename, 'w')
	fo.write("SourceCount,SourceRate,BGCount,BGRate,BackedSpectrum,BackedSpectrumError,Flux,Luminosity,PercentTotal,ExposureTime,MJD,Date,Probability,Confidence\n")
	for d in data:
		param = ""
		param += str(d.sourceCount)+","
		param += str(d.sourceRate)+","
		param += str(d.bgCount)+","
		param += str(d.bgRate)+","
		param += str(d.backedSpectrum)+","
		param += str(d.backedSpectrumError)+","
		param += str(d.flux)+","
		param += str(d.luminosity)+","
		param += str(d.percentTotal)+","
		param += str(d.exposureTime)+","
		param += str(d.mjd)+","
		param += str(d.date_readable)+","
		param += str(d.detection_prob)+","
		param += str(d.confidence)+"\n"
		fo.write(param)
	fo.close();
	cmd = "mv " + filename + " ../."
	subprocess.call(cmd, shell=True)

def SNInfoFind(inname, file="/home/swyatt/SN_info.dat"):
        inname = inname.lower()
        if os.path.exists(file):
                data = []
                info = open(file, "r")
                for line in info:
                        sl = line.split(',')
                        dname = sl[0] if len(sl) > 0 else ""
                        ra = sl[1] if len(sl) > 1 else ""
                        dec = sl[2] if len(sl) > 2 else ""
                        dist = sl[3] if len(sl) > 3 else ""
                        name = dname.lower()
                        if inname in name or inname is name or name in inname:
                                return SNInfo(dname, ra, dec, dist)
		info.close()
        printv("SN Info not found, input parameters manually\n", True)
        return SNInfor("", "", "", "")

def FormatCoordforNH(coord):
        #if ra or dec is in deg
        #it needs to formatted as "hh mm ss.s"
        if ":" in coord:
                split_coord = coord.split(":")
                hh = split_coord[0] if len(split_coord) > 0 else "00"
                mm = split_coord[1] if len(split_coord) > 1 else "00"
                ss = split_coord[2] if len(split_coord) > 2 else "00.0"
                return (hh + " " + mm + " " + ss)
        #if ra or dec is in dec, it needs to have a decimal place (just in case)
        elif "." not in coord:
                coord += ".00"
        return coord

def GetNH(RA, DEC):
	#format ra and dec so that ftools nh command will work
	#call ftools command and spit into a file
	#read the file and find the line number with the correct value for column density
	#maybe parse through the the output ra and decs and get the best likely
	#start off with finding the average
	#return nH as a string
	#nh equinox=2000 ra="21 44 41.20" dec="38 19 18" > test
        RA = FormatCoordforNH(RA)
        DEC = FormatCoordforNH(DEC)

        cmd = "nh equinox=2000 ra\"=" + RA + "\" dec=\"" + DEC + "\" > nH.dat"
        subprocess.call(cmd, shell=True)
        nH_dict = TextFileToDictionary("nH.dat")
        return nH_dict[str((len(nH_dict)-2))].split(" ")[8].split("\n")[0]

def GetFlux(rate, nH, energy_range, powerlaw, pathToPimms):
	#nH called in main, and passed into Extract function
	#if pimms exists /usr/local/pimms/pimms, continue, else return nothing
	#create an .xco file and and call the pims command > file
	#find the line number in the output file with the correct value for the flux.
	#maybe use it with nH value and non nH value to get absorbed and unabsorbed
	#return desired valueMakePimmsCommand_CountsToFlux(rate, nH, energy_range, powerlaw)
	MakePimmsCommand_CountsToFlux(rate, nH, energy_range, powerlaw)
	cmd = pathToPimms + " @pimms_countstoflux.xco > pimms_ctf.dat"
	subprocess.call(cmd, shell=True)
	ctf_dict = TextFileToDictionary("pimms_ctf.dat")
	fluxline =  [x for x in ctf_dict["21"].split(" ") if x is not ""]
	return float(fluxline[len(fluxline)-2])

def GetLuminosity(flux, distance):
        #assumes distance is in mpc
        #converts to cm
        #returns L = (4piR^2)*f
        distance = 3.86E24 * distance
        return 4*math.pi*distance*distance*flux
def MakeGraphs(data, name, energy):
	mjd = []
	lum = []
	for d in data:
		mjd.append(float(d.mjd))
		lum.append(float(d.luminosity))
	plt.plot(mjd, lum, 'r+')
	#plt.axis([mjd.min() - 5, mjd.max() + 5, lum.min(), lum.max()])
	plt.xlabel('MJD (days)', fontsize=14)
	plt.ylabel(r'$L_{'+energy+r'}(\frac{erg}{s})$', fontsize=14)
	plt.savefig("flux_curve.png")
	cmd = "mv flux_curve.png ../."
	subprocess.call(cmd, shell=True)

def ExtractSpectrum(eventFile, sourceRA, sourceDEC, nH, energy_range, powerlaw, pathToPimms, distance, extractBG, pimmsExists):	

	#main function that calls xselect and xspec
	#also performs poisson statistics and determines if there is a detection
	eventhdr = fits.getheader(eventFile)
	mjd = eventhdr['MJD-OBS']
	date = eventhdr['DATE-OBS']

        #creates the xselect script to extract the source spectrum
        make_xselectfile("TESTxspec_script", eventFile, "source")
        cmd = "xselect @xselectSpectrum.xco > xsel_source.dat"
        subprocess.call(cmd, shell=True)
	
        xselDict_source = TextFileToDictionary('xsel_source.dat')
	boy = 1 if "completed" in xselDict_source['47'] else 0

        #output the exposure time
	printv("DATE OF OBSERVATION: "+date+"\n", verbose)
        printv("Exposure sec: "+xselDict_source[str(52+boy)], verbose)

        #output the source counts 
        printv("Source cts/s (unextracted)", verbose)
        printv(xselDict_source[str(53+boy)], verbose)

	splitSource = [x for x in xselDict_source[str(52+boy)].split(' ') if x != '']
	exposure = float(splitSource[1])

	splitSource = [x for x in xselDict_source[str(53+boy)].split(' ') if x != '']
	srcCount = float(splitSource[2])
	srcRate = float(splitSource[5])
	
	backSpectrum_error = 0.0
	backSpectrum = 0.0
	bgCount = 0
	bgRate = 0.0
	percentTotal = 0.0
	poisson_prob = 0.0
	poisson_confidence = "N/A"
	flux = GetFlux(srcRate, nH, energy_range, powerlaw, pathToPimms) if pimmsExists else 0.0
	
	if extractBG:
        	#creates the xselect script to extract the background spectrum
        	make_xselectfile("TESTxspec_script", eventFile, "bg")
        	cmd = "xselect @xselectSpectrum.xco > xsel_bg.dat"
        	subprocess.call(cmd, shell=True)
		
        	#crates the xspec script to remove the background spectrum from the source spectrum
        	make_xspecfile()
        	cmd = "xspec xspecFile.xco > xspecoutput.dat"
        	subprocess.call(cmd, shell=True)

        	xselDict_bg = TextFileToDictionary('xsel_bg.dat')
		girl = 1 if "completed" in xselDict_bg['47'] else 0

        	xspecDict = TextFileToDictionary('xspecoutput.dat')

        	#outpt the bg counts
        	printv("BG cts/s", verbose)
        	printv(xselDict_bg[str(53+girl)], verbose)

        	#output the extracted spectrum cts/s
        	printv("After extracting the Background from the Source", verbose)
        	printv(xspecDict['20'], verbose)

		splitSpec = [x for x in xspecDict[str(20)].split(' ') if x != '']

		backSpectrum = float(splitSpec[6])
		backSpectrum_error = float(splitSpec[8])
		percentTotal = float(splitSpec[9].split('(')[1]) if len(splitSpec) > 9 else 0.0
		#if there is a backed rate, then get the flux that corresponds to that
		flux = GetFlux(backSpectrum, nH, energy_range, powerlaw, pathToPimms) if pimmsExists else 0.0
        	splitBg = [x for x in xselDict_bg[str(53+girl)].split(' ') if x != '']
		bgCount = float(splitBg[2])
		bgRate = float(splitBg[5])
		
		#HERE DO POISSON STATISTICS ON BG COUNT RATE !BACK AVG
		
		bckAvg = GetBackAvg(bgCount)
		printv("Background expected rate for source: "+ str(bckAvg)[0:8], verbose)
		poisson_prob, poisson_confidence = PoissonCalc(srcCount, bckAvg)
		
		printv("Probability of source: "+ str(poisson_prob), verbose)
		printv("Confidence of detection: "+ poisson_confidence, verbose)
		printv("", verbose)
	
        	cmd = 'mv bg.pi bg_'+eventFile.split('.')[0]+'.pi'
        	subprocess.call(cmd, shell=True)

	luminosity = GetLuminosity(flux, float(distance)) if distance is not "" else 0.0

	if pimmsExists:	
		printv("Flux (erg/cm/cm/s): " + str(flux), verbose)
		printv("Luminosity (erg/s): " + str(luminosity), verbose)
		printv("", verbose) 

        cmd = 'mv source.pi source_'+eventFile.split('.')[0]+'.pi'
        subprocess.call(cmd, shell=True)

	return Spectrum(srcCount, srcRate, bgCount, bgRate, backSpectrum, backSpectrum_error, percentTotal, exposure, mjd, date, poisson_prob, poisson_confidence, flux, luminosity)

def main():
	#implmented option parsing

	description = "OptionParserTest"
	usage = "python %prog -a [run_all observations] --ra [RA] --dec [DEC]"
	
	parser = OptionParser(usage=usage, description=description, version="%prog 2.0")
	parser.add_option("-a", "--runall", dest="runall", action="store_true", default=False, help="Runs the script on all the observations.")
	parser.add_option("-o", "--run1", dest="run_one", type=str, default="", help="Select the one event file that you want to run")
	parser.add_option("--ra", dest="ra", type=str, default="", help="Object's Right Ascension")
	parser.add_option("--dec", dest="dec", type=str, default="", help="Object's Declinaiton")
	parser.add_option("-p", "--pimms", dest="pathToPimms", type=str, default="/usr/local/pimms/pimms", help="Pathway to pimms commands for flux calculations. Defaults to \"/usr/local/pimms/pimms\"")
	parser.add_option("--pl", dest="powerlaw", type=str, default="2", help="Powerlaw for flux calculations. Defaults to 2")
	parser.add_option("--el", dest="energylevel", type=str, default="0.3-10", help="Energy level range for flux calculations. Defaults to 0.3-10")
	parser.add_option("--dist", dest="distance", type=str, default="", help="Distance to the SN in mpc")
	parser.add_option("-n", "--name", dest="name", type=str, default="", help="Finds a SN_info.dat and matches the name to the Supernovae to get the ra, dec, and distance so it doesn't need to input everytime")
	parser.add_option("--infopath", dest="infopath", type=str, default="/home/swyatt/SN_info.dat", help="Pathway to the SN_info.dat file")
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose option")

	option, args = parser.parse_args()
	_runall = option.runall
	_run1 = option.run_one
	_ra = option.ra
	_dec = option.dec
	_pathToPimms = option.pathToPimms
	_powerlaw = option.powerlaw
	_energyrange = option.energylevel
	_distance = option.distance
	_name = option.name
	_infopath = option.infopath
	setGlobalVerbose(option.verbose)
	
	data = []
	if (_ra == "" and _dec == "" and _name is not ""):
		sninfo = SNInfoFind(_name, _infopath)
		_ra = sninfo.ra
		_dec = sninfo.dec
		_distance = sninfo.dist
		_name = sninfo.name

	if (not _runall and _run1 == ""):
		sys.argv.append('--help')
	option, args = parser.parse_args()	
	
	cwd = os.getcwd()
	if Check_Untar(cwd):
		printv("Files Need to be untarred\n", verbose)
		Untar(cwd)
	if Check_Clsetup(cwd):
		printv("CLDirectory needs to be setup\n", verbose)
		CLSetup(cwd)
	
	pimmsExists = os.path.exists(_pathToPimms)
	if not pimmsExists:
		printv("WARNING", True)
		printv("Pathway to pimms command does not exist!", True)
                printv("Unable to perform Flux calculations", True)
                printv("Make sure to have the correct pimms pathway passed as a param. It defaults to \"/usr/local/pimms/pimms\"", True)
                printv("Pass the parameter: -p [--pimms] /correct/pimms/path/pimms\n", True)
		printv("If pimms is not accessible, go to \"https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl\", and follow the instructions to get the flux\n", True)
	
	RA, DEC = _ra, _dec

        os.chdir(os.getcwd()+"/clfiles")

        dirFiles = os.listdir(cwd+"/clfiles")
	
        clfiles = [x for x in dirFiles if "_cl.evt.gz" in x] if _runall else  [_run1]
	clfiles.sort()
        extractBG = "bg.reg" in dirFiles

        if not extractBG:
                printv("There is no background to extract, obtaining raw source spectrum\n", True)

        if RA == "" and DEC == "":
                printv("WARNING! Grabbing Object's RA and DEC from Event File!", True)
                printv("WCS may not be perfect for each event file!", True)
                printv("Manually input the RA and DEC in the script, or pass as parameters.\n", True)

                if len(clfiles) > 0:
                        RA, DEC = AverageRADEC(clfiles)
	else:
		if "source.reg" in dirFiles:
			cmd = "rm source.reg"
			subprocess.call(cmd, shell=True)

		Region(str(RA), str(DEC), aper, "source")

	if "source.reg" not in dirFiles:
		Region(str(RA), str(DEC), aper, "source")


        printv("RA: "+ RA + " DEC: " + DEC, verbose)
        printv("", verbose)

	nH = GetNH(RA, DEC)

	if pimmsExists:
                printv("Parameters for Flux Calculations", verbose)
                printv("Power Law: " + _powerlaw, verbose)
                printv("Energy Range [keV]: " + _energyrange, verbose)
                printv("nH Column Density [cm^(-2)]: " + str(nH) + "\n", verbose)
	
        if len(clfiles) > 0:

                #if "bg.reg" not in dirfiles:
                        #Leave uncommented. Generates a random background outside the source not scientifically significant
                        #generate_background_files(RA, DEC, aper, camera_aper_size)

                for file in clfiles:
                        data.append(ExtractSpectrum(file, RA, DEC, nH, _energyrange, _powerlaw, _pathToPimms, _distance, extractBG, pimmsExists))
	
		OutputData(data, _name)
		OutputHtmlTable(data, _name)
		if CanGraph:
			MakeGraphs(data, _name, _energyrange)
		printv("type \"elinks SpectrumData.html\" to view the output html table", verbose)	

def CheckFTools():
	#I don't really like how I wrote this
	#But it works
	#Checks if the terminal can use xselect or xspec
	fo = open("test.xco", "w")
	fo.write("test\nexit\nno")
	fo.close()
	cmd = "xselect @test.xco > test.xcm"
	test = subprocess.call(cmd, shell=True)
	subprocess.call("rm test.xc*", shell=True)
	if test is 127:
        	printv("Type 'heainit' into the terminal so that you can access the xselect and xspec commands used in this script", True)
		return False
	return True
	
if __name__ == "__main__":
	if CheckFTools():	
		main()
