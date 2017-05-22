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
#       Samuel Wyatt 2017
#       samuel.wyatt@ttu.edu

#Usage:
#       python swift_xrt_analysis.py -a [bool RunsAllNights] --ra [RA] --dec [DEC]
#
#       Arguments: 
#            
#	It will unpack that data if it is a .tar file
#	Sets up a working director /clfiles/
#	loops over each night and extracts the spectrum in cts/sec
#	if there is a background region (bg.reg) to compare to
#		it will extract the background
#		and perform poisson statistics on the file to determine if there is a detection
#	Outputs data in csv and html, along with on the screen

#Requirments:
#       Have xspec and xselect in your pathway
#       You need to already have the Background region in the directory named: bg.reg

#Outputs:
#	csv and html formats
#		Date (yyyy/mm/dd/T) and mjd
#       	Exposure time
#       	Source region cts/s
#       	Background region cts/s
#       	Extracted Source from Background cts/s
#		Poisson Probability
#		Detection Confidence

#Future implementations
#       HEARSAC API to generate flux

global aper
aper = "25\""
global camera_aper_size
camera_aper_size = 0.2

class Spectrum:
	def __init__(self, sC, sR, bC, bR, bS, bSError, pT, e, mjd, dater, prob, confidence):
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
		print "clfiles Directory Does Not Exist"
		print "Making clfiles Directory"
		os.makedirs(cwd+"/clfiles")
	print("")
	for file in dirFiles:
		if "000" in file and ".tar" not in file:
                	cmd = "cp "+file+"/xrt/event/sw"+file+"xpcw3po_cl.evt.gz clfiles/."
                	print cmd
                	subprocess.call(cmd, shell=True)
		if file == "source.reg" or file == "bg.reg":
			cmd = "cp "+file+" clfiles/."
			print cmd
			subprocess.call(cmd, shell=True)
	print("")

def Untar(cwd):
	#Untars the data
        dirlis = os.listdir(cwd)
        tarboys = [x for x in dirlis if ".tar" in x]
        for tar in tarboys:
                cmd = "tar -xvf "+tar
                subprocess.call(cmd, shell=True)
	print("")

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
		
def OutputHtmlTable(data):
	#outputs the data into an htlml file
	#type elinks SpectrumData.html to view it.
	param = "<!DOCTYPE html>\n<html>\n<head>\n<style>\ntable, th, td {\n   border: 1px solid black;\n}\n</style>\n</head>\n<body>"
	param += "<table>"
	
	param += "<tr>\n<th>Date</th>\n<th>MJD</th>\n<th>Exposure</th>\n<th>Source Count</th>\n<th>Source Rate</th>\n<th>BG Count</th>\n<th>BG Rate</th>\n<th>Backed Rate</th>\n<th>Backed Rate Error</th>\n<th>Probability</th>\n<th>Confidence</th>\n"
	for d in data:
		param += "<tr>\n<td>"+str(d.date_readable)+"</td>\n<td>"+str(d.mjd)+"</td>\n<td>"+str(d.exposureTime)+"</td>\n<td>"+str(d.sourceCount)+"</td>\n<td>"+str(d.sourceRate)+"</td>\n<td>"+str(d.bgCount)+"</td>\n<td>"+str(d.bgRate)+"</td>\n<td>"+str(d.backedSpectrum)+"</td>\n<td>"+str(d.backedSpectrumError)+"</td>\n<td>"+str(d.detection_prob)+"</td>\n<td>"+d.confidence+"</td>\n</tr>\n"
		
	param += "</table>\n"
	param += "</body>\n</html>"
	fo = open("SpectrumData.html","w")
	fo.write(param)
	fo.close()	

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

def OutputData(data):
	#outputs the data in to a comma seperated variable file
	fo = open("AllSpectrum.csv", 'w')
	fo.write("SourceCount,SourceRate,BGCount,BGRate,BackedSpectrum,BackedSpectrumError,PercentTotal,ExposureTime,MJD,Date,Probability,Confidence\n")
	for d in data:
		param = ""
		param += str(d.sourceCount)+","
		param += str(d.sourceRate)+","
		param += str(d.bgCount)+","
		param += str(d.bgRate)+","
		param += str(d.backedSpectrum)+","
		param += str(d.backedSpectrumError)+","
		param += str(d.percentTotal)+","
		param += str(d.exposureTime)+","
		param += str(d.mjd)+","
		param += str(d.date_readable)+","
		param += str(d.detection_prob)+","
		param += str(d.confidence)+"\n"
		fo.write(param)
	fo.close();

def ExtractSpectrum(eventFile, sourceRA, sourceDEC):	
	#main function that calls xselect and xspec
	#also performs poisson statistics and determines if there is a detection
	eventhdr = fits.getheader(eventFile)
	mjd = eventhdr['MJD-OBS']
	date = eventhdr['DATE-OBS']
	cwd = os.getcwd()
	dirFiles = os.listdir(cwd)
	extractBG = "bg.reg" in dirFiles
	if "source.pi" in dirFiles or "bg.pi" in dirFiles:
		cmd = 'rm bg.pi source.pi'
        	subprocess.call(cmd, shell=True)


        #creates the xselect script to extract the source spectrum
        make_xselectfile("TESTxspec_script", eventFile, "source")
        cmd = "xselect @xselectSpectrum.xco > xsel_source.dat"
        subprocess.call(cmd, shell=True)
	
        xselDict_source = TextFileToDictionary('xsel_source.dat')
	boy = 1 if "completed" in xselDict_source['47'] else 0

        #output the exposure time
	print("DATE OF OBSERVATION: "+date+"\n")
        print("Exposure sec: "+xselDict_source[str(52+boy)])

        #output the source counts 
        print("Source cts/s (unextracted)")
        print(xselDict_source[str(53+boy)])

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
        	print("BG cts/s")
        	print(xselDict_bg[str(53+girl)])

        	#output the extracted spectrum cts/s
        	print("After extracting the Background from the Source")
        	print(xspecDict['20'])

		splitSpec = [x for x in xspecDict[str(20)].split(' ') if x != '']

		backSpectrum = float(splitSpec[6])
		backSpectrum_error = float(splitSpec[8])
		percentTotal = float(splitSpec[9].split('(')[1]) if len(splitSpec) > 9 else 0.0

        	splitBg = [x for x in xselDict_bg[str(53+girl)].split(' ') if x != '']
		bgCount = float(splitBg[2])
		bgRate = float(splitBg[5])
		
		#HERE DO POISSON STATISTICS ON BG COUNT RATE !BACK AVG
		
		bckAvg = GetBackAvg(bgCount)
		print("Background expected rate for source: "+ str(bckAvg)[0:8])
		poisson_prob, poisson_confidence = PoissonCalc(srcCount, bckAvg)
		
		print("Probability of source: "+ str(poisson_prob))
		print("Confidence of detection: "+ poisson_confidence)
		print("")
	
        	cmd = 'mv bg.pi bg_'+eventFile.split('.')[0]+'.pi'
        	subprocess.call(cmd, shell=True)

        cmd = 'mv source.pi source_'+eventFile.split('.')[0]+'.pi'

        subprocess.call(cmd, shell=True)
	return Spectrum(srcCount, srcRate, bgCount, bgRate, backSpectrum, backSpectrum_error, percentTotal, exposure, mjd, date, poisson_prob, poisson_confidence)

def main():
	#implmented option parsing

	description = "OptionParserTest"
	usage = "python %prog -a [run_all observations] --ra [RA] --dec [DEC]"
	
	parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
	parser.add_option("-a", "--runall", dest="runall", action="store_true", default=False, help="Runs the script on all the observations.")
	parser.add_option("-o", "--run1", dest="run_one", type=str, default="", help="Select the one event file that you want to run")
	parser.add_option("--ra", dest="ra", type=str, default="", help="Object's Right Ascension")
	parser.add_option("--dec", dest="dec", type=str, default="", help="Object's Declinaiton")
	option, args = parser.parse_args()

	_runall = option.runall
	_run1 = option.run_one
	_ra = option.ra
	_dec = option.dec
	
	if (not _runall and _run1 == ""):
		sys.argv.append('--help')
	option, args = parser.parse_args()	
	
	cwd = os.getcwd()
	if Check_Untar(cwd):
		print "Files Need to be untarred"
		Untar(cwd)
	if Check_Clsetup(cwd):
		print "CLDirectory needs to be setup"
		CLSetup(cwd)
	
	RA, DEC = _ra, _dec

        os.chdir(os.getcwd()+"/clfiles")

        dirFiles = os.listdir(cwd+"/clfiles")
	
        clfiles = [x for x in dirFiles if "_cl.evt.gz" in x] if _runall else  [_run1]
	clfiles.sort()
        extractBG = "bg.reg" in dirFiles

        if not extractBG:
                print "There is no background to extract, obtaining raw source spectrum"

        if RA == "" and DEC == "":
                print "WARNING! Grabbing Object's RA and DEC from Event File!"
                print "WCS may not be perfect for each event file!"
                print "Manually input the RA and DEC in the script, or pass as parameters."

                if len(clfiles) > 0:
                        RA, DEC = AverageRADEC(clfiles)
	else:
		cmd = "rm source.reg"
		subprocess.call(cmd, shell=True)
		Region(str(RA), str(DEC), aper, "source")
	if "source.reg" not in dirFiles:
		Region(str(RA), str(DEC), aper, "source")


        print ("RA: "+ RA + " DEC: " + DEC)
        print("")
	data = []
	
        if len(clfiles) > 0:

                #if "bg.reg" not in dirfiles:
                        #Leave uncommented. Generates a random background outside the source not scientifically significant
                        #generate_background_files(source_RA, source_DEC, aper, camera_aper_size)

                for file in clfiles:
                        data.append(ExtractSpectrum(file, RA, DEC))
	
		OutputData(data)
		OutputHtmlTable(data)

if __name__ == "__main__":
	main()
