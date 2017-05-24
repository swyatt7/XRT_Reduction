# XRT_Reduction
Quick Script that extracts an x-ray spectrum from a swift xrt event

Swift XRAY Analysis Pipe
       Uses xselect and xspec to extract a cts/s
       Uses Ftools nH and pimms to extract a flux and luminosity
       Samuel Wyatt 2017
       samuel.wyatt@ttu.edu

Usage:
       python swift_xrt_analysis.py -a [bool RunsAllNights] --ra [RA] --dec [DEC]

       Arguments: 
            
       It will unpack that data if it is a .tar file
       Sets up a working director /clfiles/
       loops over each night and extracts the spectrum in cts/sec
       if there is a background region (bg.reg) to compare to
               it will extract the background
               and perform poisson statistics on the file to determine if there is a detection
       It will convert the cts/sec into Flux
       if there is a distance (mpc) then it will calculate a luminosity
       Outputs data in csv and html, along with on the screen


       python swift_xrt_analysis.py -a -n [supernova name] --infopath [/info/path/SN_info.dat]
       You can have a file called SN_info.dat that has the structure: name,ra,dec,distance
       It will match the supernova name (can be portions of the name)
       So you don't have to manually enter the ra, dec, distance everytime.

Requirments:
       Have xspec and xselect in your pathway
       You need to already have the Background region in the directory named: bg.reg
       In order to calculate flux and luminosity:
               You will need to have the program pimms
               you can input the pathway to the program:
               Download and instructions located here.
                       https://heasarc.gsfc.nasa.gov/docs/software/tools/pimms.html#latest
       If you don't have access to the program, there is a web tool:
               https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
               converts from (swift/xrt/pc) to flux (ergs/cec)


Outputs:
       csv and html formats
               Date (yyyy/mm/dd/T) and mjd
               Exposure time
               Source region cts/s
               Background region cts/s
               Extracted Source from Background cts/s
               Flux and Luminosity
               Poisson Probability
               Detection Confidence
       Flux_light curve graph.
