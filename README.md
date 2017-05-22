# XRT_Reduction
Quick Script that extracts an x-ray spectrum from a swift xrt event

Swift XRAY Analysis Pipe
       Uses xselect and xspec to extract a cts/s and detection confidence
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
       Outputs data in csv and html, along with on the screen

Requirments:
       Have xspec and xselect in your pathway
       You need to already have the Background region in the directory named: bg.reg

Outputs:
       csv and html formats
               Date (yyyy/mm/dd/T) and mjd
               Exposure time
               Source region cts/s
               Background region cts/s
               Extracted Source from Background cts/s
               Poisson Probability
               Detection Confidence

Future implementations
       HEARSAC API to generate flux

