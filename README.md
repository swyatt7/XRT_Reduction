# XRT_Reduction
Quick Script that extracts an x-ray spectrum from a swift xrt event
Swift XRAY Analysis Pipe
       Uses xselect and xspec to extract a cts/s
       Samuel Wyatt 2017
       samuel.wyatt@ttu.edu

Usage:
       python swift_xrt_analysis.py [eventfile]

       Arguments: Eventfile - the swift (level 2) screened event file:
                  example: sw00034976001xpcw3po_cl.evt

Requirments:
       Have xspec and xselect in your pathway
       If there is already a bg.pi or a source.pi in the directory, it will complain
       You need to already have the Background region in the directory named: bg.reg

Outputs:
       Exposure time
       Source region cts/s
       Background region cts/s
       Extracted Source from Background cts/s

Future implementations
       loop over all evt files and find average cts/s
       HEARSAC API to generate flux
       develop as an object oriented process because fun
