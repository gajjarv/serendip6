# serendip6
Set of python codes for SERENDIP 6 ETFITS file analysis 

FitsBeamRFI
v1 : Basic version. Background copied from analyedFITS.py (Kyle's code). Just made few modification to plot each beam with respect to each other. 

v2 : Input output of fits file is done with the pyFits, while in the v1 fitsio was used. So structure of the entire program is changed here (active) 

v3 : With histograms for each beam as stacked histograms to identify RFI lines in the histograms (dropped and not done)
v4 : copied from v2, All beams are in a single panel but with different colour 

SimulFits.py
This prog simulates ETHITS file that I generate using a sample input file. It creates similar file with same length and other parameters.
v1 : Basic version with lots of trials, this version is for backup. 
v2 : More advance vesion with user defined signal and changes 
v3 : This version has posibility of adding signal in just one beam (Not finished)

analyseFits
v1: Basic version I got
v2 : Trial with few changes 
v3 : With color histogram and waterfall diagram and possiblity of bin2frequency, however needs more works to implement
     Here color represent SNR
v4 : With color histogram and waterfall diagram
     Here color represent Fine Chan number
v5 : Copied from v4, this version has choice of SNR range selection to be plotted in the data
v6 : Copied from v5, along with the hits histogram, It also plots expected noise distribution from a given number of hits and number of channels 

compBeam.py

v1 : Basic version but there a bug, just copied it for backup
v2 : Working version. It removes hits present in more than one beam. Although, fits headers are not yet fixed
v3 : This version does RFI removal from near by frequency channels across all beams (working version it does remove things which are 
present near to each other across different beams)
v4 : This version does RFI removal in time and fixed corase channel which are known to show fixed RFI lines (done but does not remove in time)
v5 : This version will remove RFI lines if they present at different RA and DEC level (test version but not used)
v6 : copied from v5, This version will remove RFI averaging over time (failed)
v7 : again copied from v5 as the completely different logic was tried. This is the real version which works very well. (active)
v8 : Same as v7 but only things now it needs to do is to take command line input for length of the time and frequency window  
