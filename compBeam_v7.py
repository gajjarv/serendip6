import pyfits
from astropy.io import fits
import os,sys
import numpy as np
from itertools import *
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import operator

# 0.8 Hz resolution so frequency window is NearFreq*0.8
NearFreq = 50000 

#Time window is base_event - BL : base_event + BL 
BL = 5

def update_progress(progress):
            barLength = 10 # Modify this to change the length of the progress bar
            status = ""
            if isinstance(progress, int):
                progress = float(progress)
            if not isinstance(progress, float):
                progress = 0
                status = "error: progress var must be float\r\n"
            if progress < 0:
                progress = 0
                status = "Halt...\r\n"
            if progress >= 1:
                progress = 1
                status = "Done...\r\n"
            block = int(round(barLength*progress))
            text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
            sys.stdout.write(text)
            sys.stdout.flush()


def RADEC(RA,DEC):
#This function detects RFI lines across wide range of RA and DEC from the data
	#print zip(RA,DEC)
	#plt.plot(RA,DEC)	
	#plt.show()
	c = SkyCoord(ra=RA,dec=DEC,unit='deg',frame='icrs')
	minra = min(enumerate(c.ra.deg),key=operator.itemgetter(1))
	maxra = max(enumerate(c.ra.deg),key=operator.itemgetter(1))
	mindec = min(enumerate(c.dec.deg),key=operator.itemgetter(1))
	maxdec = max(enumerate(c.dec.deg),key=operator.itemgetter(1))

	sep = max(enumerate([maxra[1]-minra[1],maxdec[1]-mindec[1]]),key=operator.itemgetter(1))
	if(sep):
		print "Max difference in  RA"
		print minra[0],maxra[0]
	else: 
		print "Min difference in Dec"
		print mindec[0],maxdec[0]
	#print enumerate([maxra[1]-minra[1],maxdec[1]-mindec[1]]),sep

	#print minra[1],maxra[1],mindec[1],maxdec[1]

	#sep = []
	#for b in combinations(range(len(c)),2):
		#print c[b[0]].to_string('hmsdms'),c[b[1]].to_string('hmsdms')
		#sep.append(c[b[0]].separation(c[b[1]]))
	
	#print "Maximum motion in sky : %f" % (max(sep).arcmin)	

def compBeam(BEAM,NS,FS,coarseID):

	#Compare all beams and extract bin number which occurs in all beams
	#commbeam = list(set(BEAM[0]+BEAM[1])&set(BEAM[2]+BEAM[3])&set(BEAM[4]+BEAM[5])&set(BEAM[6]+BEAM[7])&set(BEAM[8]+BEAM[9])&set(BEAM[10]+BEAM[11])&set(BEAM[12]+BEAM[13]))
	
	commbeam1 = []

	#Compare all beams and extract hits which occurs in more than one beam
	#for i in combinations(range(0,14,2),2):
		#print i[0],i[0]+1,i[1],i[1]+1,list(set(BEAM[i[0]]+BEAM[i[0]+1])&set(BEAM[i[1]]+BEAM[i[1]+1]))
		#commbeam1.append(list(set(BEAM[i[0]]+BEAM[i[0]+1])&set(BEAM[i[1]]+BEAM[i[1]+1])))

	#commbeam1 = sum(commbeam1,[])
	BEAM = np.array(BEAM)
	flagchain = [np.array([]) for _ in range(14)]
	flag = [np.array([]) for _ in range(14)]

	#Compares near by frequency bins in all beams if a hit is found in one bin.  If near by bins are hit then reject it.  
	#NearFreq = 5*FS
	#NearFreq = 50000
	#First combine data from all beam except the beams of interest 
	for b in permutations(range(0,14,2),2):
		flagchain[b[0]] = list(chain(flagchain[b[0]],BEAM[b[1]],BEAM[b[1]+1]))
		flagchain[b[0]+1] = list(chain(flagchain[b[0]+1],BEAM[b[1]],BEAM[b[1]+1]))
		#print i[0],i[0]+1,i[1],i[1]+1

	#To check 	
	#for i,j in enumerate(flagchain):
		#print i,len(j),len(BEAM[i])
		#print i,np.array(j)-5

	#Then subtract beam in question's each hit chan from all the other beams hit chan
	for i in range(0,14):
		for j in BEAM[i]:
			minval =  min(abs(k-j) for k in flagchain[i])
			if(minval<=NearFreq):
				commbeam1.append(j)
			flag[i] = np.append(flag[i], minval)
	
	
	#Comment following for speedy test 
	# To remove fixed line from the ADC at coarse chan bins [0,512,1024,1536,1880,2048,2560,3072]
	try:
		ADCchan = ([1024,1536,2048,2560,3072][[1006,1326,1966,2286,2926].index(coarseID)] - coarseID)*FS
		if(ADCchan not in commbeam1):
			commbeam1.append(ADCchan)
	except:
		pass

	#To remove known RFI lines for a given coarse chan 
	FixChan = ([[],[],[1880,1885,1886],[],[],[],[],[]][[1006,1326,1646,1966,2286,2606,2926,3246].index(coarseID)])
	if(any(FixChan)):
		FixChan =  np.array(FixChan) - coarseID
		FixChan = [i*FS+j for i in FixChan for j in range(FS)]	
		FixChanFlag = []
		for i in range(0,14):
			if(any(list(set(BEAM[i])&set(FixChan)))):
				FixChanFlag = list(set(BEAM[i])&set(FixChan))
				#print FixChanFlag		
				commbeam1 = list(chain(commbeam1,FixChanFlag))
	

	return commbeam1

def indexing(DATABEAM2,DATABEAM,NS,FS,coarseID):
	BEAM = [[] for _ in range(14)]
	BEAM2 = [[] for _ in range(14)]

	#print DATABEAM
	for i in range(14):
		if(DATABEAM[i] is not None):
			for j in DATABEAM[i]:	# Extract chan number of compare from each beam 
				#print j
				BEAM[i].append(int(j[2]*FS + j[3]))

	for i in range(14):
		if(DATABEAM2[i] is not None):
			for j in DATABEAM2[i]:
				BEAM2[i].append(int(j[2]*FS + j[3]))

	commbeam = compBeam(BEAM,NS,FS,coarseID)
	#print commbeam

	indx = [[] for _ in range(14)]

	for i in range(0,14):
		for j in commbeam:	
			if j in BEAM2[i]:
				indx[i].append(BEAM2[i].index(j))
		indx[i] = list(set(indx[i]))
	
	#for i in range(0,14):
		#print len(indx[i])
		#print DATABEAM[i][0].header['NHITS'],len(indx[i])		


	return indx


if __name__ == "__main__":
	filename = sys.argv[1]
	print "Input: %s" % (filename)
	fileHandle = fits.open(filename)
	#ofileHandle = pyfits.open(ofilename,mode='ostream')
	ofileHandle = fileHandle
	DATABEAM = [[] for _ in range(14)]
	beamindx = [i for i in range(14)]
	delim = ['GBTSTATUS', 'AOSCRAM'] # Tables delimiting timesteps
	nSteps = 0
	ntable = 2 # First two tables are extra header information they are being ignored  
	NS = fileHandle[0].header['NSUBBAND']
	FS = fileHandle[0].header['NCHAN']
	coarseID = fileHandle[1].header['COARCHID']
	RA = []
	DEC = []

	flen = len(fileHandle)
	BLDATA = [[] for _ in range(14)]

	#Each elemenent is a time stamp, each time stamp has 14 array of each beam, each beam array has number of hits array, hits array has 4 entry
	#FULLDATA[TIMESTAMP][BEAM][HITnum][4 elements] 
	FULLDATA = []

	print "Starting first scan..."

	for ind, table in enumerate(fileHandle):
		if 'EXTNAME' in table.header:
			if 'ETHITS' == table.header['EXTNAME']:
				#print table.header['EXTNAME']
				if(table.header['BEAMPOL'] in beamindx):
					bp = int(table.header['BEAMPOL'])
					for j in table.data:
						DATABEAM[bp].append(j)
				if(ind >= flen-1):
					FULLDATA.append(DATABEAM)
			elif ((any(elem in table.header['EXTNAME'] for elem in delim) or (ind >= flen-1)) and ind>1):	
				#print len(DATABEAM[0])
				FULLDATA.append(DATABEAM)
				DATABEAM = [[] for _ in range(14)]
				ntable=ind
				nSteps+=1

		if (ind+1)/float(flen)*100%1 < .001:
                        update_progress((ind+1)/float(flen))

		
	#print FULLDATA[0][0:4]

	print "First scan over, Starting second scan..."

	FL = len(FULLDATA)

	nSteps = 0
	ntable = 2
	#BL = 5

	DATABEAM2 = [[] for _ in range(14)]

	for ind, table in enumerate(ofileHandle):
		if (ind+1)/float(flen)*100%1 < .001:
			update_progress((ind+1)/float(flen))

		if 'EXTNAME' in table.header:
			if ((any(elem in table.header['EXTNAME'] for elem in delim) or (ind >= flen-1)) and ind>1):
				try:
					DATABEAM = [[] for _ in range(14)]
					if(nSteps > BL and nSteps < FL-BL):
						for s in range(nSteps-BL,nSteps+BL+1):
							#print s,nSteps
							for b in range(0,14):
								DATABEAM[b].extend(FULLDATA[s][b])
								#print FULLDATA[s][b]
							#print nSteps,DATABEAM
							#DATABEAM = DATABEAM + FULLDATA[nSteps]
					else:
						#print nSteps
						DATABEAM = FULLDATA[nSteps]
					#print DATABEAM
					DATABEAM2 = FULLDATA[nSteps]
					indx = indexing(DATABEAM2,DATABEAM,NS,FS,coarseID)
					for i in range(0,15):
						try:
							if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
								bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
								ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
								ofileHandle[ntable+i].header['NHITS'] = len(ofileHandle[ntable+i].data)
								ofileHandle[ntable+i].header['NAXIS2'] = len(ofileHandle[ntable+i].data)
								#print x,ofileHandle[ntable+i].header['NHITS'],len(ofileHandle[ntable+i].data)
						except:
							pass
				except:
					pass
				nSteps+=1
				ntable=ind

	print "Second scan over. Writing FITS file..."

	ofileHandle.writeto("test.fits",clobber=True)
					
