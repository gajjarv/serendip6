import pyfits
from astropy.io import fits
import os,sys
import numpy as np
from itertools import *
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import operator

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
	#NearFreq = 3*FS
	NearFreq = 10000
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
	#To remove fixed line from the ADC at coarse chan bins [0,512,1024,1536,1880,2048,2560,3072]
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

def indexing(DATABEAM,NS,FS,coarseID):
	BEAM = [[] for _ in range(14)]
	for i in range(14):
		if(DATABEAM[i] is not None):
			#print i,DATABEAM[i][0]
			for j in DATABEAM[i][0]:	# Extract chan number of compare from each beam 
				#print i,j
				BEAM[i].append(int(j[2]*FS + j[3]))

	commbeam = compBeam(BEAM,NS,FS,coarseID)
	#print commbeam

	indx = [[] for _ in range(14)]
	for i in range(0,14):
		for j in commbeam:	
			if j in BEAM[i]:
				indx[i].append(BEAM[i].index(j))
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
	data = []

	#Block length
	BL = 3

	flen = len(fileHandle)
	
	
	for ind, table in enumerate(fileHandle):
		if 'EXTNAME' in table.header:
			#print ind,table.header['EXTNAME']
			if 'ETHITS' == table.header['EXTNAME']:
				if(table.header['BEAMPOL'] in beamindx):
					bp = int(table.header['BEAMPOL'])
					if(bp==0):
						RA.append(table.header['RA'])
						DEC.append(table.header['DEC'])
					
					BLdata = [[] for _ in range(14)]
					'''
					if(ind>BL*15 and ind<flen-BL*15):
						#print ind,table.header['EXTNAME']  
						for bltable in fileHandle[ind-BL*14:ind+BL*14]:
							if 'EXTNAME' in bltable.header:
								if 'ETHITS' == bltable.header['EXTNAME']:
									if(bltable.header['BEAMPOL'] in beamindx):
										bp = int(bltable.header['BEAMPOL'])
										BLdata[bp].append(bltable.data)
						#print BLdata[bp]
						#DATABEAM[bp].append(BLdata[bp])		
					else:
					'''
					DATABEAM[bp].append(table.data)	
					#For the last ETHITS time stamp table  
					if(ind >= flen-1):
						#print ind,ntable,DATABEAM
						indx = indexing(DATABEAM,NS,FS,coarseID)
						for i in range(0,15):
							try:
								if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
									bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
									ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
							except:
								pass									
			elif any(elem in table.header['EXTNAME'] for elem in delim):
				if(ind>1):
					indx = indexing(DATABEAM,NS,FS,coarseID)
					for i in range(0,15):
						#print ind,ntable,i,table.header['EXTNAME']
						try:
							if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
								#print ind,ntable,table.header['EXTNAME'] 
								bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
								#print ntable,indx[bp]
								ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
						except:
							pass
					ntable=ind	
					DATABEAM = [[] for _ in range(14)]
					nSteps+=1
			'''
			elif(ind >= flen-1):
				indx = indexing(DATABEAM,NS,FS,coarseID)
				for i in range(0,15):
					try:
						if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
							bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
							ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
					except:
						pass
			'''
		
	
	#RADEC(RA,DEC)		
	ofileHandle.writeto("test.fits",clobber=True)
	
