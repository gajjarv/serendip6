import pyfits
from astropy.io import fits
import os,sys
import numpy as np

def compBeam(BEAM):

	
	#Compare all beams and extract bin number which occurs in all beams
	commbeam = list(set(BEAM[0]+BEAM[1])&set(BEAM[2]+BEAM[3])&set(BEAM[4]+BEAM[5])&set(BEAM[6]+BEAM[7])&set(BEAM[8]+BEAM[9])&set(BEAM[10]+BEAM[11])&set(BEAM[12]+BEAM[13]))
	
	commbeam1 = []

	#Compare all beams and extract hits which occurs in more than one beam
	for i in range(0,14,2):
		#for j in range(1,15,2):
		for j in range(0,14,2):
			if(set(BEAM[i]+BEAM[i+1])!=set(BEAM[j]+BEAM[j+1])):
				#print i,i+1,j,j+1,list(set(BEAM[i]+BEAM[i+1])&set(BEAM[j]+BEAM[j+1]))
				commbeam1.append(list(set(BEAM[i]+BEAM[i+1])&set(BEAM[j]+BEAM[j+1])))		
				
	commbeam1 = sum(commbeam1,[])
	
	#print commbeam
	#Compare all beams and extract bin number which occurs just once

	return commbeam1

def indexing(DATABEAM):
	BEAM = [[] for _ in range(14)]
	for i in range(14):
		if(DATABEAM[i][0].data is not None):
			for j in DATABEAM[i][0].data:	# Extract chan number of compare from each beam 
				BEAM[i].append(int(j[2]*2**19 + j[3]))
			

	commbeam = compBeam(BEAM)
	#print indepbeam

	indx = [[] for _ in range(14)]
	for i in range(0,14):
		for j in commbeam:	
			if j in BEAM[i]:
				indx[i].append(BEAM[i].index(j))

	
	#for i in range(0,14):
		#print len(indx[i])
		#print DATABEAM[i][0].header['NHITS'],len(indx[i])		


	return indx

if __name__ == "__main__":
	#filename = '/Users/vishalgajjar/SETI/serendip6_eth2_AO_327MHz_1006_20160124_172014.fits' # Small file
	#filename = '/Users/vishalgajjar/SETI/data/serendip6_eth2_AO_327MHz_1646_20151128_171443.fits' # Big file 
	filename = '/Users/vishalgajjar/SETI/simulTest.fits'
	#ofilename = '/Users/vishalgajjar/SETI/test.fits'
	#fileHandle = pyfits.open(filename) # This technique was giving wrong output 
	fileHandle = fits.open(filename)
	#ofileHandle = pyfits.open(ofilename,mode='ostream')
	ofileHandle = fileHandle
	DATABEAM = [[] for _ in range(14)]
	beamindx = [i for i in range(14)]
	delim = ['GBTSTATUS', 'AOSCRAM'] # Tables delimiting timesteps
	nSteps = 0
	ntable = 2 # First two tables are extra header information they are being ignored  
	flen = len(fileHandle)

	for ind, table in enumerate(fileHandle):
		if 'EXTNAME' in table.header:
			if 'ETHITS' == table.header['EXTNAME']:
				#print table.header['EXTNAME']
				if(table.header['BEAMPOL'] in beamindx):
					bp = int(table.header['BEAMPOL'])
					DATABEAM[bp].append(table)

					#For the last ETHITS time stamp table  
					if(ind >= flen-1):
						#print ind,ntable,DATABEAM
						indx = indexing(DATABEAM)
						for i in range(0,15):
							try:
								if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
									bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
									ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
							except:
								pass
			elif any(elem in table.header['EXTNAME'] for elem in delim):
				if(ind>1):
					indx = indexing(DATABEAM)
					for i in range(0,15):
						#print ind,ntable,i,table.header['EXTNAME']
						try:
							if(ofileHandle[ntable+i].header['BEAMPOL'] in beamindx):
								#print ind,ntable,table.header['EXTNAME'] 
								bp = int(ofileHandle[ntable+i].header['BEAMPOL'])
								#print ntable,bp
								ofileHandle[ntable+i].data = np.delete(ofileHandle[ntable+i].data,indx[bp])
						except:
							pass
					ntable=ind	
					DATABEAM = [[] for _ in range(14)]
					nSteps+=1
		
		#else:
			#pass	
	ofileHandle.writeto("test.fits",clobber=True)
