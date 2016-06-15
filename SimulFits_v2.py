'''
This prog simulates fits file. A fits file should be given as an input. It uses that fits 
files leangth and header information to create a new simulated fits file 
'''

import pyfits
from astropy.io import fits
import os,sys
from random import randint
import numpy as np

def hitrandom(n,NS,FS,scaledata):
	
	#Randomly creates 1 or 2 hits per time step
	n = np.random.randint(1,3)
	detpow = np.random.uniform(20,100,n) 
	meanpow = np.random.uniform(1,20,n)
	
	bzero3 = scaledata[0]
	bscale3 = scaledata[1]
	bzero4 = scaledata[2]
	bscale4 = scaledata[3]

	corchan = np.int16(np.array([i-bzero3 for i in np.random.randint(1,NS,size=n)]))
	finchan = np.int32(np.array([i-bzero4 for i in np.random.randint(1,FS,size=n)]))

#	ET Signal
	detpow = np.append(detpow,1200.0)
	detpow = np.append(detpow,1300.0)
	meanpow = np.append(meanpow,10.0)
	meanpow = np.append(meanpow,10.0)
	corchan = np.append(corchan,np.int16(125-bzero3))
	corchan = np.append(corchan,np.int16(225-bzero3))
	finchan= np.append(finchan,np.int32(0-bzero4))
	finchan= np.append(finchan,np.int32(0-bzero4))

	c1 = pyfits.Column(name='DETPOW',format='1E',array=detpow)
	c2 = pyfits.Column(name='MEANPOW',format='1E',array=meanpow)
#	c3 = pyfits.Column(name='COARCHAN',format='1I',array=corchan,bzero=bzero3,bscale=bscale3)
	c3 = pyfits.Column(name='COARCHAN',format='1I',array=corchan)
#	c4 = pyfits.Column(name='FINECHAN',format='1J',array=finchan,bzero=bzero4,bscale=bscale4)
	c4 = pyfits.Column(name='FINECHAN',format='1J',array=finchan)
	tbhdu = pyfits.new_table([c1, c2, c3, c4])
	return  tbhdu

if __name__ == "__main__":
#	filename = '/Users/vishalgajjar/SETI/serendip6_eth2_AO_327MHz_1006_20160124_172014.fits' # Small file
	filename = '/Users/vishalgajjar/SETI/data/serendip6_eth2_AO_327MHz_1646_20151128_171443.fits' # Big file 
	fileHandle = pyfits.open(filename)
	of = fileHandle

	delim = ['GBTSTATUS', 'AOSCRAM'] # Tables delimiting timesteps
	nSteps = -1
	beamindx = [i for i in range(14)]

	outfile = "simulTest.fits"

	NS = fileHandle[0].header['NSUBBAND']
	FS = fileHandle[0].header['NCHAN']


	for ind, table in enumerate(of):
		try:
			if 'ETHITS' == table.header['EXTNAME']:
				if(table.header['BEAMPOL'] in beamindx):
					hdr = table.header 
					try: 
						bzero3 = int(table.header['TZERO3']) 
						bscale3 = int(table.header['TSCAL3'])
						bzero4 = int(table.header['TZERO4'])
						bscale4 = int(table.header['TSCAL4'])
						scaledata = [bzero3,bscale3,bzero4,bscale4]
					except:
#						scaledata = [32768,1,2147483648,1]
						scaledata = [0,1,0,1]

					data1 =	hitrandom(table.header['NHITS'],NS,FS,scaledata)
					data2 = data1.data
					of[ind].header = hdr
					of[ind].data = data2
#					pyfits.append(outfile,data2,hdr)
#					pyfits.update(outfile,data2,hdr,ind)
				else:
					print "Error1"		
			elif any(elem in table.header['EXTNAME'] for elem in delim):
				hdr1 = pyfits.BinTableHDU()
				hdr1 = table
#				pyfits.append(outfile,hdr1.data,hdr1.header)
#				print table.header['EXTNAME']
				nSteps+=1
			else:
				print "Error2"
		except:
			print "Error3"
			pass

	
	of.verify()
	of.writeto(outfile,clobber=True)
	of.close()

#	To test
	test = 0
	of1 = pyfits.open(outfile,uint=True)
	if(test):
		for ind,table in enumerate(of1):
			try:
				if 'ETHITS' == table.header['EXTNAME']:
					for j in table.data:
#						pass
						print ind,table.header['BEAMPOL'],j
#						print j[2],of[ind].data[j][2]
			except:
				pass
