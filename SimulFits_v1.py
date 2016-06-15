'''
This prog simulates fits file. A fits file should be given as an input. It uses that fits 
files leangth and header information to create a new simulated fits file 
'''

import pyfits
from astropy.io import fits
import os,sys
from random import randint
import numpy as np

def hitrandom(hitpoint):
#	pass
#	hitpoint = list(hitpoint)
#	print hitpoint
#	of[ind].data[j][0] = of[ind].data[j][1] + randint(0,1000)
	hitpoint[1] = randint(1,100)
#	hitpoint[1] = 1000
	hitpoint[0] = hitpoint[1] + randint(0,1000)
	hitpoint[2] = np.unint16(-32768)
#	hitpoint[2] = np.array(-32768,dtype='>i4')
#	hitpoint = tuple(hitpoint)
#	print hitpoint
#	hitpoint[2] = int(randint(0,314))
#	print temp,hitpoint[2]
#	hitpoint[3] = randint(0,2.6e8)
#	print hitpoint	
#	hitpoint1 = np.zeros_like(hitpoint)
#	hitpoint = np.zeros(1,

def timehitrandom(n,NS,FS,scaledata):
	detpow = np.random.uniform(20,100,n) 
	meanpow = np.random.uniform(1,20,n)
	
	bzero3 = scaledata[0]
	bscale3 = scaledata[1]
	bzero4 = scaledata[2]
	bscale4 = scaledata[3]

#	corchan = np.uint16(np.random.randint(1,NS,size=n))	
	corchan = np.int16(np.array([i-bzero3 for i in np.random.randint(1,NS,size=n)]))
#	corchan = np.zeros(n)
#	finchan = np.uint32(np.random.randint(1,FS,size=n))
#	finchan = np.zeros(n)
	finchan = np.int32(np.array([i-bzero4 for i in np.random.randint(1,FS,size=n)]))

	c1 = pyfits.Column(name='DETPOW',format='1E',array=detpow)
	c2 = pyfits.Column(name='MEANPOW',format='1E',array=meanpow)
#	c3 = pyfits.Column(name='COARCHAN',format='1I',array=corchan,bzero=bzero3,bscale=bscale3)
	c3 = pyfits.Column(name='COARCHAN',format='1I',array=corchan)
#	c4 = pyfits.Column(name='FINECHAN',format='1J',array=finchan,bzero=bzero4,bscale=bscale4)
	c4 = pyfits.Column(name='FINECHAN',format='1J',array=finchan)
	tbhdu = pyfits.new_table([c1, c2, c3, c4])
#	return	tbhdu.data
	return  tbhdu

if __name__ == "__main__":
#	filename = '/Users/vishalgajjar/SETI/serendip6_eth2_AO_327MHz_1006_20160124_172014.fits' # Small file
	filename = '/Users/vishalgajjar/SETI/data/serendip6_eth2_AO_327MHz_1646_20151128_171443.fits' # Big file 
	fileHandle = pyfits.open(filename)
	of = fileHandle
#	of = fileHandle[0:4]
#	Header information	
#	of[0] = fileHandle[0]
#	of[1] = fileHandle[1]
#	header = 


	delim = ['GBTSTATUS', 'AOSCRAM'] # Tables delimiting timesteps
	nSteps = -1
	beamindx = [i for i in range(14)]

	outfile = "simulTest.fits"

	hdr0  = pyfits.PrimaryHDU()
	hdr0.header = fileHandle[0].header
	hdr1 = pyfits.BinTableHDU()
	hdr1 = fileHandle[1]
	
	NS = fileHandle[0].header['NSUBBAND']
	FS = fileHandle[0].header['NCHAN']

	thdr = pyfits.HDUList([hdr0, hdr1])

#	Write necessary headers
	thdr.writeto(outfile,clobber=True)
	
#	fits.writeto(outfile,of[0].header,clobber=True)

#	for ind, table in enumerate(fileHandle[0:31]):
	for ind, table in enumerate(of):
		try:
			if 'ETHITS' == table.header['EXTNAME']:
				if(table.header['BEAMPOL'] in beamindx):
#					print ind,table.header['BEAMPOL'],table.header['NHITS']
#					data = table.data
#					hdr = fileHandle[5].header
					hdr = table.header 
#					print hdr
					try: 
						bzero3 = int(table.header['TZERO3']) 
						bscale3 = int(table.header['TSCAL3'])
						bzero4 = int(table.header['TZERO4'])
						bscale4 = int(table.header['TSCAL4'])
						scaledata = [bzero3,bscale3,bzero4,bscale4]
					except:
#						scaledata = [32768,1,2147483648,1]
						scaledata = [0,1,0,1]

					data1 =	timehitrandom(table.header['NHITS'],NS,FS,scaledata)
#					hdr = data1.header
					data2 = data1.data
#					for j in range(int(table.header['NHITS'])):							
#						data = table.data
#						hdr = table.hdr
#						pyfits.append(outfile,data,hdr)
#						pass
#					 	hitrandom(of[ind].data[j])
#						of[ind].data.scale
#						print temp,of[ind].data[j][2]
#						print ind,of[ind].data[j]
#					print nSteps,ofileHandle[nSteps].data
#					print ind,table.header['BEAMPOL'],data2[3][2]
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
				print table.header['EXTNAME']
				nSteps+=1
			else:
				print "Error2"
		except:
			print "Error3"
			pass

	
	of.verify()
	of.writeto(outfile,clobber=True)
#	of.flush(outfile)
#	of.close()
	of1 = pyfits.open(outfile,uint=True)
#	of1 = fits.open(outfile)

#	To test
	for ind,table in enumerate(of1):
		try:
			if 'ETHITS' == table.header['EXTNAME']:
				for j in table.data:
#					pass
					print ind,table.header['BEAMPOL'],j
#					print j[2],of[ind].data[j][2]
		except:
			pass
