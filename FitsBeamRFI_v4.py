"""
Analyzes and summarizes FITS file data

Use: Python analyzeFits.py <filename> <output directory>
"""
import os, sys, glob
import fitsio
import pyfits
import matplotlib as mplt
mplt.use('Agg') # Must be called before importing pyplot if X server is not allowed. Using this we can still be able to print in a file
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np
import time

# Constants and parameters
tStep = 1.43


def getFitsData(fileHandle):
	"""
	Extracts FITS data for other functions
	"""

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


	# Seperate data by tables
	ccpwrs = [] # Table seperated ccpwrs
	data = [] # Table seperated data

	# Seperate data by time steps
	nSteps = -1 # Number of time steps
	delim = ['GBTSTATUS', 'AOSCRAM'] # Tables delimiting timesteps
	rows = [] # Single data vector
	times = [] # Times for rows
	BEAMPOL = 0
	ABEAMPOL = [[] for _ in range(14)]
	DATABEAM = []
	RADEC = []

	# Each element is an ETHITS table
	# Columns: DETPOW, MEANPOW, COARCHAN, FINECHAN
	# Parse data
	for ind, table in enumerate(fileHandle):
		try: 
			if 'ETHITS' == table.header['EXTNAME']:
#				print table.header['BEAMPOL']
#				if(table.header['BEAMPOL']):
				BEAMPOL = int(table.header['BEAMPOL'])
#				if(BEAMPOL==0): print nSteps
				if(table.header['RA']):
				        RA = float(table.header['RA'])
				if(table.header['DEC']):
				        DEC = float(table.header['DEC'])
				DATABEAM.append(table.data)	
				for row in table.data:					
#					times.append(tStep * nSteps)
#					if row[2] == 0:
#						rows.append(row[3])
#					else:
#						rows.append(row[2]*2**19 + row[3])
					ABEAMPOL[BEAMPOL].append(float(row[0]/row[1]))							
#				if(tStep*nSteps<=3):
#					print tStep*nSteps,RA,DEC,BEAMPOL			
			elif 'CCPWRS' in table.header['EXTNAME']:
				ccpwrs.append(table.data)
			elif any(elem in table.header['EXTNAME'] for elem in delim):
				outblock = compBeam(DATABEAM) # Compare all beam for a given time stamp	
#				pyfits.writeto("test.fits",outblock,clobber='True')
				DATABEAM = [] # Clear all beam array 
				nSteps += 1
		except:
			pass


		# Progress bar
		if (ind+1)/float(len(fileHandle))*100%1 < .001: 
			update_progress((ind+1)/float(len(fileHandle)))


	return ABEAMPOL

	

def PLOTBEAMPOL(ABEAMPOL,outpdf,date):
	
	print "Ploting beams"
	datep = str("Date :"+ date[0:10])
	timep = str("Time :"+ date[11:19])

	mplt.rcParams['axes.linewidth'] = 1
	#plt.subplots_adjust(hspace = .5)
	#plt.subplots_adjust(wspace = .3)
	maxbin = int(max([max(j) for j in ABEAMPOL]))
	print maxbin
	totbin = 50
	histbins = np.linspace(20,maxbin,totbin)
	print histbins
	
	params = {'legend.fontsize': 10,'legend.linewidth': 2}
	plt.rcParams.update(params)

	#ax1 = plt.subplot2grid((6,3), (2,1), rowspan=2,colspan=1)
	ax1 = plt.subplot2grid((2,2), (1,1),rowspan=1,colspan=1)
	ax1.set_ylabel('Number of Hits',fontsize=10, fontweight='bold')
	ax1.set_xlabel('Relative Power (DETPOW/MEANPOW)',fontsize=10, fontweight='bold')
#	ax1.set_title("Beam 0 & 1")
#	plt.plot(ABEAMPOL[1],label="Beam 1", color='b',marker='.',linestyle='None')	
#	plt.plot(ABEAMPOL[0],label="Beam 0", color='r',marker='.',linestyle='None')
	plt.yscale('log',nonposy='clip')
	#plt.xlim('auto')
	#plt.ylim('auto')
	#plt.xlim(0,800)
	plt.xlim(20,40)
	#plt.ylim(1,100)

	y,x = np.histogram(ABEAMPOL[0],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[0],totbin)
	plt.step(x[:-1],y,color='b',label="Beam 1",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[2],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[2],totbin)
	plt.step(x[:-1],y,color='r',label="Beam 2",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[4],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[4],totbin)
	plt.step(x[:-1],y,color='c',label="Beam 3",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[6],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[6],totbin)
	plt.step(x[:-1],y,color='m',label="Beam 4",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[8],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[8],totbin)
	plt.step(x[:-1],y,color='k',label="Beam 5",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[10],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[10],totbin)
	plt.step(x[:-1],y,color='CadetBlue',label="Beam 6",alpha=0.6)

	y,x = np.histogram(ABEAMPOL[12],bins=histbins)
	#print max(x[:-1])
	#y,x = np.histogram(ABEAMPOL[12],totbin)
	plt.step(x[:-1],y,color='DarkGoldenRod',label="Beam 7",alpha=0.6)

	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)

	fig = plt.gcf()
	fig.set_size_inches(15,12)
	plt.legend(frameon=False)
	plt.savefig(outpdf, bbox_inches='tight')

	plt.show()


def main(filename,outpdf):
	"""
	Summarizes and analyzes FITS file
	"""

	# Get data
#	print('Parsing Data')
#	fileHandle = fitsio.FITS(f) # File isn't opened until fileHandle is used
	fileHandle = pyfits.open(f) # Open Fits file 
	date = fileHandle[0].header['DATE']	
	ABEAMPOL = getFitsData(fileHandle) # This opens fileHandle, will take a few minutes	

#	Compare All beams
	

	PLOTBEAMPOL(ABEAMPOL,outpdf,date)

	# End
	fileHandle.close
	print('Done')

if __name__ == "__main__":
	start = time.time()
	print('Starting timer')
	f = sys.argv[1]
	filename = f[f.rindex('/')+1:f.index('.fits')]
	filepath = f[0:f.rindex('/')]
	outfile = filename + "mBeamRFIClean.fits"
	outpdf = filename + ".beamSNR.pdf"
	outDir = sys.argv[2]
	if outDir[-1] != '/':
		sys.exit('Exiting: Output directory must end with /')
	elif os.path.isdir(outDir) == False:
		sys.exit('Exiting: Output directory does not exist')
	main(filename,outpdf)
	print('Total Time: %d Seconds' %(time.time()-start))
