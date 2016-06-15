"""
Analyzes and summarizes FITS file data

Use: Python analyzeFits.py <filename> <output directory>
"""
import os, sys, glob
import fitsio
import matplotlib as mplt
mplt.use('Agg') # Must be called before importing pyplot if X server is not allowed. Using this we can still be able to print in a file
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np
import time
#from tabulate import tabulate
#import matplotlib as mplt
#import matplotlib.pyplot as plt


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
	RADEC = []

	# Each element is an ETHITS table
	# Columns: DETPOW, MEANPOW, COARCHAN, FINECHAN
	# Parse data
	for ind, table in enumerate(fileHandle):
		try: 
			if 'ETHITS' in table.read_header()['EXTNAME']:
				data.append(table.read())
				if(table.read_header()['BEAMPOL']):
					BEAMPOL = int(table.read_header()['BEAMPOL'])
#				if(BEAMPOL==0): print nSteps
				if(table.read_header()['RA']):
				        RA = float(table.read_header()['RA'])
				if(table.read_header()['DEC']):
				        DEC = float(table.read_header()['DEC'])
				for row in table:
					times.append(tStep * nSteps)
#					if row[2] == 0:
#						rows.append(row[3])
#					else:
#						rows.append(row[2]*2**19 + row[3])
					ABEAMPOL[BEAMPOL].append(float(row[0]/row[1]))		
					
#				if(tStep*nSteps<=3):
#					print tStep*nSteps,RA,DEC,BEAMPOL			
			elif 'CCPWRS' in table.read_header()['EXTNAME']:
				ccpwrs.append(table.read())
			elif any(elem in table.read_header()['EXTNAME'] for elem in delim):
				nSteps += 1
		except:
			pass


		# Progress bar
		if (ind+1)/float(len(fileHandle))*100%1 < .001: 
			update_progress((ind+1)/float(len(fileHandle)))

	coarseID = fileHandle[1].read_header()['COARCHID']
	ccpwrsPoints = [[row[0], row[1], tStep*ind] for ind, table in enumerate(ccpwrs) for row in table]
	hitPoints = zip(rows, times)

	return data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID, ABEAMPOL

def PLOTBEAMPOL(ABEAMPOL,outpdf,date):
	
	print "Ploting beams"
	datep = str("Date :"+ date[0:10])
	timep = str("Time :"+ date[11:19])

	mplt.rcParams['axes.linewidth'] = 2
	plt.subplots_adjust(hspace = .5)
	plt.subplots_adjust(wspace = .3)

	ax1 = plt.subplot2grid((6,3), (2,1), rowspan=2,colspan=1)
	ax1.set_xlabel('# Hits',fontsize=14, fontweight='bold')
	ax1.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
#	ax1.set_title("Beam 0 & 1")
	plt.plot(ABEAMPOL[1],label="Beam 1", color='b',marker='.',linestyle='None')	
	plt.plot(ABEAMPOL[0],label="Beam 0", color='r',marker='.',linestyle='None')

	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)

	ax2 = plt.subplot2grid((6,3), (4,1), rowspan=2,colspan=1)
	ax2.set_xlabel('# Hits',fontsize=14, fontweight='bold')
	ax2.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
#	ax2.set_title("Beam 2 & 3")
	plt.plot(ABEAMPOL[2],color='b',marker='.',linestyle='None',label="Beam 2")
        plt.plot(ABEAMPOL[3],color='r',marker='.',linestyle='None',label="Beam 3")
	handles, labels = ax2.get_legend_handles_labels()
	ax2.legend(handles, labels)
	
	ax3 = plt.subplot2grid((6,3), (3,2), rowspan=2,colspan=1)
	ax3.set_xlabel('# Hits',fontsize=14, fontweight='bold')
	ax3.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
#	ax3.set_title("Beam 4 & 5")
	plt.plot(ABEAMPOL[4],color='b',marker='.',linestyle='None',label="Beam 4")
        plt.plot(ABEAMPOL[5],color='r',marker='.',linestyle='None',label="Beam 5")
	handles, labels = ax3.get_legend_handles_labels()
	ax3.legend(handles, labels)

	ax4 = plt.subplot2grid((6,3), (1,2), rowspan=2,colspan=1)
	ax4.set_xlabel('# Hits',fontsize=14, fontweight='bold')
        ax4.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
#	ax1.set_title("Beam 0 & 1")
	plt.plot(ABEAMPOL[6],color='b',marker='.',linestyle='None',label="Beam 6")
        plt.plot(ABEAMPOL[7],color='r',marker='.',linestyle='None',label="Beam 7")
	handles, labels = ax4.get_legend_handles_labels()
	ax4.legend(handles, labels)

	ax5 = plt.subplot2grid((6,3), (0,1), rowspan=2,colspan=1)
	ax5.set_xlabel('# Hits',fontsize=14, fontweight='bold')
        ax5.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
	plt.plot(ABEAMPOL[8],color='b',marker='.',linestyle='None',label="Beam 8")
        plt.plot(ABEAMPOL[9],color='r',marker='.',linestyle='None',label="Beam 9")
	plt.text(1.5, 0.80, datep, horizontalalignment='center',fontsize=20,transform = ax5.transAxes)
	plt.text(1.5, 0.65, timep, horizontalalignment='center',fontsize=20,transform = ax5.transAxes)
#	ax5.text(0,0,"Date and Time")
	handles, labels = ax5.get_legend_handles_labels()
	ax5.legend(handles, labels)

	ax6 = plt.subplot2grid((6,3), (1,0), rowspan=2,colspan=1)
	ax6.set_xlabel('# Hits',fontsize=14, fontweight='bold')
        ax6.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
	plt.plot(ABEAMPOL[10],color='b',marker='.',linestyle='None',label="Beam 10")
        plt.plot(ABEAMPOL[11],color='r',marker='.',linestyle='None',label="Beam 11")
	handles, labels = ax6.get_legend_handles_labels()
	ax6.legend(handles, labels)

	ax7 = plt.subplot2grid((6,3), (3,0), rowspan=2,colspan=1)
	ax7.set_xlabel('# Hits',fontsize=14, fontweight='bold')
        ax7.set_ylabel('DETPOW/MEANPOW',fontsize=14, fontweight='bold')
	plt.plot(ABEAMPOL[12],color='b',marker='.',linestyle='None',label="Beam 12")
	plt.plot(ABEAMPOL[13],color='r',marker='.',linestyle='None',label="Beam 13")
	handles, labels = ax7.get_legend_handles_labels()
	ax7.legend(handles, labels)

	fig = plt.gcf()
	fig.set_size_inches(15,12)
	plt.savefig(outpdf, bbox_inches='tight')

	plt.show()


def main(filename,outpdf):
	"""
	Summarizes and analyzes FITS file
	"""

	# Get data
#	print('Parsing Data')
	fileHandle = fitsio.FITS(f) # File isn't opened until fileHandle is used
	date = fileHandle[0].read_header()['DATE']	
	data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID, ABEAMPOL = getFitsData(fileHandle) # This opens fileHandle, will take a few minutes	

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
