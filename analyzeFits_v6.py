"""
Analyzes and summarizes FITS file data

Use: Python analyzeFits.py <filename> <output directory>
"""

import os
import fitsio
import sys
import matplotlib
matplotlib.use('Agg') # Must be called before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages # Save plots to PDF
import matplotlib.cm as cm
from pyPdf import PdfFileWriter, PdfFileReader# Merge PDF pages
import numpy as np
# For drawing paragraphs in PDF
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.lib.pagesizes import letter #, landscape
import time
import math
from scipy.stats import chi2
from scipy.stats import ncx2

# Constants and parameters
tStep = 1.43
lowthr = 20
highthr = 2200

def bin2if(CCbin,FCbin,TELID):
        instance = 1
#       instance = 2    
        cc = CCbin
        fc = np.uint32(FCbin)
	#print fc.dtype

	TELID = TELID.strip()
	if TELID == 'GBTSTATUS':
  		clock = 3000
        	cc_per_sys = 4096
        	cc_per_instance = 512
        	fc_per_cc = pow(2.0, 19);
	elif TELID == 'AOSCRAM':
		clock = 896
		cc_per_sys = 4096
		cc_per_instance = 320
		fc_per_cc = pow(2.0, 17);
	else:
		print "Unknown telescope !!"
		return birdibin

	band_width      = clock/2
        cc_bin_width    = band_width/cc_per_sys
        fc_bin_width    = band_width/(cc_per_sys * fc_per_cc)
        resolution = fc_bin_width * 1000000
        #sys_cc  = instance*cc_per_instance + cc
	sys_cc = cc	

	#print fc.dtype
	#print FCbin,CCbin,np.int32((0xFFFE0000 if(fc & 0x00010000) else 0) | (fc & 0x0001FFFF))

	#remember that lower half of a cc is the higher in freq!
	if TELID == 'GBTSTATUS':
		signed_fc = np.int32((0xFFF80000 if(fc & 0x00040000) else 0) | (fc & 0x0007FFFF))
	elif TELID == 'AOSCRAM':
		signed_fc = np.int32((0xFFFE0000 if(fc & 0x00010000) else 0) | (fc & 0x0001FFFF))
	else:
		print "Unknown telescope !!"
		return birdibin
	
	sys_fc = fc_per_cc * sys_cc + signed_fc	
	IFbin = ((sys_fc + fc_per_cc) * resolution)/(1000000)

	'''
	#old technique
	if(fc > fc_per_cc/2):
                corrected_fc = fc - fc_per_cc/2
        else:
                corrected_fc = fc + fc_per_cc/2

        corrected_fc -= fc_per_cc/2
	sys_fc = fc_per_cc * sys_cc + corrected_fc
	print sys_fc,corrected_fc,sys_cc
        IFbin = (sys_fc * resolution) / 1000000
	'''

	'''
        IfDiff = float(ifbirdi - IfCenter)
        SSB = SSB.strip()

        if(SSB == "lower"):
                Skyfreq = SkyCenter - IfDiff
        elif(SSB == "upper"):
                Skyfreq = SkyCenter + IfDiff
#       else:
#               print "Error in SSB"
#               Skyfreq = 0

        return(Skyfreq)
	'''

	return(IFbin)

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
	CC = []
	IF = []

	NS = fileHandle[0].read_header()['NSUBBAND']
	FS = fileHandle[0].read_header()['NCHAN']
	coarseID = fileHandle[1].read_header()['COARCHID']
	TELID = fileHandle[1].read_header()['EXTNAME']

	print TELID

	# Each element is an ETHITS table
	# Columns: DETPOW, MEANPOW, COARCHAN, FINECHAN
	# Parse data
	for ind, table in enumerate(fileHandle):
		try: 
			if 'ETHITS' in table.read_header()['EXTNAME']:
				data.append(table.read())
				BEAMPOL = int(table.read_header()['BEAMPOL'])
				for row in table:
					if(lowthr<float(row[0]/row[1])<highthr):
						times.append(tStep * nSteps)
						if row[2] == 0:
							rows.append(row[3])
						else:
							#rows.append(row[2]*2**19 + row[3])
							rows.append(row[2]*FS + row[3])
						#print row[2]*FS + row[3]
						CC.append(coarseID+row[2])	
						IFbin = bin2if(coarseID+row[2],row[3],TELID)
						print row[0],row[1],IFbin,tStep*nSteps,BEAMPOL
						IF.append(IFbin)
						#print IFbin
						#if(coarseID+row[2]==1885):
							#print coarseID+row[2],row[0]/row[1],row[3]
						#IF.append(bin2if
			elif 'CCPWRS' in table.read_header()['EXTNAME']:
				ccpwrs.append(table.read())
			elif any(elem in table.read_header()['EXTNAME'] for elem in delim):
				nSteps += 1
		except:
			pass

		# Progress bar
		#if (ind+1)/float(len(fileHandle))*100%1 < .001: 
			#update_progress((ind+1)/float(len(fileHandle)))

	
	#print max(rows),(NS)*2**19
	coarseID = fileHandle[1].read_header()['COARCHID']
	ccpwrsPoints = [[row[0], row[1], tStep*ind] for ind, table in enumerate(ccpwrs) for row in table]
	hitPoints = zip(rows, times)

	return data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID,NS,FS,CC,IF


def plotPDF(fileHandle, data, hitPoints, ccpwrsPoints, ccpwrs, coarseID,NS,FS,CC,IF):
	"""
	Prints plots to PDF
	"""
	
	def powerHist():
		"""
		Plots number of hits occured at each power level
		"""
	
		# Get DETPOW/MEANPOW
		relpow = [row[0]/row[1] for table in data for row in table]

		#TEST 
		#relpow = [row[0] for table in data for row in table]

		#relpow = relpow/np.mean(relpow)

		#Number of bins are equal to maximum SNR
		bins = int(max(relpow))
		#bins = 1e2

		#Number of expected random variable 
		N = NS*14*FS*(hitPoints[-1][1]/tStep)
		#print hitPoints[-1][1],hitPoints[-1][1]/tStep
		# Plot histogram
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		#y,x  = np.histogram(relpow,1000) # arbitrary bin amount
		#print "Normalization begin"
		#y,x,c = ax.hist(relpow,bins,color='white',normed=True)
		y,x,c = ax.hist(relpow,bins,color='white')
		#print "Normalized histogram"
		#print max(y)
		#print np.sum(y)
		#ax.plot(x[:-1],N*chi2.pdf(x[:-1],2),color='green',linewidth=1)
		#ax.plot(x[:-1],N*chi2.pdf(x[:-1],4),color='blue',linewidth=1)
		ax.plot(x[:-1],N*np.exp(-x[:-1]),color='red',linewidth=2)
		#ax.plot(x[:-1],np.sum(y)*ncx2.pdf(x[:-1],2,3),color='blue',linewidth=1)
		#ax.hist(relpow, bins,color='green')
		ax.axvline(lowthr, color='k',linestyle='--')
		ax.axvline(highthr, color='k',linestyle='--')
		plt.title('Power Histogram')
		plt.xlabel('Relative Power (DETPOW/MEANPOW)')
		#plt.xlabel('Norm detected Power (DETPOW)')
		plt.ylabel('Number of Hits')
		plt.autoscale(enable=True, axis='x', tight=True)
		plt.autoscale(enable=True, axis='y', tight=True)
		#ax1.set_yticks([0, 100, 1000, 10000])
		plt.yscale('log', nonposy='clip')
		#plt.xscale('log',nonposx='clip')
		#ax.set_yticks([1,10,100,1000])
		#pyplot.locator_params(axis='y',nbins=4)
		plt.ylim(ymin=1,ymax=max(y)+(max(y)-min(y))/10)
		#plt.ylim(ymin=1,ymax=10)
		plt.xlim(xmin=20,xmax=bins)
		#plt.plot(x,y,color='red',linewidth=3)
		#plt.show(block = False)

	def coarseHist():
		"""
		Plots the number of hits occured in each bin
		"""

		# Get COARSCHAN
		coarse = [coarseID + row[2] for table in data for row in table]

		# Plot histogram
		plt.figure()
		#plt.hist(coarse, bins = 1e3) # arbitrary bin amount
		#plt.hist(CC,bins=1e3)
		plt.hist(IF,bins=1e3)
		plt.title('Coarse Bin Histogram')
		#plt.xlabel('Coarse Bin Number')
		plt.xlabel('IF')
		plt.ylabel('Number of Hits')
		#plt.autoscale(enable=True, axis='x', tight=True)
		#plt.xlim(xmin=coarseID,xmax=coarseID+NS)
		plt.yscale('log', nonposy='clip')
		# plt.show(block = False)

	def coarseSpectrum():
		"""
		Plots power in each coarse channel (pole)
		"""

		# Get powers
		xpol = [row[0] for table in ccpwrs for row in table]
		ypol = [row[1] for table in ccpwrs for row in table]

		# Plot histogram
		plt.figure()
		plt.plot(xpol[0:512], '-r', label = 'XPOL')
		plt.plot(ypol[0:512], '-b', label = 'YPOL')
		plt.title('Coarse Spectrum')
		plt.legend(loc = 'center right')
		plt.xlabel('Coarse Bin Number')
		plt.autoscale(enable=True, axis='x', tight=True)
		# plt.xlim([0,511])
		plt.ylabel('Power')
		plt.yscale('log')
		# plt.show(block = False)

	def waterfallHits():
		"""
		Plots the time at which a signal was received vs its frequency
		"""

		plt.figure()
		#axes = plt.gca()
		plt.plot(*zip(*hitPoints), rasterized=True, linestyle='', color='black', marker='o', markersize=1)
		#print hitPoints
		plt.title('Waterfall Hits')
		plt.xlabel('Fine Channel Number (Starting at channel ID)')
		plt.ylabel('Time (s)')
		#plt.autoscale(enable=True, axis='x', tight=True)
		plt.xlim(xmin=0,xmax=FS*NS)
		#plt.xlim(xmin=0,xmax=NS*2**19)
		#axes.set_xlim([0,1e8])
		plt.ylim(ymin=0)

	def waterfallCoarse():
		"""
		Plots CCPWRS over time. Color indicates Power.
		"""

		plt.figure()
		plt.figure(figsize=(10,10))
		plt.subplot(2,1,1)
		xpol = [row[0] for table in ccpwrs for row in table]
		imgX = np.array(xpol)
		imgX = imgX.reshape(len(ccpwrs),512)
		plt.imshow(imgX.astype(int), origin='lower', aspect='auto', cmap = cm.hot)
		plt.title('X-Pole CCPWRS')
		plt.ylabel('No. Time Steps (Time/Time Step)')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.xlabel('Coarse Channel ID')


		plt.subplot(2,1,2)
		ypol = [row[1] for table in ccpwrs for row in table]
		imgY = np.array(ypol)
		imgY = imgY.reshape(len(ccpwrs),512)
		plt.imshow(imgY.astype(int), origin='lower', aspect='auto', cmap = cm.hot)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.title('Y-Pole CCPWRS')
		plt.ylabel('No. Time Steps (Time/Time Step)')
		plt.xlabel('Coarse Channel ID')
		plt.subplots_adjust(hspace=0.4)

	# Create PDF
	pdf = PdfPages('%s%s_plots.pdf' %(outDir, filename))
	
	# Generate plots and save to PDF
	powerHist()
	pdf.savefig()
	plt.close()
	coarseHist()
	pdf.savefig()
	plt.close()
	waterfallHits()
	pdf.savefig()
	plt.close()
	# coarseSpectrum()
	if ccpwrs:
		waterfallCoarse()
		pdf.savefig()
		plt.close()
	# Close PDF and figures
	pdf.close()
	plt.close('all') # Just to be sure

def metaSummary(fileHandle, data, nSteps, hitPoints):
	"""
	Prints meta data to cover page as a summary of the FITS file
	"""

	def drawPage(meta):
		"""
		Creates cover page
		"""

		def coverPage(canvas, doc):
			"""
			Cover page format
			"""
			canvas.saveState()
			canvas.setFont('Times-Bold',16)
			canvas.drawCentredString(PAGE_WIDTH/2.0, PAGE_HEIGHT-108, Title)
			canvas.setFont('Times-Roman',9)
			canvas.restoreState()

		# PDF Parameters
		PAGE_HEIGHT=defaultPageSize[1]; PAGE_WIDTH=defaultPageSize[0]
		styles = getSampleStyleSheet()
		Title = 'FITS Summary'

		# Create cover page
		doc = SimpleDocTemplate('%s%s_meta.pdf' %(outDir, filename))
		content = [Spacer(1,2*inch)]
		style = styles["Normal"]
		for key in sorted(meta.keys()):
			text = ("%s: %s \n" % (key, meta[key]))
			p = Paragraph(text, style)
			content.append(p)
		doc.build(content, onFirstPage = coverPage)

	def getMeta():
		"""
		Gather select meta data
		"""

		def getnHits():
			"""
			Add up NHITS key value in all tables
			"""

			nHits = len(hitPoints)

			return nHits

		def getDuration():
			"""
			Duration of file
			"""

			duration = str(round((tStep*nSteps)/60, 2)) + ' Minutes'

			return duration

		def getTime():
			"""
			Get time of first table
			"""

			time = fileHandle[0].read_header()['DATE']

			return time

		def getAvgHits():
			"""
			Get mean and median hit counts for entire file
			For GBT: 16 means and medians (16 subands)
			For AO: 14 means and medians (14 Beampols)
			"""

			# if 'GBTSTATUS' in fileHandle[1].read_header()['EXTNAME']:
			# 	for i = range(len()

		def getFileInfo():
			"""
			Gets fileinfo of file on disk
			e.g. filesize on disk
			"""

			statInfo = os.stat(f)
			fileSize = round(statInfo.st_size/(1024.0**2), 2) # Bytes to MB
			fileSize = str(fileSize) + ' MB'

			return fileSize


		meta = {
				'FILENAME': filename,
				'NHITS': getnHits(),
				'TIME': getTime(),
				'DURATION': getDuration(),
				'FILE SIZE': getFileInfo() 
				}

		return meta


	drawPage(getMeta())

def pdfMerge():
	"""
	Merges generated PDFs into one called <filename>_summary
	Deletes the individual consituent PDFs
	"""

	def append_pdf(input, output):
		"""
		Combines PDF pages to be  merged
		"""

		[output.addPage(input.getPage(page_num)) for page_num in range(input.numPages)]

	# Merge PDFs
	output = PdfFileWriter()
	print outDir
	append_pdf(PdfFileReader(file('%s%s_meta.pdf' %(outDir, filename), 'rb')), output)
	append_pdf(PdfFileReader(file('%s%s_plots.pdf' %(outDir, filename), 'rb')), output)

	outputFile = file('%s%s_summary.pdf' %(outDir, filename), 'wb')
	output.write(outputFile)
	outputFile.close()

	# Delete PDFs
	os.remove('%s%s_plots.pdf' %(outDir, filename))
	os.remove('%s%s_meta.pdf' %(outDir, filename))

def main(filename):
	"""
	Summarizes and analyzes FITS file
	"""

	# Get data
	print('Parsing Data')
	fileHandle = fitsio.FITS(f) # File isn't opened until fileHandle is used
	data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID,NS,FS,CC,IF = getFitsData(fileHandle) # This opens fileHandle, will take a few minutes	
	# Print select meta data to PDF
	print('Printing Metadata to PDF')
	metaSummary(fileHandle, data, nSteps, hitPoints)

	# Plot data to PDF
	print('Plots printing, this may take a while')
	plotPDF(fileHandle, data, hitPoints, ccpwrsPoints, ccpwrs, coarseID,NS,FS,CC,IF)

	# Merge preceeding PDFs
	print('Cleaning up...')
	pdfMerge()

	# End
	fileHandle.close
	print('Done')

if __name__ == "__main__":
	start = time.time()
	print('Starting timer')
	f = sys.argv[1]
	filename = f[f.rindex('/')+1:f.index('.fits')]
	filepath = f[0:f.rindex('/')]
	outDir = sys.argv[2]
	if outDir[-1] != '/':
		sys.exit('Exiting: Output directory must end with /')
	elif os.path.isdir(outDir) == False:
		sys.exit('Exiting: Output directory does not exist')
	main(filename)
	print('Total Time: %d Seconds' %(time.time()-start))
