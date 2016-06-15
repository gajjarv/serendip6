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
	medstep = [] # Median for each time stamp
	temp = []
	flag = [] # Flag timestamps 	
	hitpow = [] # Hits power
	# Each element is an ETHITS table
	# Columns: DETPOW, MEANPOW, COARCHAN, FINECHAN
	# Parse data
	for ind, table in enumerate(fileHandle):
		try: 
			if 'ETHITS' in table.read_header()['EXTNAME']:
				data.append(table.read())	
				for row in table:
					times.append(tStep * nSteps)
					if(row[0]!=0):
						if row[2] == 0:
							rows.append(row[3])
						else:
							rows.append(row[2]*2**19 + row[3])
#						temp.append(row[0])
						hitpow.append(row[0])
			elif 'CCPWRS' in table.read_header()['EXTNAME']:
				ccpwrs.append(table.read())
			elif any(elem in table.read_header()['EXTNAME'] for elem in delim):
				nSteps += 1
				if temp: 
					medstep.append(np.mean(temp))
				del temp[:]
		except:
			pass
#		print table.read_header()	
		# Progress bar
		if (ind+1)/float(len(fileHandle))*100%1 < .001: 
			update_progress((ind+1)/float(len(fileHandle)))

#	for i in medstep: print i		
	print np.median(medstep),mad(medstep)
	print np.mean(medstep),np.std(medstep)
	for i,j in enumerate(medstep):
		if(medstep[i] >= (np.median(medstep) + mad(medstep))):
			flag.append(i)
			print medstep[i] 
	flag=sorted(flag,reverse=True)	
	
#	for i in flag:
#		del data[i]
#		del times[i]	

	coarseID = fileHandle[1].read_header()['COARCHID']
	ccpwrsPoints = [[row[0], row[1], tStep*ind] for ind, table in enumerate(ccpwrs) for row in table]
	hitPoints = zip(times,rows,hitpow)

	return data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def plotPDF(fileHandle, data, hitPoints, ccpwrsPoints, ccpwrs, coarseID):
	"""
	Prints plots to PDF
	"""
	
	def powerHist():
		"""
		Plots number of hits occured at each power level
		"""

		# Get DETPOW/MEANPOW
		relpow = [row[0]/row[1] for table in data for row in table]

		# Plot histogram
		plt.figure()
		plt.hist(relpow, bins = 1e3) # arbitrary bin amount
		plt.title('Power Histogram')
		plt.xlabel('Relative Power (DETPOW/MEANPOW)')
		plt.ylabel('Number of Hits')
		plt.autoscale(enable=True, axis='x', tight=True)
		plt.yscale('log', nonposy='clip')
		# plt.show(block = False)

	def coarseHist():
		"""
		Plots the number of hits occured in each bin
		"""

		# Get COARSCHAN
		coarse = [coarseID + row[2] for table in data for row in table]

		# Plot histogram
		plt.figure()
		plt.hist(coarse, bins = 1e3) # arbitrary bin amount
		plt.title('Coarse Bin Histogram')
		plt.xlabel('Coarse Bin Number')
		plt.ylabel('Number of Hits')
		plt.autoscale(enable=True, axis='x', tight=True)
		# plt.xlim(xmin=0)
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

	def waterfallcolor():
		"""
		It is waterfall in color intensity 
		"""
		plt.figure()
		plt.xlabel('Fine Channel Number (Starting at channel ID)')
                plt.ylabel('Time (s)')
                plt.autoscale(enable=True, axis='x', tight=True)
              # plt.xlim(xmin=0, xmax=max(x))
                plt.ylim(ymin=0)
			
	

	def waterfallHits():
		"""
		Plots the time at which a signal was received vs its frequency
		"""

		plt.figure()
		plt.plot(*zip(*hitPoints), rasterized=True, linestyle='', color='black', marker='o', markersize=1)
		plt.title('Waterfall Hits')
		plt.xlabel('Fine Channel Number (Starting at channel ID)')
		plt.ylabel('Time (s)')
		plt.autoscale(enable=True, axis='x', tight=True)
		# plt.xlim(xmin=0, xmax=max(x))
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
#	waterfallcolor()
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
	data, ccpwrs, hitPoints, ccpwrsPoints, nSteps, coarseID = getFitsData(fileHandle) # This opens fileHandle, will take a few minutes	
	# Print select meta data to PDF
	print('Printing Metadata to PDF')
	metaSummary(fileHandle, data, nSteps, hitPoints)

	# Plot data to PDF
	print('Plots printing, this may take a while')
	plotPDF(fileHandle, data, hitPoints, ccpwrsPoints, ccpwrs, coarseID)

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
