"""
Module contains spectrum processing routines used on all files, consists of using Lieber method to remove background, smoothing with Savitzky-Golay and computing normalized and non-normalized results.

Assumes tab-delimited files where first column is spectrum units and second column is intensity values.
"""



""""
Begin function definitions
"""



#--------------------------------------------------------------------------------------
def onepartfit(Data, Order=0):

	"""
	Single, contiguous modified polynomial fit
	"""

	from ModPolyFit import lieberfit


# Perform modified polynomial fit over entire data range
	Curvefit = lieberfit(Data, Order)
	
	SmoothCurve = np.copy(Curvefit)		


# Subtract background and shift above zero
	np.copyto(SmoothCurve, Data, where = SmoothCurve > Data)

	SubtractedResult=Data-SmoothCurve

	while any(SubtractedResult<=0):
		SubtractedResult+=10

	return(SubtractedResult, SmoothCurve)



#--------------------------------------------------------------------------------------
def threepartfit(Data, Order=5):

	"""
	Piecewise modified polynomial fit
	"""

	import statsmodels.api as sm
	from ModPolyFit import lieberfit
	

	IdxLeft = Divide1
	IdxRight = Divide2


# Perform piecewise modified polynomial fit over middle 4/5 and 1/3 extreme ranges
	Curvefit = np.zeros(len(Data))
# 	Middle = np.zeros(len(Data))
	Curvefit[:IdxLeft] = lieberfit(Data[:IdxLeft],Order-1)
	Curvefit[IdxRight:] = lieberfit(Data[IdxRight:],Order-1)
#	Middle[IdxOuter:-IdxOuter] = lieberfit(Data[IdxOuter:-IdxOuter],Order)
	Curvefit[IdxLeft:IdxRight] = lieberfit(Data[IdxLeft:IdxRight],Order)
			
			
# Heavy smoothing of background curve to remove discontinuities at dividers
	SmoothCurve=np.copy(Curvefit)
	StartSmooth = max(int(IdxLeft-10),0)
	TruncLeft = max(StartSmooth,IdxLeft)
# 	Use locally-weighted linear regression to smooth curves over appropriate region
	lowess = sm.nonparametric.lowess
	x = range(TruncLeft,Pixels)
	filt = lowess(Curvefit[TruncLeft:],x,frac=0.1)
	SmoothCurve[TruncLeft:] = filt[:,1]


# 	Subtract background and shift above zero
	np.copyto(SmoothCurve, Data, where = SmoothCurve > Data)
	
	SubtractedResult=Data-SmoothCurve
	
	while any(SubtractedResult<=0):
		SubtractedResult+=10


# 	SmoothCurve[StartSmooth:]=signal.savgol_filter(SmoothCurve[StartSmooth:],11,1)	
		
		
	return(SubtractedResult, SmoothCurve)



#--------------------------------------------------------------------------------------
def specproc(FileNameList, FileBase):

	"""
	Heavy lifting for spectrum processing
	"""


# Remove files (typically SkipIdx=1, the first file) due to instrument issues
	if SkipIdx > 0:
		FileIdx=[]
		Suffix='_'+str(SkipIdx)+FileExt
		for Idx,Name in enumerate(FileNameList):
			if not Suffix in Name:
				FileIdx.append(Idx)
		SkippedList=[FileNameList[i] for i in FileIdx]
	elif SkipIdx == 0:
		SkippedList = FileNameList
	else:
		raise Exception('Skipped file must be indexed with a positive integer.')


	numfiles=len(SkippedList)
	currfiles = [s.replace(ReadFolder,'').replace(FileExt,'') for s in SkippedList]
		
		
# Initiate matrices for spectrum and specunits data.
	FullName = ['' for i in range(numfiles)]
	Raw = np.zeros((Pixels,2,numfiles))   # 3d array
	SmoothData = np.zeros((Pixels,numfiles))   # 2d arrays here on
	Curve = np.zeros((Pixels,numfiles))
	SmoothCurve = np.zeros((Pixels,numfiles))
	SubtractedResult = np.zeros((Pixels,numfiles))


# Background file processing
	BGsub = np.zeros(Pixels)
	if args.bg:
		BGfileCount = 3   # must be set manually
		BGspectrum = np.zeros((Pixels,BGfileCount))
		for i in range(BGfileCount):
			BGdata = np.loadtxt(ReadFolder+args.bg+'_'+str(i+1+SkipIdx)+FileExt)
			BGspectrum[:,i] = BGdata[:,1]
		BGsub = np.mean(BGspectrum, axis=1)
			
			
# Load files and perform background fit
	for i in range(numfiles):
		Raw[:,:,i]=np.loadtxt(SkippedList[i],delimiter="\t")
		
	
# Check for instrument error where consecutive spectra share many identical values
		if i > 0: 
			DupChk = sum(Raw[:,1,i] == Raw[:,1,i-1])   # integer-valued, tol=1
			if DupChk > Pixels/2:
				print('')
				import warnings
				warnmsg='Possible duplicate spectra, {:} identical pixels'.format(DupChk)
				warnings.warn(warnmsg, UserWarning)
				print('Files:')
				print(SkippedList[i-1])
				print(SkippedList[i])
		
		
# Apply Savitzky-Golay filtering (11-point window, 1st order polynomial)
# Note: larger frames and lower order polynomial fits have stronger smoothing
		SmoothData[:,i]=signal.savgol_filter(Raw[:,1,i] - 0.8*BGsub,11,1) 


# Perform modified polynomial fit
		if str.lower(args.src)=='tweez':
			SubtractedResult[:,i],SmoothCurve[:,i] = threepartfit(SmoothData[:,i],fitordr)
		else:   # args.src is libs or opto
			SubtractedResult[:,i],SmoothCurve[:,i] = onepartfit(SmoothData[:,i],fitordr)

	MeanSmooth=np.mean(SmoothData,axis=1)
	MeanCurve=np.mean(SmoothCurve,axis=1)
	
	
# 	Perform peak fitting, mean, std, and output all to tab-delimited files

	def outputRslt(SaveExt):
		
		import pandas as pd
		
		if SaveExt == '-SF':   # No normalization
			RsltData = np.copy(SubtractedResult)
		elif SaveExt == '-SFN':   # p-norm normalization
			RsltData = normalize(SubtractedResult,norm='l1',axis=0)
		elif SaveExt == '-SFM':   # Min-max normalization
			MMscale = MinMaxScaler()
			RsltData = MMscale.fit_transform(SubtractedResult)
		else:
			raise Exception('Unknown processing method specified by programmer')


# 	Calculate mean & std dev of smoothed/fitted data and background removal curves for visualization and create dataframe
		Mean = np.mean(RsltData,axis=1)
		Std = np.std(RsltData,axis=1)
		specmin = np.amin(RsltData,axis=1)
		specmax = np.amax(RsltData,axis=1)
		
		dsstat = list(zip(specunits, Mean, Std, specmin, specmax))
		colstat = [xlbl, 'Mean', 'Std dev', 'Min', 'Max']
		dfstat = pd.DataFrame(data=dsstat, columns=colstat)


#	Find peaks of mean data with built in function that fits wavelets to data
		PeakWidths=np.arange(1,30)
		PkIdx = signal.find_peaks_cwt(Mean,PeakWidths)


# 	Compute peak locations and intensities for non-normalized and normalized results and create dataframe
		pkloc = specunits[PkIdx]
		pkint = Mean[PkIdx]
		
		dspks = list(zip(pkloc, pkint))
		colpks = ['Peak Locations', 'Peak Intensity']
		dfpks = pd.DataFrame(data=dspks, columns=colpks)
		dfpks.sort_values(by='Peak Intensity', ascending=False, inplace=True)
		dfpks.index = range(dfpks.shape[0])


# 		Use pandas to save results as an Excel file, include each spectrum
		dfrslt = pd.DataFrame(data=RsltData, columns=currfiles)

		dfwrite = pd.concat([dfstat, dfpks, dfrslt], axis=1)
		namewrite = WriteFolder+FileBase+SaveExt+'.xlsx'
		dfwrite.to_excel(namewrite, header=True, index=False)
	
	
# 		Visualize results just once not three times
		if SaveExt == '-SF':
			plt.rc('font', family = 'Arial', size='14')
			plt.ion()
	
			plt.figure()
			plt.plot(specunits,Mean,color='blue')
			plt.plot(pkloc,pkint+max(Mean)*0.03,"o",color="green",markersize=6)
			plt.title(FileBase+' - Subtracted with peaks')

			plt.figure()
			plt.plot(specunits,MeanCurve,color='red')
			plt.plot(specunits,MeanSmooth,color='blue')
			plt.title(FileBase+' - Smoothed & curve fit')
			plt.axvline(x=specunits[Divide1],color='k',ls='--',lw=1)
			plt.axvline(x=specunits[Divide2],color='k',ls='--',lw=1)

# 			Create min/max plot, first row is not shown in plot (fn expects header row)
			from plotting.makeplot import fillbtwn
# 			plt.ioff()
# 			fillbtwn(Wdata[:,:4], savename=FileBase)
# 			fillbtwn(Wdata[:,:4], savename=None)
# 			plt.ion()
	
	
# 	No need to print out results from list comprehension	
	NoDisplay = [outputRslt(s) for s in ['-SFM','-SFN','-SF']]
	
	
	return(MeanCurve)



"""
End function definitions
"""



#--------------------------------------------------------------------------------------
# Main

# As of 2015-09-29, this file is not configured for use as an importable module
if __name__ == '__main__':


# 	Define a few global variables that may need infrequent modification
	WriteFolder = '/Volumes/TRANSFER/Analysis/'
	CalibPath = '/Users/dbricare/Desktop/CBST Lab/Labwork/'


# 	Import libraries
	import re, argparse, time
	import numpy as np
	from scipy import signal
	import matplotlib.pyplot as plt
	from sklearn.preprocessing import MinMaxScaler, normalize

	import sys
	sys.path.append('/Users/dbricare/Documents/Python/mypyutil')
	from FileOpen import openfdiag
# 	from . import mypyutil


# 	Parse the arguments passed by the user
	parser = argparse.ArgumentParser(description='process Raman spectra')

	parser.add_argument("src", metavar='SOURCE', 
	help='either libs or tweez or opto', action='store')
	
	parser.add_argument('--sample', 
	help='average spectra by sample instead of by measurement', action='store_true')

	parser.add_argument("--bg", metavar='BGFILENAME', help='subtract BGFILENAME from spectra, exclude file extension and spectrum id in BGFILENAME', action='store')

	args = parser.parse_args()
	
	
# 	Provide feedback to user and establish arg rules
	print('')
	
	if args.bg:
		print('All spectra will be subtracted from indicated background.')
		print('')

	print('Each sample has several measurements which consist of several spectra.')

	if args.sample:
		print('Results will be averaged across each sample.')
		Suffix = re.compile('[-]\d+_\d+[.]\w+$') # Remove Measurement & Spectrum ID	
	else:
		print('Results will be averaged across each measurement.')
		Suffix = re.compile('_\d+[.]\w+$') # Remove Spectrum ID
		
	print('')
		
	if str.lower(args.src) == 'libs':
		ReadFolder = '/Volumes/TRANSFER/LIBS/'
		CalibFile = 'Pixel-Wavelength-LIBS.xls'
		Divide1 = 255
		Divide2 = 1725
		Pixels = 3648
		fitordr = 0
		SkipIdx = 0
		xlbl = 'Wavelength(nm)'
	elif str.lower(args.src) == 'tweez':
		ReadFolder = '/Volumes/TRANSFER/Raman/'
		CalibFile = 'Pixel-WaveNumber-Grating600-866.1nm.xls'
		# 	For SERS suggest Divide1=190, for normal Raman suggest Divide1=130
		Divide1 = 190
		Divide2 = 1090
		Pixels = 1340
		SkipIdx = 1
		xlbl = 'Wave number(1/cm)'
	elif str.lower(args.src) == 'opto':	
		ReadFolder = '/Volumes/TRANSFER/Raman/'
		CalibFile = 'Pixel-WaveNumber-Optofluidics.xls'
		Divide1 = 80
		Divide2 = 1330
		Pixels = 1340
		SkipIdx = 0
		fitordr = 5
		xlbl = 'Wave number(1/cm)'
	else:
		raise Exception('data source must be either libs or tweez or opto')


# 	Load specunits calibration from file
	specunits = np.loadtxt(CalibPath+CalibFile, delimiter="\t")


# 	Read in file names, Hierarchy: SampleID -> MeasurementID -> FileID
	FileList = list(openfdiag(DialogCaption='Select spectrum files', \
	ReadFolder=ReadFolder, FileFilter=None))


# 	Use natural sorting for file list so it looks neater
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	FileList.sort(key=alphanum_key)


# 	Identify file extension, typically .xls or .txt, check that all are the same
	FileExt = ''.join(re.findall('[.]\w+$',FileList[0]))   # Convert list to str
	if not all(FileExt in s for s in FileList):
		raise Exception('selected files must share a single file extension')


# 	Start timer to monitor program run time
	startTime = time.time()


# 	Begin filename processing
	RemoveDir = [s.replace(ReadFolder,'') for s in FileList]
	RemoveSuffix = [Suffix.sub('',s) for s in RemoveDir]
	Spectra = list(set(RemoveSuffix))


# 	For loop used to process files and return mean curve for visual check
	for i in range(len(Spectra)):
		print('')
		print('Processing:', Spectra[i])
				
# 		File names to pass to specproc are different for sample vs measurement averaging
		FileNames = [s for s in FileList if Spectra[i] in s]

		Curvefit = specproc(FileNames, Spectra[i])


	print('')			
	print('Elapsed time:', round(time.time()-startTime,1), 'seconds')


	# import ipdb; ipdb.set_trace() # Breakpoint	
