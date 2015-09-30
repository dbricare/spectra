"""
Module contains spectrum processing routines used on all files, consists of using Lieber method to remove background, smoothing with Savitzky-Golay and computing normalized and non-normalized results.
"""



""""
Begin function definitions
"""



#--------------------------------------------------------------------------------------
def onepartfit(Data, Order=5):

	"""
	Single, contiguous modified polynomial fit
	"""

	from ModPolyFit import lieberfit


# Perform modified polynomial fit over entire data range
	Curvefit = lieberfit(Data, Order)
	
	SmoothCurve = np.copy(Curvefit)		


# Subtract background and shift above zero
	np.copyto(SmoothCurve, Data-10, where = SmoothCurve > Data)

	SubtractedResult=Data-SmoothCurve

	while any(SubtractedResult<0):
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
	np.copyto(SmoothCurve, Data-10, where = SmoothCurve > Data)
	SubtractedResult=Data-SmoothCurve
	while any(SubtractedResult<0):
		SubtractedResult+=10


	SmoothCurve[StartSmooth:]=signal.savgol_filter(SmoothCurve[StartSmooth:],11,1)	
		
		
	return(SubtractedResult, SmoothCurve)



#--------------------------------------------------------------------------------------
def specproc(FileNameList, FileBase):

	"""
	Heavy lifting for spectrum processing
	"""

# Remove files (typically SkipIdx=1, the first file) due to instrument issues
	if SkipIdx != 0 and SkipIdx > 0:
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


	NumberFiles=len(SkippedList)
		
		
# Initiate matrices for spectrum and wavenumber data.
	FullName = ['' for i in range(NumberFiles)]
	Raw = np.zeros((Pixels,2,NumberFiles))   # 3d array
	SmoothData = np.zeros((Pixels,NumberFiles))   # 2d arrays here on
	Curve = np.zeros((Pixels,NumberFiles))
	SmoothCurve = np.zeros((Pixels,NumberFiles))
	SubtractedResult = np.zeros((Pixels,NumberFiles))


# Background file processing
	BGsub = np.zeros(Pixels)
	if args.bg:
		BGfileCount = 3
		BGspectrum = np.zeros((Pixels,BGfileCount))
		for i in range(BGfileCount):
			BGdata = np.loadtxt(ReadFolder+args.bg+'_'+str(i+1+SkipIdx)+FileExt)
			BGspectrum[:,i] = BGdata[:,1]
		BGsub = np.mean(BGspectrum, axis=1)
			
			
# Load files and perform background fit
	for i in range(NumberFiles):
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
	
		if args.fit:
			SubtractedResult[:,i], SmoothCurve[:,i] = onepartfit(SmoothData[:,i],fitorder)
		else:
			SubtractedResult[:,i], SmoothCurve[:,i] = threepartfit(SmoothData[:,i],Order)

	MeanSmooth=np.mean(SmoothData,axis=1)
	MeanCurve=np.mean(SmoothCurve,axis=1)
	
	
# 	Perform peak fitting, mean, std, and output all to tab-delimited files

	def outputRslt(SaveExt):
		
		if SaveExt == '-SF':   # No normalization
			RsltData = np.copy(SubtractedResult)
		elif SaveExt == '-SFN':   # p-norm normalization
			RsltData = normalize(SubtractedResult,norm='l1',axis=0)
		elif SaveExt == '-SFM':   # Min-max normalization
			MMscale = MinMaxScaler()
			RsltData = MMscale.fit_transform(SubtractedResult)
		else:
			raise Exception('Unknown processing method specified by programmer')

# 	Calculate mean & std dev of smoothed/fitted data and background removal curves for visualization
		Mean = np.mean(RsltData,axis=1)
		Std = np.std(RsltData,axis=1)

#	Find peaks of mean data with built in function that fits wavelets to data
		PeakWidths=np.arange(1,30)
		PkIdx = signal.find_peaks_cwt(Mean,PeakWidths)

# 	Compute peak locations and intensities for non-normalized and normalized results
		PkLoc = WaveNumber[PkIdx]
		PkInt = Mean[PkIdx].reshape((len(PkIdx),1))
	
# 	Extend peak data to the same length as other data in order to add to file output then sort in descending order.
		WritePeak = np.zeros((Pixels,2))
		WritePeak[:len(PkIdx),:] = np.hstack((PkLoc, PkInt))
		WritePeak = np.flipud(WritePeak[WritePeak[:,1].argsort()])

# 	Save results in a tab-delimited file with a .xls extension for import into Excel	
		Wname = WriteFolder+FileBase+SaveExt+'.xls'
		Wdata = np.hstack((WaveNumber, Mean.reshape((Pixels,1)), Std.reshape((Pixels,1)), WritePeak, RsltData))
		np.savetxt(Wname, Wdata, fmt='%g', delimiter='\t', header=Wheader, comments='')
	
		if SaveExt == '-SF':
	# Visualize results
			plt.rc('font', family = 'Arial', size='14')
			plt.ion()
	
			plt.figure()
			plt.plot(WaveNumber,Mean,color='blue')
			plt.plot(PkLoc,PkInt+max(Mean)*0.03,"o",color="green",markersize=6)
			plt.title(FileBase+' - Subtracted with peaks')

			plt.figure()
			plt.plot(WaveNumber,MeanCurve,color='red')
			plt.plot(WaveNumber,MeanSmooth,color='blue')
			plt.title(FileBase+' - Smoothed & curve fit')
			if not args.fit:
				plt.axvline(x=WaveNumber[Divide1],color='k',ls='--',lw=1)
				plt.axvline(x=WaveNumber[Divide2],color='k',ls='--',lw=1)
			else:
				plt.axvline(x=WaveNumber[0],color='k',ls='--',lw=1)
				plt.axvline(x=WaveNumber[-1],color='k',ls='--',lw=1)

	
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
	ReadFolder = '/Volumes/TRANSFER/Raman/'
	WriteFolder = '/Volumes/TRANSFER/Analysis/'
	CalibPath = '/Users/dbricare/Desktop/Chan Lab/Labwork/'
# 	CalibFile = 'Pixel-Wavenumber-Grating600-866.1nm.xls'
	CalibFile = 'Pixel-Wavenumber-Optofluidics.xls'
	SkipIdx = 1
	Order = 5
	Pixels = 1340
# 	For SERS suggest Divide1=190, for normal Raman suggest Divide1=130
	Divide1 = 80
	Divide2 = 1330
# 	StartSmooth = 50


# 	Import libraries
	import re, argparse, time
	from ui_python.FileOpen import openFileDialog
	import numpy as np
	from scipy import signal
	import matplotlib.pyplot as plt
	from sklearn.preprocessing import MinMaxScaler, normalize


# 	Parse the arguments passed by the user
	parser = argparse.ArgumentParser(description='process Raman spectra')

	parser.add_argument('--sample', \
	help='average spectra by sample instead of by measurement', action='store_true')

	parser.add_argument("--bg", help='subtract BGFILENAME from spectra, exclude file extension and spectrum id', action='store')
	
	parser.add_argument("--fit", help='Perform one and three part modified polynomial fit', action='store')

	args = parser.parse_args()
	
	
# 	Provide feedback to user and set
	print('')
	
	if args.bg:
		print('All spectra will be subtracted from indicated background.')

	if args.sample:
		print('Spectra will be averaged across each SAMPLE.')
		Suffix = re.compile('[-]\d+_\d+[.]\w+$') # Remove Measurement & Spectrum ID	
	else:
		print('Spectra will be averaged across each MEASUREMENT.')
		Suffix = re.compile('_\d+[.]\w+$') # Remove Spectrum ID
		
	if args.fit:
		print('Performing one part polynomial fit with indicated order.')
		fitorder = int(args.fit)
	else:
		print('Performing three part fit with default polynomial orders.')	
	


# 	Load wavenumber calibration from file
	WaveNumber = np.loadtxt(CalibPath+CalibFile, delimiter="\t")
	WaveNumber = WaveNumber.reshape((Pixels,1))
	Wheader='Wave number\t'+'Average\t'+'Std dev\t'+'Peak locations\t'+ \
	'Peak intensity\t'+'Smoothed data'	


# 	Read in file names, Hierarchy: SampleID -> MeasurementID -> FileID
	FileList = list(openFileDialog(DialogCaption='Select spectrum files', \
	ReadFolder=ReadFolder, FileFilter=None))
	# Sort file list
	FileList.sort()


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
		
		FileNames = [s for s in FileList if Spectra[i] in s]
	
		MeanCurve = specproc(FileNames, Spectra[i])


	print('')			
	print('Elapsed time:', round(time.time()-startTime,1), 'seconds')


	# import ipdb; ipdb.set_trace() # Breakpoint	
