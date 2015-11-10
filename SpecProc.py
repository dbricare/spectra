"""
Module contains spectrum processing routines used on all files, consists of using Lieber method to remove background, smoothing with Savitzky-Golay and computing normalized and non-normalized results.

Assumes tab-delimited files where first column is spectrum units and second column is intensity values.
"""



""""
Begin function definitions
"""

#--------------------------------------------------------------------------------------
# from numba import jit
# jit decorator tells Numba to compile this function.
# The argument types will be inferred by Numba when function is called.
# @jit
def modpolyfit(Data,Order=5):

	"""
This function contains the background removal method detailed by Lieber et al. and returns the modified polynomial curvefit for the supplied data.

'Data' is the spectrum to be corrected, it contains a vector of intensity values. 
'Order' is the order of the modified polynomial fit to be performed.

This code uses a convergence criteria and can converge automatically.
	"""

	NewCurve = np.ones(Data.shape[0])*(Data.min()-1)

	if Order > 0:
		OldCurve = np.copy(Data)
		Diff = OldCurve - NewCurve
		Convergence = np.dot(Diff,Diff)

		m = 0
		while Convergence > 1: # Suggest convergence criteria == intensity resolution
			P = np.polyfit(range(len(Data)),OldCurve,Order)
			NewCurve = np.polyval(P,range(len(Data)))
# 			Replace OldCurve values with NewCurve values that are smaller
			np.copyto(OldCurve, NewCurve, where = NewCurve < OldCurve)
			m+=1
# 			Stop when the squared difference is below convergence criteria
			Diff = OldCurve - NewCurve
			Convergence = np.dot(Diff,Diff)
# 		print('Iterations needed for convergence: ',m,sep='')


	CurveFit=np.copy(NewCurve)


	return (CurveFit)
	

#--------------------------------------------------------------------------------------
def onepartfit(Data, Order=0):

	"""
	Single, contiguous modified polynomial fit
	"""

# Perform modified polynomial fit over entire data range
	Curvefit = modpolyfit(Data, Order)
	
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

# Perform piecewise modified polynomial fit over middle 4/5 and 1/3 extreme ranges
	Curvefit = np.zeros(len(Data))
	Curvefit[:divide1] = modpolyfit(Data[:divide1],Order-1)
	Curvefit[divide2:] = modpolyfit(Data[divide2:],Order-1)
	Curvefit[divide1:divide2] = modpolyfit(Data[divide1:divide2],Order)
			
			
# Heavy smoothing of background curve to remove discontinuities at dividers
	SmoothCurve=np.copy(Curvefit)
	StartSmooth = max(int(divide1-10),0)
	TruncLeft = max(StartSmooth,divide1)
# 	Use locally-weighted linear regression to smooth curves over appropriate region
	lowess = sm.nonparametric.lowess
	x = range(TruncLeft,pixels)
	filt = lowess(Curvefit[TruncLeft:],x,frac=0.1)
	SmoothCurve[TruncLeft:] = filt[:,1]


# 	Subtract background and shift above zero
	np.copyto(SmoothCurve, Data, where = SmoothCurve > Data)
	
	SubtractedResult=Data-SmoothCurve
	
	while any(SubtractedResult<=0):
		SubtractedResult+=10

		
	return(SubtractedResult, SmoothCurve)



#--------------------------------------------------------------------------------------
def specproc(FileNameList, FileBase):

	"""
	Heavy lifting for spectrum processing
	"""
	import pandas as pd

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
		
		
# Initiate matrices for spectrum data.
	FullName = ['' for i in range(numfiles)]
	SmoothData = np.zeros((pixels,numfiles))   # 2d arrays here on
	Curve = np.zeros((pixels,numfiles))
	SmoothCurve = np.zeros((pixels,numfiles))
	SubtractedResult = np.zeros((pixels,numfiles))
	
	if args.src == 'libs':
		Raw = np.zeros((fullrng,2,numfiles))   # 3d array
	else:
		Raw = np.zeros((pixels,2,numfiles))   # 3d array


# Background file processing
	if args.bg:
		BGfileCount = 3   # must be set manually
		BGspectrum = np.zeros((pixels,BGfileCount))
		for i in range(BGfileCount):
			BGdata = np.loadtxt(ReadFolder+args.bg+'_'+str(i+1+SkipIdx)+FileExt)
			BGspectrum[:,i] = BGdata[:,1]
		BGsub = np.mean(BGspectrum, axis=1)
	else:
		BGsub = np.zeros(pixels)
			
			
# Load files and perform background fit
	for i in range(numfiles):
		Raw[:,:,i]=np.loadtxt(SkippedList[i],delimiter="\t")
	
# Check for instrument error where consecutive spectra share many identical values
		if i > 0: 
			DupChk = sum(Raw[:,1,i] == Raw[:,1,i-1])   # integer-valued, tol=1
			if DupChk > pixels/2:
				print('')
				import warnings
				warnmsg='Possible duplicate spectra, {:} identical pixels'.format(DupChk)
				warnings.warn(warnmsg, UserWarning)
				print('Files:')
				print(SkippedList[i-1])
				print(SkippedList[i])
		
		
# Apply Savitzky-Golay filtering (11-point window, 1st or 2nd order polynomial)
# Note: larger frames and lower order polynomial fits have stronger smoothing
		if args.src == 'libs':
			SmoothData[:,i] = signal.savgol_filter(Raw[divide1:divide2,1,i] -
			 0.8*BGsub,11,1) 
		else:
			SmoothData[:,i] = signal.savgol_filter(Raw[:,1,i] - 0.8*BGsub,11,1) 


# Perform modified polynomial fit
		if str.lower(args.src)=='tweez':
			import statsmodels.api as sm
			SubtractedResult[:,i],SmoothCurve[:,i] = threepartfit(SmoothData[:,i],fitordr)
		else:   # args.src is libs or opto
			SubtractedResult[:,i],SmoothCurve[:,i] = onepartfit(SmoothData[:,i],fitordr)

	MeanSmooth=np.mean(SmoothData,axis=1)
	MeanCurve=np.mean(SmoothCurve,axis=1)
	
# 	
# 	Perform peak fit, mean, std, and output all to xlsx files using loop
# 	
	SaveExts = ('-SFP','-SFM','-SF')  # p-norm (l1), min-max and none normalizations

	mmscale = MinMaxScaler()

	RsltDatas = (normalize(SubtractedResult,norm='l1',axis=0),
	mmscale.fit_transform(SubtractedResult), 
	np.copy(SubtractedResult))

	for SaveExt, RsltData in zip(SaveExts, RsltDatas):

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


# 		Use pandas to save results as an Excel file, include all processed spectra
		dfrslt = pd.DataFrame(data=RsltData, columns=currfiles)

		dfwrite = pd.concat([dfstat, dfpks, dfrslt], axis=1)
		namewrite = WriteFolder+FileBase+SaveExt+'.xlsx'
		dfwrite.to_excel(namewrite, header=True, index=False)
	
	
# 		Visualize results just once not three times for each normalization
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
			if args.src == 'tweez':
				plt.axvline(x=specunits[divide1],color='k',ls='--',lw=1)
				plt.axvline(x=specunits[divide2],color='k',ls='--',lw=1)

# 			Create min/max plot, first row is not shown in plot (fn expects header row)
# 			from plotting.makeplot import fillbtwn
# 			plt.ioff()
# 			fillbtwn(Wdata[:,:4], savename=FileBase)
# 			fillbtwn(Wdata[:,:4], savename=None)
# 			plt.ion()
	
	
	return(MeanCurve)


"""
End function definitions
"""


#--------------------------------------------------------------------------------------
# Main

# As of 2015-09-29, this file is not configured for use as an importable module, functions defined above call global variables from __main__
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


# 	Parse the arguments passed by the user
	parser = argparse.ArgumentParser(description='process Raman spectra')

	parser.add_argument("src", metavar='SOURCE', type=str.lower, 
	help='either libs or tweez or opto', action='store')
	
	parser.add_argument('--sample', 
	help='average spectra by sample instead of by measurement', action='store_true')

	parser.add_argument("--bg", metavar='BGFILENAME', help='subtract BGFILENAME from spectra, exclude file extension and spectrum id in BGFILENAME', action='store')

	args = parser.parse_args()
	
		
	if args.src == 'libs':
		ReadFolder = '/Volumes/TRANSFER/LIBS/'
		CalibFile = 'Pixel-Wavelength-LIBS.xls'
		divide1 = 255
		divide2 = 1725
		pixels = divide2-divide1
		fullrng = 3648
		SkipIdx = 0
		fitordr = 0
		xlbl = 'Wavelength(nm)'
	elif args.src == 'tweez':
		ReadFolder = '/Volumes/TRANSFER/Raman/'
		CalibFile = 'Pixel-WaveNumber-Grating600-866.1nm.xls'
		# 	For SERS suggest divide1=190, for normal Raman suggest divide1=130
		divide1 = 190
		divide2 = 1090
		pixels = 1340
		SkipIdx = 1
		xlbl = 'Wave number(1/cm)'
	elif args.src == 'opto':	
		ReadFolder = '/Volumes/TRANSFER/Raman/'
		CalibFile = 'Pixel-WaveNumber-Optofluidics.xls'
		divide1 = 80
		divide2 = 1330
		pixels = 1340
		SkipIdx = 0
		fitordr = 5
		xlbl = 'Wave number(1/cm)'
	else:
		raise Exception('data source must be either libs or tweez or opto')


# 	Provide feedback to user and establish arg rules
	print('')
	
	print('File organization hierarchy: Sample -> Measurement -> Spectra.')
	
	if args.bg:
		print('')
		print('All spectra will be subtracted from indicated background.')

	print('')

	if args.sample:
		print('Results will be averaged across each sample.')
		Suffix = re.compile('[-]\d+_\d+[.]\w+$') # Remove Measurement & Spectrum ID	
	else:
		print('Results will be averaged across each measurement.')
		Suffix = re.compile('_\d+[.]\w+$') # Remove Spectrum ID
		
	print('')
	

# 	Load spectral units calibration from file (wavelength or wavenumber)
	specunits = np.loadtxt(CalibPath+CalibFile, delimiter="\t")
	if args.src == 'libs':
		specunits = specunits[divide1:divide2]


# 	Read in file names, Hierarchy: SampleID -> MeasurementID -> FileID
	print('Select files for analysis from dialog window...')
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
	for spectrum in Spectra:
		print('')
		print('Processing:', spectrum)
				
# 		File names to pass to specproc are different for sample vs measurement averaging
		FileNames = [s for s in FileList if spectrum in s]

		Curvefit = specproc(FileNames, spectrum)


	print('')			
	print('Elapsed time:', round(time.time()-startTime,1), 'seconds')


	# import ipdb; ipdb.set_trace() # Breakpoint	
