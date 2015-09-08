"Module contains spectrum processing routines used on all files, consists of using Lieber method to remove background, smoothing with Savitzky-Golay and computing normalized and non-normalized results."


# Define a few global variables that may need infrequent modification
ReadFolder = '/Volumes/TRANSFER/Raman/'
WriteFolder = '/Volumes/TRANSFER/Analysis/'
SkipIdx = 1
Order = 5
Pixels = 1340
# For SERS suggest Divide1=190, for normal Raman suggest Divide1=130
Divide1 = 90
Divide2 = 1090
# StartSmooth = 50


# Parse the arguments passed by the user
import argparse
parser = argparse.ArgumentParser(\
description='to process spectra collected with the Raman trap')

parser.add_argument('--sample', \
help='average spectra by sample instead of by measurement', action='store_true')

parser.add_argument("--bg", metavar='BGFILENAME', help='subtract BGFILENAME from spectra, exclude file extension and spectrum id', action='store')

args = parser.parse_args()

print('')
if args.bg:
	print('All spectra will be subtracted from indicated background.')

if args.sample:
	print('Spectra will be averaged across each SAMPLE.')
else:
	print('Spectra will be averaged across each MEASUREMENT.')
	
	
Avg = args.sample
bg = args.bg


## Start function definitions
#
#--------------------------------------------------------------------------------------

def threepartfit(Data, Order=5):

	import statsmodels.api as sm
	from ModPolyFit import lieberfit
	
		
# 	IdxLeft = int(len(Data)/20)
# 	IdxRight = int(len(Data)/10)
# 	IdxOuter = int(len(Data)/15)
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
	StartSmooth = int(IdxLeft-10)
	# Use locally-weighted linear regression to smooth curves over appropriate region
	lowess = sm.nonparametric.lowess
	x = range(StartSmooth,len(SmoothCurve))
	filt = lowess(Curvefit[StartSmooth:],x,frac=0.1)
	SmoothCurve[StartSmooth:] = filt[:,1]


# Subtract background and shift above zero
	np.copyto(SmoothCurve, Data-10, where = SmoothCurve > Data)
	SubtractedResult=Data-SmoothCurve
	while any(SubtractedResult<0):
		SubtractedResult+=10


	SmoothCurve[StartSmooth:]=signal.savgol_filter(SmoothCurve[StartSmooth:],11,1)	
		
		
	return(SubtractedResult, SmoothCurve)



#--------------------------------------------------------------------------------------
# This function does the bulk of the spectrum processing

def specproc(FileNameList, FileBase):

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
	if bg:
		BGfileCount = 3
		BGspectrum = np.zeros((Pixels,BGfileCount))
		for i in range(BGfileCount):
			BGdata = np.loadtxt(ReadFolder+bg+'_'+str(i+1+SkipIdx)+FileExt)
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
	
		SubtractedResult[:,i], SmoothCurve[:,i] = threepartfit(SmoothData[:,i], Order)

	MeanSmooth=np.mean(SmoothData,axis=1)
	MeanCurve=np.mean(SmoothCurve,axis=1)
	
	
	# Perform peak fitting, mean, std, and output all to tab-delimited files

	def outputRslt(SaveExt):
		
		if SaveExt == '-SF':   # No normalization
			RsltData = np.copy(SubtractedResult)
		elif SaveExt == '-SFN':   # p-norm normalization
			RsltData = normalize(SubtractedResult,norm='l1',axis=0)
# 			import ipdb; ipdb.set_trace() # Breakpoint	
		elif SaveExt == '-SFM':   # Min-max normalization
			MMscale = MinMaxScaler()
			RsltData = MMscale.fit_transform(SubtractedResult)
		else:
			raise NameError('Unknown processing method specified')

	# Calculate mean & std dev of smoothed/fitted data and background removal curves for visualization
		Mean = np.mean(RsltData,axis=1)
		Std = np.std(RsltData,axis=1)

	# Find peaks of mean data with built in function that fits wavelets to data
		PkIdx = signal.find_peaks_cwt(Mean,PeakWidths)

	# Compute peak locations and intensities for non-normalized and normalized results
		PkLoc = WaveNumber[PkIdx]
		PkInt = Mean[PkIdx].reshape((len(PkIdx),1))
	
	# Extend peak data to the same length as other data in order to add to file output then sort in descending order.
		WritePeak = np.zeros((Pixels,2))
	# 	import ipdb; ipdb.set_trace() # Breakpoint	
		WritePeak[:len(PkIdx),:] = np.hstack((PkLoc, PkInt))
		WritePeak = np.flipud(WritePeak[WritePeak[:,1].argsort()])

	# Save results in a tab-delimited file with a .xls extension for import into Excel	
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
			plt.axvline(x=WaveNumber[Divide1],color='k',ls='--',lw=1)
			plt.axvline(x=WaveNumber[Divide2],color='k',ls='--',lw=1)
			plt.title(FileBase+' - Smoothed & curve fit')
		
	NoDisplay = [outputRslt(s) for s in ['-SFM','-SFN','-SF']]
	
	return(MeanCurve)


#--------------------------------------------------------------------------------------
# Main

import re
from FileOpen import openFileDialog
import time
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler, normalize


# Set a few more constant values
PeakWidths=np.arange(1,30)
WaveNumber = np.loadtxt('Pixel-Wavenumber-Grating600-866.1nm.xls', delimiter="\t")
WaveNumber = WaveNumber.reshape((Pixels,1))
Wheader='Wave number\t'+'Average\t'+'Std dev\t'+'Peak locations\t'+ \
'Peak intensity\t'+'Smoothed data'	


# Read in file names, Hierarchy: SampleID -> MeasurementID -> FileID
FileList = list(openFileDialog(DialogCaption='Select spectrum files', \
ReadFolder=ReadFolder, FileFilter=None))
# Sort file list
FileList.sort()


# Identify file extension, typically .xls or .txt, check that all are the same
#reExt = re.compile('[.]\w+$')
global FileExt
FileExt = ''.join(re.findall('[.]\w+$',FileList[0]))   # Convert list to str
if not all(FileExt in s for s in FileList):
	raise Exception('the selected files do not share a single file extension')


# Start timer to track program run time
startTime = time.time()


# Identify sample and measurement/spectrum IDs
if Avg:
	Suffix = re.compile('[-]\d+_\d+[.]\w+$') # Remove Measurement & Spectrum ID		
else:
	Suffix = re.compile('_\d+[.]\w+$') # Remove Spectrum ID


RemoveDir = [s.replace(ReadFolder,'') for s in FileList]
RemoveSuffix = [Suffix.sub('',s) for s in RemoveDir]
Spectra = list(set(RemoveSuffix))


for i in range(len(Spectra)):
	print('')
	print('Processing:', Spectra[i])
	#print('')
	FileNames = [s for s in FileList if Spectra[i] in s]
	
	
	MeanCurve = specproc(FileNames, Spectra[i])


print('')			
print('Elapsed time:', round(time.time()-startTime,1), 'seconds')


# import ipdb; ipdb.set_trace() # Breakpoint	
