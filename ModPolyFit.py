"This module contains the background removal method detailed by Lieber et al."

# The lieberfit function returns the modified polynomial curvefit for the supplied data.
#
# 'Data' is the spectrum to be corrected, it contains a single vertical column of 
# intensity values.
# 'Order' is the order of the modified polynomial fit to be performed.
# This code uses a convergence criteria and can converge automatically.

def lieberfit(Data,Order=5):
	import numpy as np

	NewCurve = np.zeros(shape=(Data.shape[0]))
	OldCurve = np.array(Data)
	Diff = NewCurve-OldCurve
	Convergence = np.dot(Diff,Diff)

	m = 0
	while Convergence > 1: # Suggest setting convergence criteria == pixel resolution
		P = np.polyfit(range(len(Data)),OldCurve,Order)
		NewCurve = np.polyval(P,range(len(Data)))
		np.copyto(OldCurve, NewCurve, where = NewCurve < OldCurve)
		m+=1
		Diff = NewCurve - OldCurve
		Convergence = np.dot(Diff,Diff)
	#print('Iterations needed for convergence: ',m,sep='')

	CurveFit=np.copy(NewCurve)

	return (CurveFit)
	

#--------------------------------------------------------------------------------------	

def threepartfit(Data, Order=5):
	import numpy as np
	from scipy import signal
	import statsmodels.api as sm
		
	IdxLeft = int(len(Data)/10)
	IdxRight = int(len(Data)/10)
# 	IdxOuter = int(len(Data)/15)

	# Perform piecewise modified polynomial fit over middle 4/5 and 1/3 extreme ranges
	Curvefit = np.zeros(len(Data))
# 	Middle = np.zeros(len(Data))
	Curvefit[:IdxLeft] = lieberfit(Data[:IdxLeft],Order-1)
	Curvefit[-IdxRight:] = lieberfit(Data[-IdxRight:],Order-1)
#	Middle[IdxOuter:-IdxOuter] = lieberfit(Data[IdxOuter:-IdxOuter],Order)
	Curvefit[IdxLeft:-IdxRight] = lieberfit(Data[IdxLeft:-IdxRight],Order)

	# In overlapping regions select average of two curves
# 	np.copyto(Curvefit, Middle, where = abs(Data - Middle) < abs(Data - Curvefit))
			
	# Heavy smoothing of background curve to remove discontinuities at dividers.
	SmoothCurve=np.copy(Curvefit)
# 	SmoothCurve[int(IdxOuter/2):-int(IdxOuter/2)]=signal.savgol_filter(
# 	SmoothCurve[int(IdxOuter/2):-int(IdxOuter/2)],9,1)
	# Use locally-weighted linear regression to smooth curves over appropriate region
# 	lowess = sm.nonparametric.lowess
# 	StartSmooth = int(IdxLeft/2)
# 	x = range(StartSmooth,len(SmoothCurve))
# 	filt = lowess(Curvefit[StartSmooth:],x,frac=0.1)
# 	SmoothCurve[StartSmooth:] = filt[:,1]

	# Subtract background and shift above zero.
	SubtractedResult=Data-SmoothCurve
	while any(SubtractedResult<0):
		SubtractedResult+=10
		# print('Files remaining ')
		
# 	import ipdb; ipdb.set_trace() # Breakpoint	
		
	return(SubtractedResult, SmoothCurve)
