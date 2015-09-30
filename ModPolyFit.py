"""
This module contains the background removal method detailed by Lieber et al.
"""


# The lieberfit function returns the modified polynomial curvefit for the supplied data.
#
# 'Data' is the spectrum to be corrected, it contains a single vertical column of 
# intensity values.
# 'Order' is the order of the modified polynomial fit to be performed.
# This code uses a convergence criteria and can converge automatically.



#--------------------------------------------------------------------------------------
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
	
