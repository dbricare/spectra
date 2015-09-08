# Commands used to load in data, organize, perform LDA classification of 2D PCA of data

import numpy as np
import matplotlib.pyplot as plt
import sklearn, itertools
from sklearn import lda, preprocessing, cross_validation, cross_decomposition
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.metrics import f1_score, accuracy_score, mean_squared_error
from sklearn.cross_decomposition import PLSRegression, PLSCanonical, PLSSVD
from FileOpen import openFileDialog

# Global variables that may need modification
PixelStart = 220 #145 #220
PixelEnd = 1104 #930 #1104
ClassSubDiv = 2

minPLS = 2
maxPLS = 20
intervalPLS = 1


#--------------------------------------------------------------------------------------
# Create label names
def genLblList(FileList):
	RemoveDir = [s.replace(ReadFolder,'') for s in FileList]
	Suffix = '-LNS.xls'
	LblList = [s.replace(Suffix,'') for s in RemoveDir]
	return(LblList)


#--------------------------------------------------------------------------------------
# Define function to load files, remove first 5 columns and organize data
def loaddata(FileList, LblList):
	ExampleClasses = []
	ClassInts = range(len(LblList))
	for i in range(len(FileList)):
		Temp = np.loadtxt(FileList[i],delimiter='\t',skiprows=1)
		Temp = Temp[:,5:]
		if i == 0:
			Store = np.array(Temp.T)
			ExamplesPerClass = Temp.shape[1]
		else:
			Store = np.vstack((Store,Temp.T))
		CurrInt = [ClassInts[i]]*ExamplesPerClass
		ExampleClasses.extend(CurrInt)
	return(Store, ExampleClasses)


#--------------------------------------------------------------------------------------
# Variable importance in projection
# VIP is not appropriate if number of features >> number of examples

def pls1vip():
	# PLS 1 (one-block) regression VIP here y is a vector
	# Calculation is from 2005 paper by Chong and Jun, Performance of some variable selection methods when multicollinearity is present
	# PLS regression with h latent variables expressed as:
	# X = T*P.T + E, y = T*b + f
	# T is X scores, P is X loadings, b is y loadings

	b = pls.y_loadings_
	T = pls.x_scores_
	h = pls.n_components
	W = pls.x_weights_
	p = pls.x_weights_.shape[0]   # number of features/original variables

	# Calculation of SS vector for all h latent variable, vector with length = # components
	SS = np.dot(b**2,np.dot(T.T,T))

	# Square of the quantity W divided by norm of each column but first create norm vector
	NormVec = np.linalg.norm(W,2,axis=1).reshape(-1,1)
	Square = np.power(np.divide(W,NormVec),2)

	# Compute VIP
	Num = np.multiply(Square,SS.reshape(1,-1))
	Vip = np.sqrt(p / np.sum(SS) * np.sum(Num, axis=1))
	
	return(Vip)


#--------------------------------------------------------------------------------------
# Visualization

def vis():
	plt.ion()
	plt.rc('font', family = 'Arial', size='18')
	
	nClss = len(LblList)

	# For customized labels not derived from file names
	#DataLbls = ['Mutant', 'Parent']
# 	DataLbls = ['X', 'Y']

	#LblList = ['Supernatant - no wash', 'Supernatant - wash 4x', 'Citrate']
	#PCAsort = np.zeros(np.shape(PCACoeff[:,:2]))
	#CustomClss = []
	#for i in range(nClss):
	#	CustomClss.extend([i]*int(len(Data[:,0])/nClss))
	#y = CustomClss

	# Custom colors and markers
#	ColorOpts = ['CornflowerBlue', 'Crimson', 'DarkSeaGreen', 'DarkOrchid', 'Coral', \
#	'Gold', 'LightSeaGreen', 'BurlyWood', 'SaddleBrown']
	ColorOpts = ['CornflowerBlue','Coral', 'MediumPurple', 'DarkCyan', 'Gold',  'Crimson', 'YellowGreen', 'BurlyWood', 'SaddleBrown']
	ColorSeq = ColorOpts[:int(len(LblList)/ClassSubDiv)]*ClassSubDiv
	MarkerOpts = ['o','s','^','D']
	#MarkerSeq = MarkerOpts[:ClassSubDiv-1]


	### Plot first PLS component, X and Y scores
# 	fig = plt.figure(figsize=(12,8))
# 	ax = plt.subplot(111)
# 	j = 0
# 
# 	# Set range(nComp/nComp) to plot first component only
# 	plsScores = np.hstack((Xpls,Ypls))
# 	for i,c,class_name in zip(range(int(nComp/nComp)),ColorSeq,LblList):
# 		ax.scatter(plsScores[:,i],plsScores[:,i+nComp], c=c, label=class_name, s=100, marker=MarkerOpts[j])
# 		if i >= ((j+1)*len(y)/ClassSubDiv)-1:
# 			j+=1
# 	# Add trendline to aid visualization of correlation between X & Y scores
# 	xmin, xmax, _, _ = plt.axis()
# 	Xtrend = range(int(xmin),int(xmax))
# 	ax.plot(Xtrend,Xtrend,'k-')
# 
# 	# Shrink current axis by 20% and put legend to right of current axis
# # 	box = ax.get_position()
# # 	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# # 	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1)
# 
# 	if i != 0:
# 		plt.legend(loc='best', scatterpoints=3)
# 
# 	plt.ylabel(DataLbls[1]+' Scores')   # Y scores
# 	plt.xlabel(DataLbls[0]+' Scores')   # X scores


	### Plot first two discriminant axes from LDA
# 	fig = plt.figure(figsize=(14,8))
# 	ax = plt.subplot(111)
# 	j = 0
# 	for i,clr,class_name in zip(range(nClss),ColorSeq,LblList):
# 		IdxLda = np.where(np.array(ExampleClasses)==i)
# 		ax.scatter(Xlda[IdxLda,0], Xlda[IdxLda,1], c=clr, label=class_name, s=100, marker=MarkerOpts[j])
# 		if i >= ((j+1)*int(nClss/ClassSubDiv))-1:
# 			j+=1
# 
# 	# Shrink current axis by 20% and put legend to right of current axis
# 	box = ax.get_position()
# 	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# 	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1)
# 
# 	plt.ylabel('Discriminant Axis 2')
# 	plt.xlabel('Discriminant Axis 1')
# 	plt.title('Two Component LDA with {:} PLS Components'.format(nPLS))
# # 	plt.legend(loc='best')
	

	### Plot PLS2 MSEP as a function of number of PLS components
	### Plot PLS2 MSEP from k-folds CV as a function of number of PLS components
	### Plot LDA misclassification as a function of number of PLS components
	plt.figure(figsize=(8,6))
	plt.plot(PLSrange,StoreAccLDA,'b-',label='LDA misclassification')
	plt.plot(PLSrange,StoreAccPLS,'g-',label='MSEP of PLS2')
	plt.plot(PLSrange,StorekAccPLS,'r-',label='MSEP from k-folds CV of PLS2')
	plt.xlabel('Number of PLS components')
	plt.ylabel('Error')
	plt.legend(loc='upper right')
# 	plt.ylim(0,1)

# 	### Plot Weights squared for first component
	plt.figure(figsize=(9,6))
	plt.plot(WaveNumber,pls.x_weights_[:,0]**2,'b-')
	plt.xlabel(r'Wave Number ($\mathregular{cm}^{-1}$)')
	plt.ylabel('Weight squared for first PLS component')


# 	### Plot Variables important in projection
# 	plt.figure(figsize=(8,6))
# 	plt.plot(WaveNumber,Vip,'b-')
# 	plt.xlabel(r'Wave Number ($\mathregular{cm}^{-1}$)')
# 	plt.ylabel('Variable importance in projection')
	

#--------------------------------------------------------------------------------------
# Main
#
# Read in file names, files must have -LNS label indicating normalization
DialogCaption = 'Select spectrum files containing (-LNS.xls)'
ReadFolder = '/Volumes/TRANSFER/Analysis/'
FileFilter = '*-LNS.xls'
FileList = list(openFileDialog(DialogCaption, ReadFolder, FileFilter))

# Sort file list, generate label names and integers corresponding to each examples class
# The number of unique elements in ExamplesClasses == len(LblList)
FileList.sort()
LblList = genLblList(FileList)
Data, ExampleClasses = loaddata(FileList,LblList)

# Trim data window if necessary
WaveNumber = np.loadtxt('Pixel-Wavenumber-Grating600-866.1nm.xls',delimiter="\t",ndmin=2)
WaveNumber = WaveNumber[PixelStart:PixelEnd]
Data = Data[:,PixelStart:PixelEnd]

# Additionally, it may be possible to remove collinear features from the data and run LDA directly on the reduced data set. Pearson's correlation coefficient may be calculated.
	
# Perform PLS
Data = preprocessing.StandardScaler().fit_transform(Data)   # mean=0, variance=1
#Xdata = Data[:(Data.shape[0]/2),:]
#Ydata = Data[(Data.shape[0]/2):,:]

# Create matrix of dummy variables with class for each dependent variable (output)

# LblList = ['Mutant','Parent']   # Put non-filename derived labels here

m = Data.shape[0]
n = len(LblList)
y = -1*np.ones((m,n),dtype='float')   # in-class examples should be +1, out -1
ExampleClasses = []
for i in range(len(LblList)):
	y[range(i*int(m/n),(i+1)*int(m/n)),i]+=2
	ExampleClasses.extend([i]*int(m/n))   # 1D version of dummy variable matrix y
ExampleClasses = np.array(ExampleClasses,dtype='int')

y = np.copy(ExampleClasses)   # Uncomment for PLS1

PLSrange = intervalPLS*np.arange(minPLS/intervalPLS,1+maxPLS/intervalPLS, dtype='int')
StoreAccLDA = np.zeros((1+(maxPLS-minPLS)/intervalPLS,1))
StoreAccPLS = np.zeros((1+(maxPLS-minPLS)/intervalPLS,1))
StorekAccPLS = np.zeros((1+(maxPLS-minPLS)/intervalPLS,1))

# Fit LDA for a varying range of PLS components up to maxPLS
for nPLS in PLSrange:
	pls = PLSRegression(n_components=nPLS)
	TrnsfrmPls = pls.fit(Data,y).transform(Data,y)


	# Compute correlation between two datasets for first and second 
	#Xpls = pls.x_scores_
	#Ypls = pls.y_scores_
	#CorrCoef = np.corrcoef(Xpls,Ypls,rowvar=0)
	#print('')
	#print('Correlation between the two datasets in component 1: {:.3}'.format(CorrCoef[2,0]))
	#print('Correlation between the two datasets in component 2: {:.3}'.format(CorrCoef[1,3]))


	### Determine cross-validation scores using k-folds repeated n_iter times with a new random sorting
	cvPLS = cross_validation.StratifiedShuffleSplit(y, n_iter=10, test_size=0.2, random_state=None)   # Stratified k-folds of 1/test_size or 5 typically


	### Find CV scores using root means square error for PLS to help determine appropriate number of components
	print('')

	predPLS = np.array(pls.predict(Data), dtype='int')

	msepPLS = mean_squared_error(predPLS,y)
	print('PLS MSEP with {:} PLS components: {:.2e}'.format(nPLS, msepPLS))

	msePLSScores = cross_validation.cross_val_score(
	pls, Data, y, cv=cvPLS, scoring='mean_squared_error') # bug- returns negative values
	print('k-folds PLS MSEP: {:.2e}'.format(abs(np.mean(msePLSScores))))


	### Perform classification then transform PLS data to LDA basis
	nLDA = 2
	clfLDA = lda.LDA(n_components = nLDA)
	Xlda = clfLDA.fit_transform(TrnsfrmPls[0],ExampleClasses)
	
	# Predict and calculate misclassification rate
	cvLDA = cross_validation.StratifiedShuffleSplit(ExampleClasses, n_iter=10, test_size=0.2, random_state=None)   # Stratified k-folds of 1/test_size or 5 typically
	predLDA = np.array(clfLDA.predict(TrnsfrmPls[0]), dtype='int')
	AccLDAmean = np.mean(abs(predLDA-ExampleClasses))
# 	print('')
	print('LDA misclassification rate with {:} PLS components: {:.2%}'.format(nPLS, AccLDAmean))

	# Calculate accuracy and f1 score for k-folds cross validation
# 	AccLDAScores = cross_validation.cross_val_score(clfLDA, TrnsfrmPls[0], ExampleClasses, cv=cvLDA)
# 	#print(AccScores)
# 	### Identify the incorrect predictions from LDA and print them
# 	# 	print('')
# # 	if AccLDAmean <= 0.1 and AccLDAmean > 0.0:
# # 		Idx = np.where(ExampleClasses != predLDA)
# # 		print('Incorrect Example indices ({}):'.format(len(Idx[0])), np.array(Idx[0], dtype='int'))
# # 		print('Predicted:', [LblList[i] for i in predLDA[Idx[0]]])
# # 		print('True:', [LblList[i] for i in ExampleClasses[Idx[0]]])
# # 		print('')
# # 	print("k-folds accuracy: %0.2f (+/- %0.2f)" % (np.mean(AccLDAScores), np.std(AccLDAScores) * 2))
# # 	F1LDAScores = cross_validation.cross_val_score(clfLDA, TrnsfrmPls[0], ExampleClasses, scoring='f1_micro', cv=cvLDA)
# # 	print("k-folds F1 score: %0.2f (+/- %0.2f)" % (np.mean(F1LDAScores), np.std(F1LDAScores) * 2))
# # 	
	Idx = int((nPLS-minPLS)/intervalPLS)
	StoreAccLDA[Idx] = AccLDAmean
	StoreAccPLS[Idx] = msepPLS
	StorekAccPLS[Idx] = abs(np.mean(msePLSScores))

#Vip = pls1vip()

# Calculate variance explained by each PC
SS = np.dot(pls.y_loadings_**2,np.dot(pls.x_scores_.T,pls.x_scores_))
ExplVar = SS/np.sum(SS)
np.set_printoptions(precision=2, suppress=True)
print('Variance Explained for each component:')
print(ExplVar*100)

# Run visualization
vis()


# import ipdb; ipdb.set_trace() # Breakpoint	
