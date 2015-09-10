# Commands used to load in data, organize, perform LDA classification of 2D PCA of data

import itertools
import numpy as np
import matplotlib.pyplot as plt
import sklearn, itertools
from sklearn import lda, preprocessing, cross_validation
from sklearn.decomposition import PCA
from sklearn.metrics import f1_score, accuracy_score, pairwise
from FileOpen import openFileDialog
# from SignificanceTest import transtestsig, mahdist
# from scipy.stats import levene
from sklearn.cluster import KMeans

# Global variables that may need modification
PixelStart = 75 #154 #145 #220
PixelEnd = 1200 #930 #930 #1104
ClassSubDiv = 1
ExplVar = 0.99


#--------------------------------------------------------------------------------------
# Create label names

def genLblList(FileList):
	RemoveDir = [s.replace(ReadFolder,'') for s in FileList]
	Suffix = '-LNS.xls'
	LblList = [s.replace(Suffix,'') for s in RemoveDir]

	return(LblList)


#--------------------------------------------------------------------------------------
# Load files, remove first 5 columns and organize data into single array with corresponding example class in a vector

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
# Find unique items in list while preserving order

def uniquefun(items):
    found = set([])
    keep = []

    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)

    return(keep)


#--------------------------------------------------------------------------------------
# Compute and print metrics for clustering when ground truths are known

def clustrmetrics():

	from sklearn import metrics
	
	print('Homog','Complt','ARI','AMI',sep='\t')
	print('{:.3} \t  {:.3}  \t {:.3} \t {:.3}'.format(metrics.homogeneity_score(y,km.labels_), metrics.completeness_score(y,km.labels_), metrics.adjusted_rand_score(y,km.labels_), metrics.adjusted_mutual_info_score(y,km.labels_)))


#--------------------------------------------------------------------------------------
# Compute and print Euclidean separation between each pair of cluster centroids

def eucldistpairs():

	# For each combination of 2 clusters, compute the euclidean separation between cluster centroids in PCA space, doesn't make sense to compute mahalanobis distances in PCA space
	DistEuclid = pairwise.pairwise_distances(km.cluster_centers_, metric='euclidean')
	ClustrPairs = list(itertools.combinations(list(range(nClusters)),2))

	print(''*2)
	print('Euclidean distance between cluster centroids in 2D PCA:')

# Match centroid labels to data group-labels then print distances
	for i in range(len(ClustrPairs)):
		print('')
		Ctroid = [LblList[uniqueLbls[s]] for s in ClustrPairs[i]]
		print(Ctroid[0], 'vs.', Ctroid[1])
		print('Separation in {:.0%} variance PCA space:'.format(ExplVar), '{:.2e}'.format(DistEuclid[ClustrPairs[i]]))



#--------------------------------------------------------------------------------------
# Visualization

def vis():

	from matplotlib.patches import Ellipse
	import ConfEllipse

	plt.ion()
	plt.rc('font', family = 'Arial', size='16')

	# Custom colors and markers
	ColorOpts = ['RoyalBlue', 'FireBrick', 'OliveDrab', 'DarkCyan', 'Coral', 'CornflowerBlue', 'SaddleBrown', 'YellowGreen', 'DarkViolet', 'DarkTurquoise']
# 	ColorOpts = ['CornflowerBlue','YellowGreen','MediumPurple','Coral', 'DarkCyan',  'Crimson', 'Gold', 'BurlyWood', 'HotPink', 'SaddleBrown', 'Aqua', 'Lime']
	ColorSeq = ColorOpts[:int(len(LblList)/ClassSubDiv)]*ClassSubDiv
	MarkerOpts = ['x','o','s','^','D','x']
	#LgdLoc = 1+np.arange(4)
	#MarkerSeq = MarkerOpts[:ClassSubDiv-1]

	nClss = len(LblList)
	
	# For customized labels not derived from file names
	#LblList = ['Mutant','Parent']
	#LblList = ['Supernatant - no wash', 'Supernatant - wash 4x', 'Citrate']
	#PCAsort = np.zeros(np.shape(PCACoeff[:,:2]))
	#CustomClss = []
	#for i in range(nClss):
	#	CustomClss.extend([i]*int(len(Data[:,0])/nClss))
	#y = CustomClss



### Plot first two principal components
	fig = plt.figure(figsize=(12,7))
	ax = fig.add_subplot(111)
	j = 0
#	for i, c, class_name in zip(range(len(LblList)), ColorSeq, LblList):
#		ax.scatter(PCACoeff[np.where(np.array(y)==i),0], #PCACoeff[np.where(np.array(y)==i),1], c=c, label=class_name, s=80, marker=MarkerOpts[j])
#		if i >= ((j+1)*len(LblList)/ClassSubDiv)-1:
#			j+=1
	for i,c,class_name in zip(range(nClss),ColorSeq,LblList):
		PointIdx = np.where(np.array(y)==i)
		ax.scatter(PCACoeff[PointIdx,0], PCACoeff[PointIdx,1], c=c, label=class_name, s=100, marker=MarkerOpts[j], linewidth=2)   #, edgecolors = '#333333')
		if i >= ((j+1)*len(y)/ClassSubDiv)-1:
			j+=1
# 	ax.scatter(km.cluster_centers_[:,0], km.cluster_centers_[:,1], c='k', s=200, marker='x', linewidths=(4,))   # Centroid

	# Shrink axes box width by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put legend inside or to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1)
# 	plt.legend(loc='upper right', scatterpoints=3)

	plt.ylabel(
	'Principal Component 2 '+'({:.1%})'.format(pca.explained_variance_ratio_[1]))
	plt.xlabel(
	'Principal Component 1 '+'({:.1%})'.format(pca.explained_variance_ratio_[0]))
# 	plt.title('First two principal components')
	PCAmax = np.amax(np.absolute(PCACoeff[:,:2]))
	limLog = np.floor(np.log10(PCAmax))
	limInt = np.ceil(PCAmax/(10**np.floor(np.log10(PCAmax))))
	lim = limInt*(10**limLog)
	plt.xlim(-1*lim,lim)
	plt.ylim(-1*lim,lim)

	plt.subplots_adjust(bottom=0.10,left=0.10,right=0.75,top=0.95,wspace=0.1,hspace=0.1)



### Plot 95% statistical hypothesis testing limit for each centroid, if there is no overlap then clusters are significant for p<0.05 (Worley et al, Anal Biochem, 2013)
# 	for i in range(nClss):
# 		StdClustr = np.sqrt(np.var(TrnsfrmK[np.where(km.labels_==i),:2],axis=1)).reshape(2)
# 		ConfInt = Ellipse(xy=km.cluster_centers_[uniqueLbls[i],:], height=4*StdClustr[1], width=4*StdClustr[0], color=ColorOpts[i], alpha=0.2)
# 		ax.add_artist(ConfInt)
# 	for i in range(nClss):
# 		points = TrnsfrmK[km.labels_==i,:2]
# 		pos = km.cluster_centers_[uniqueLbls[i],:]
# 		cov = np.cov(points, rowvar=0)
# 		ConfEllipse.plot_cov_ellipse(cov, pos, nstd=3, color=ColorOpts[i], alpha=0.3)
	
	
	
### Plot loadings
	nPC = 2   # number of PCs to include on plot
	LoadgClr = ['Blue', 'Crimson', 'Green']
	plt.figure(figsize=(12,8))
	for i in range(nPC):
		plt.plot(WaveNumber, pca.components_[i,:], color=LoadgClr[i], linestyle='-',
		label='PC '+str(i+1)+' Loading')
	plt.xlim(WaveNumber[0],WaveNumber[PixelEnd-PixelStart-1])
	plt.legend(loc='lower left')
	#plt.title('PC Loadings')
# 	plt.ylabel('Loading Value (a.u.)')
	plt.xlabel(r'Wave Number ($\mathregular{cm}^{-1}$)')



### Plot first two discriminant axes from LDA
	fig = plt.figure(figsize=(12,8))
	ax = plt.subplot(111)
	j = 0
#	for i, c, class_name in zip(range(len(LblList)), ColorSeq, LblList):
#		P=plt.scatter(Xlda[np.where(np.array(y)==i), 0], Xlda[np.where(np.array(y)==i), #1], c=c, label=class_name, s=80, marker=MarkerOpts[j])
#		if i >= ((j+1)*len(LblList)/ClassSubDiv)-1:
#			j+=1
			
	for i,c,class_name in zip(range(nClss),ColorSeq,LblList):
		ax.scatter(Xlda[np.where(np.array(y)==i),0], Xlda[np.where(np.array(y)==i),1], c=c, label=class_name, s=100, marker=MarkerOpts[j])
		if i >= ((j+1)*len(y)/ClassSubDiv)-1:
			j+=1

	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1)

	plt.ylabel('Discriminant Axis 2')
	plt.xlabel('Discriminant Axis 1')
	plt.title('LDA with two discriminant axes')


#--------------------------------------------------------------------------------------
# Main


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


# Load wavenumber calibration and trim data window if necessary
WaveNumber = np.loadtxt('Pixel-Wavenumber-Grating600-866.1nm.xls',delimiter="\t",ndmin=2)
WaveNumber = WaveNumber[PixelStart:PixelEnd]
Data = Data[:,PixelStart:PixelEnd]


# Perform PCA to reduce dimensionality and select features responsible for variance between examples.
# For Raman data it is not necessary to standardize (mean=0,variance=1). PCA will mean center the normalized data. Also, PCA will not necessarily remove multicollinearity from the training features (i.e., make features independent).
#DataStd = preprocessing.StandardScaler().fit_transform(Data)
pca = PCA(n_components=None)
PCACoeff = pca.fit_transform(Data)   # Fit standard or non-standardized data
iVar = 0
CurrVar = 0


# Use enough PCs to describe original data given percentage of the variance
while CurrVar <= ExplVar:
	CurrVar = CurrVar + pca.explained_variance_ratio_[iVar]
	iVar+=1
print('')
print('Number of features in original data:', Data.shape[1])
print('Number of PCs explaining {:.0%} of variance:'.format(ExplVar), iVar)
X = PCACoeff[:,:iVar]
y = ExampleClasses


# Quantify 2D group separation and determine statistical significance, kMeans works well for the first few PCs (Yeung and Ruzzo, Bioinformatics, 2001)
# Compute groups centroids with k-means for first two principal components and transform
nClusters = len(LblList)
km = KMeans(n_clusters=nClusters)
TrnsfrmK = km.fit_transform(X[:,:2])


# Compute and print Euclidean separation between each pair of cluster centroids
uniqueLbls = uniquefun(km.labels_)
# eucldistpairs()

# Compute and print metrics for clustering when ground truths are known
clustrmetrics()

#Perform classification then transform PCA data to LDA basis
clf = lda.LDA(n_components=None)
Xlda = clf.fit_transform(X,ExampleClasses)


# Determine cross-validation scores using k-folds repeated n_iter times with a new random sorting
#cv = cross_validation.StratifiedShuffleSplit(ExampleClasses, n_iter=10, test_size=0.2, random_state=None)
#print('')
#AccScores = cross_validation.cross_val_score(clf, X, y, cv=cv)
#print(AccScores)
#print('')
#print("Accuracy: %0.2f (+/- %0.2f)" % (np.mean(AccScores), np.std(AccScores) * 2))
#F1Scores = cross_validation.cross_val_score(clf, X, y, scoring='f1_micro', cv=cv)
#print("F1 score: %0.2f (+/- %0.2f)" % (np.mean(F1Scores), np.std(F1Scores) * 2))

# Run visualization
vis()


#import ipdb; ipdb.set_trace() # Breakpoint	
