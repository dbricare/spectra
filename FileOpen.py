"This module uses Qt to open a file dialog box."

def openFileDialog(DialogCaption=None, ReadFolder=None, FileFilter=None):
	from PySide import QtCore, QtGui
	import sys

	app = QtGui.QApplication.instance()
	if app is None:
		app = QtGui.QApplication(sys.argv)	
    
	FileList, _ = QtGui.QFileDialog.getOpenFileNames(None, \
	caption=DialogCaption, dir=ReadFolder, filter=FileFilter)
		
	if not FileList:
		raise ValueError('File list is empty.') 
		
	return (FileList)
