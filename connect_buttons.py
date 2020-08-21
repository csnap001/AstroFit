
import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys



def connect(self):
        self.imgbut.clicked.connect(self.file2d)
        self.specbut.clicked.connect(self.file1d)
        self.rm1dPlot.clicked.connect(self.removePlot)
        self.tbl.clicked.connect(self.fileTab)
        self.Dcolorbut.clicked.connect(self.set1dFluxCol)
        self.EcolorBut.clicked.connect(self.set1dErrCol)
        self.ewBut.clicked[bool].connect(self.fitting_1d)#TODO: change this
        self.contBut.clicked[bool].connect(self.Cont_Fit)
        self.getFitsBut.clicked.connect(self.showPDFS)
        self.addReg.clicked.connect(self.whichPlot)
        self.remReg.clicked.connect(self.removeRegPlot)
        return