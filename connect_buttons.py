
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
        self.ewBut.clicked.connect(self.fitting_1d)
        self.nonPEW.clicked.connect(self.nonParamEW)
        self.lincent.clicked.connect(self.LineCenter)
        self.contBut.clicked.connect(self.Cont_Fit)
        self.getFitsBut.clicked.connect(self.showPDFS)
        self.addReg.clicked.connect(self.whichPlot)
        self.remReg.clicked.connect(self.removeRegPlot)
        self.fluxToLum.clicked.connect(self.Flux_to_Lum)
        return