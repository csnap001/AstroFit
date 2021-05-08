
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
        self.addReg.clicked.connect(self.whichPlot)
        self.remReg.clicked.connect(self.removeRegPlot)
        self.fluxToLum.clicked.connect(self.Flux_to_Lum)
        self.Math.clicked.connect(self.Data_operations)
        self.aztrace.clicked.connect(self.arviz_trace)
        self.azdensity.clicked.connect(self.arviz_density)
        self.azautocorr.clicked.connect(self.arviz_autocorrelation)
        self.azenergy.clicked.connect(self.arviz_energy)
        self.azforest.clicked.connect(self.arviz_forest)
        self.azjoint.clicked.connect(self.arviz_joint)
        self.azparallel.clicked.connect(self.arviz_parallel)
        self.azposterior.clicked.connect(self.arviz_plot_posterior)
        self.save.clicked.connect(self.save_data)
        self.artgauss.clicked.connect(self.artificial_gauss)
        self.guessZ.clicked.connect(self.Show_lines)
        self.dLines.clicked.connect(self.Remove_lines)
        self.r2D.clicked.connect(self.remove2D)
        self.sm.clicked.connect(self.start_smoothing)
        self.pypeit.clicked.connect(self.openFileNamesDialog)
        return