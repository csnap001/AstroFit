import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys

def Lay(self):
        self.toolbar.addWidget(self.imgbut)
        self.toolbar.addWidget(self.specbut)
        self.toolbar.addWidget(self.rm1dPlot)
        self.toolbar.addWidget(self.addReg)
        self.toolbar.addWidget(self.remReg)
        self.toolbar.addWidget(self.tbl)
        self.toolbar.addWidget(self.ewBut)
        self.toolbar.addWidget(self.nonPEW)
        self.toolbar.addWidget(self.lincent)
        self.toolbar.addWidget(self.addImgbut)
        self.toolbar.addWidget(self.Dcolorbut)
        self.toolbar.addWidget(self.EcolorBut)
        self.toolbar.addWidget(self.contBut)
        self.toolbar.addWidget(self.fluxToLum)
        self.toolbar.addWidget(self.Math)
        self.toolbar.addWidget(self.aztrace)
        self.toolbar.addWidget(self.azdensity)
        self.toolbar.addWidget(self.azautocorr)
        self.toolbar.addWidget(self.azenergy)
        self.toolbar.addWidget(self.azforest)
        self.toolbar.addWidget(self.azjoint)
        self.toolbar.addWidget(self.azparallel)
        self.toolbar.addWidget(self.azposterior)
        self.toolbar.addWidget(self.save)
        self.toolbar.addWidget(self.artgauss)
        self.toolbar.addWidget(self.guessZ)
        self.toolbar.addWidget(self.dLines)
        self.toolbar.addWidget(self.r2D)
        self.toolbar.addWidget(self.sm)
        self.toolbar.addWidget(self.pypeit)
        self.toolbar.addWidget(self.coadd1d)
        self.toolbar.addWidget(self.s_coadd)