import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys

def Lay(self):
        self.dtool.addWidget(self.imgbut,row=0,col=0)
        self.dtool.addWidget(self.specbut,row=0,col=1)
        self.dtool.addWidget(self.rm1dPlot,row=0,col=2)
        self.dtool.addWidget(self.tbl,row=0,col=3)
        self.dtool.addWidget(self.ewBut,row=0,col=4)
        self.dtool.addWidget(self.addImgbut,row=0,col=5)
        self.dtool.addWidget(self.Dcolorbut,row=0,col=6)
        self.dtool.addWidget(self.EcolorBut,row=0,col=7)
        self.dtool.addWidget(self.contBut,row=0,col=8)
        self.dtool.addWidget(self.getFitsBut,row=0,col=9)
        self.dtool.addWidget(self.addReg,row=1,col=0)
        self.dtool.addWidget(self.remReg,row=1,col=1)
        self.dtool.addWidget(self.nonPEW,row=1,col=2)
        self.dtool.addWidget(self.lincent,row=1,col=3)
        self.dtool.addWidget(self.fluxToLum,row=1,col=4)
        self.dtool.addWidget(self.Math,row=1,col=5)
        self.dtool.addWidget(self.save,row=1,col=6)
        self.dtool.addWidget(self.aztrace,row=2,col=0)
        self.dtool.addWidget(self.azdensity,row=2,col=1)
        self.dtool.addWidget(self.azautocorr,row=2,col=2)
        self.dtool.addWidget(self.azenergy,row=2,col=3)
        self.dtool.addWidget(self.azforest,row=2,col=4)
        self.dtool.addWidget(self.azjoint,row=2,col=5)
        self.dtool.addWidget(self.azparallel,row=2,col=6)
        self.dtool.addWidget(self.azposterior,row=2,col=7)
        self.dtool.addWidget(self.artgauss,row=2,col=8)