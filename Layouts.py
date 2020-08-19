import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys

def Lay(self):
        self.dtool.addWidget(self.imgbut)
        self.dtool.addWidget(self.specbut)
        self.dtool.addWidget(self.rm1dPlot)
        self.dtool.addWidget(self.tbl)
        self.dtool.addWidget(self.ewBut)
        self.dtool.addWidget(self.addImgbut)
        self.dtool.addWidget(self.Dcolorbut)
        self.dtool.addWidget(self.EcolorBut)
        self.dtool.addWidget(self.contBut)
        self.dtool.addWidget(self.getFitsBut)