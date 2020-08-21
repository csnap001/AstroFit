

import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys

def buttons(self):
        self.imgbut = qt.QPushButton("Show 2d fits",self)
        self.specbut = qt.QPushButton("Add 1d Plot",self)
        self.rm1dPlot = qt.QPushButton("Remove 1d Plot",self)
        self.addReg = qt.QPushButton("Add region",self)
        self.remReg = qt.QPushButton("Remove region",self)
        self.tbl = qt.QPushButton("Fits table",self)
        self.ewBut = qt.QPushButton("EW",self)
        self.ewBut.setCheckable(True)
        self.addImgbut = qt.QPushButton("Add Images",self)
        self.Dcolorbut = qt.QPushButton("Data color",self)
        self.EcolorBut = qt.QPushButton("Err color", self)
        self.contBut = qt.QPushButton("Fit Cont", self)
        self.getFitsBut = qt.QPushButton("Get Fits",self)
        return