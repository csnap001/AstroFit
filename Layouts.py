import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys

def Lay(self):
        #laying out buttons and image spaces
        self.mainLay = qt.QGridLayout(self)# creates 2x2 grid layout of widgets
        self.pltimgLay = qt.QGridLayout(self)
        self.specButLay1d = qt.QVBoxLayout(self)
        self.specButLay2d = qt.QVBoxLayout(self)
        self.butLay = qt.QHBoxLayout(self)
        self.butLay.addWidget(self.imgbut)#these show up in this order
        self.butLay.addWidget(self.tbl)
        self.butLay.addWidget(self.specbut)

        #setting upmost level layout
        self.mainLay.addLayout(self.butLay,0,0,1,1) #can extend the 2x2 structure of gridlayout this way
        self.mainLay.addWidget(self.img_view,0,1,2,1)
        self.mainLay.addLayout(self.pltimgLay,1,0,1,1)

        #setting embedded layout
        self.pltimgLay.addWidget(self.Gwin1d,0,0,1,1)
        self.pltimgLay.addWidget(self.tabWid,1,0,1,1)
        self.pltimgLay.addLayout(self.specButLay1d,0,1,1,1)
        self.pltimgLay.addLayout(self.specButLay2d,1,1,1,1)

        #setting 1d embedded layout
        self.specButLay1d.addWidget(self.ewBut)
        self.specButLay1d.addWidget(self.Dcolorbut)
        self.specButLay1d.addWidget(self.EcolorBut)
        self.specButLay1d.addWidget(self.contBut)
        self.specButLay1d.addWidget(self.getFitsBut)
        self.specButLay1d.addWidget(self.fitProgress)
        

        #setting 2d embedded layout
        self.specButLay2d.addWidget(self.addImgbut)