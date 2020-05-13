import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys



def Views(self):
        self.Gwin1d = pg.GraphicsWindow() #Better to work with a graphics window
        self.text = pg.LabelItem(justify='right')
        self.Gwin1d.addItem(self.text) #Must add Label BEFORE plot, otherwise get unstable behavior
        self.plot_view = self.Gwin1d.addPlot(row=1,col=0) #row and col might be unnecessary, but may as well leave it for now
        self.img_view = pg.GraphicsLayoutWidget()
        self.tabWid = qt.QTableWidget()
        self.plot_view.setCursor(core.Qt.CursorShape(2)) #sets the cursor for this widget
        self.img_view.setCursor(core.Qt.CursorShape(2))
        