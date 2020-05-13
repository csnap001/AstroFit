import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys





def doesNothing():
    pass

def Pos(self,Graphics_item):
    Graphics_item.scene().sigMouseClicked.connect(self.onClick)
    
    return 

def off(self,Graphics_item):
    Graphics_item.scene().sigMouseClicked.connect(doesNothing)#needed because button is pressed on startup
    Graphics_item.scene().sigMouseClicked.disconnect()


