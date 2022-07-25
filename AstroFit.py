"""
Created Thursday 8th, 2019 by Christopher Snapp-Kolas

#TODO: possible overhaul of the structure of the code. It seems like it would be better to set up "dock layers" where
each layer is performing a specific task (i.e. dock 1-3 each are working on 1d spectra and dock 4 has a table open)
Then for each layer displaying on the left-hand(or right?) side displaying the data files being used, fits performed
and any additional elements relevant to the display. I think it would be best to create in place a folder that "saves"
all this information from which the user can access. These files should be deleted upon shut-down of the layer.
This list or lists would include the pickle fits, corner plots, and original fits files. Perhaps place user adjustable
slots on the righthand side (i.e. What parameterization, how many iterations of MCMC, number of masking regions, etc.).
buttons, or other operations, relevant to that kind of object (i.e. 1d spectra have differernt functions than 2d)
can perhaps go above or below the plot.

#TODO:break into sections. 1d spectra (and functions), 2d spectra (and functions), and table (with plotting/data manipulation
# functions)


TODO: add "clear screen" option for clearing displays
TODO: Add threading of functions to allow continued functionality of other QtGui objects. Although there
should be checks to ensure the user isn't overworking their computer.
TODO: Add progress bar for fits: self.progress = QtGui.ProgressBar(self) -> self.progress.setValue(**some increasing number**)
This is a difficult task, as pymc3 holds all this info internally and doesn't appear to allow for access
TODO: Re-organize storage of data. There should be meta-data fxns that feed to data analysis fxns
Should write out and organize my thoughts on this. 
TODO: Should save data to files as well (.dat?), currently using pickle, works for me but not generally what people
will want
TODO: Add operations to data elements (mostly to account for inverse variance stuff)

TODO: Add in Stellar population synthesis model plotting from Hoki (which is a wrapper to BPASS)
This likely needs more work with data containers to help out with the workflow of the code

TODO: make plots histograms? only need to set stepMode="center" in plot call, but requires
len(x) = len(y) + 1. might be able to acomplish this using numpy hist

TODO: create capability of multiplot over single axes.

TODO: create file select that saves selected data products to PDF

TODO: create capability to guess redshift and plot known emission/absorption lines
perhaps create using pg.infitelines(...movable=True), then force them to all be
at the same redshift (would be nice to link positions, needs helper update function to be done)

TODO: Plot parameter space (image_view?) and allow user to move through space and choose points
then plot using those parameters over original data set. This is good for visualization.
This will make clear whether a given fit behaved as expected.

TODO: Add in capability to choose what fitting method to use (i.e. NUTS, Hamiltonian, MCMC)
as these pre-exist in the pymc3 framework. (some consideration should be given to upgrading
to pymc4)

TODO: add in slider for smoothing plots (need to hold original data set to be able to reset)
Created Slider class to handle this. Connect Slider widget to update plot with gaussian smoothing
function. example: Slide = Slider(lower_bound,Upper_bound), Slide.valueChanged.connect(smooth_plot)

TODO: Add in ability to remove various added plotItems (i.e. fits and whatnot), there exists a 
removeItem(item) function. Would be nice to hold these items in memory and allow the user to add
and remove at will.

TODO: Grab names from table (all keywords) and allow user to select and then move view to there
and highlight. Would be nice to make operations on table elements possible. Both the ability to
edit in place, and the ability to change entire rows or columns, or perform operations between 
rows and columns.

Module for GUI spectroscopic fitting environment based on pymc3, pyqtgraph,
and astropy. (Possibly, desired) This module will also have basic image arithmatic capabilities.

Meta data:

Class Functions:

"""

from re import L
import PyQt5.QtWidgets as qt
from astropy.io.fits.hdu.image import ImageHDU
from astropy.io.fits.hdu.table import BinTableHDU, TableHDU
from pyqtgraph.Qt import QtCore,QtGui
import pyqtgraph as pg
import pandas as pd
from astropy.io import fits,ascii
from astropy.table import Table
from astropy import units as u
from astropy.modeling import models
from pyqtgraph.dockarea import *
import multiprocessing as multi
import numpy as np
import sys

import create_buttons as cb
import connect_buttons as conb
from Layouts import Lay
import time
import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
from functools import partial
import os
import scipy as sc
from scipy.optimize import curve_fit
from IPython import embed
from astropy import cosmology
from astropy.modeling.functional_models import Voigt1D
import pickle
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, \
    QVBoxLayout, QWidget, QToolBar
from scipy.ndimage import gaussian_filter1d
from spectres import spectres
import theano.tensor as t
import theano
from multipledispatch import dispatch

def conf_low(data,half=0.34):
    mean = np.mean(data)
    sd = np.std(data)
    #embed()
    kde = sc.stats.gaussian_kde(data,0.05)
    i = 0
    conf = 0
    while conf < half:
        i += 1
        conf = kde.integrate_box_1d(mean-i*0.01*sd,mean)
        if i > 500:
            #print("timeout: distribution has problems")
            break
    #print("finished low")
    return mean - i*0.01*sd

def conf_high(data,half=0.34):
    mean = np.mean(data)
    sd = np.std(data)
    #embed()
    kde = sc.stats.gaussian_kde(data,0.05)
    i = 0
    conf = 0
    while conf < half:
        i += 1
        conf = kde.integrate_box_1d(mean,mean+i*0.01*sd)
        if i > 500:
            #print("timeout: distribution has problems")
            break
    #print('finished high')
    return mean + i*0.01*sd

def Powpgauss(x,A0,x0,b0,f0,a0,cent):
    return f0*(x/cent)**a0 + abs(A0)*np.exp(-((x-x0)/(np.sqrt(2)*b0))**2)

def Pow(x,f0,a0,cent):
    return f0*(x/cent)**a0

def gauss0(x,A0,x0,b0):
    return abs(A0)*np.exp(-(((x-x0)/b0)**2.))

def polygauss(x,Ag,xg,bg,b,m):
    return polynomial(x,b,m) + gauss0(x,Ag,xg,bg)

def piece(x,b,xb,xc,xg,A0,Ag,a,bg,m):
    return pm.math.switch(x <= xb,polygauss(x,Ag,xg,bg,b,m),Pow(x,A0,a,xc)).eval()

@theano.compile.ops.as_op(itypes=[t.dvector, t.dscalar,t.dscalar, t.dscalar,t.dscalar, t.dscalar,t.dscalar, t.dscalar,t.dscalar, t.dscalar],otypes=[t.dvector])
@dispatch(object,object,object,object,object,object,object,object,object,object)
def piecewise(x,xb,xc,xg,A0,Ag,a,bg,b,m):
    #y = t.zeros(x.shape)
    #params = t.switch(x <= xb,[Ag,xg,bg,b,m],[A0,xc,a])
    #func = t.switch(x <= xb,polygauss,Pow)
    #y = t.switch(x <= xb,polynomial(x,b,m) + gauss0(x,Ag,xg,bg),Pow(x,A0,xc,a))
    return pm.math.switch(x <= xb,polygauss(x,Ag,xg,bg,b,m),Pow(x,A0,a,xc))

@dispatch(np.ndarray,float,float,float,float,float,float,float,float,float)
def piecewise(x,xb,xc,xg,A0,Ag,a,bg,b,m):
    print("numpy array")
    y = np.zeros_like(x)
    mask = x <= xb
    y[mask] = polygauss(x[mask],Ag,xg,bg,b,m)
    y[~mask] = Pow(x[~mask],A0,a,xc)
    return y

@dispatch(float,float,float,float,float,float,float,float,float,float)
def piecewise(x,xb,xc,xg,A0,Ag,a,bg,b,m):
    print("All floats")
    y = np.zeros_like(x)
    mask = x <= xb
    y[mask] = polygauss(x[mask],Ag,xg,bg,b,m)
    y[~mask] = Pow(x[~mask],A0,xc,a)
    return y

def Voigt(x,A0,Nl,x0,sig):
    return abs(A0)*np.exp(-0.7580*(Nl/10e13)*(10/sig)*gauss0(3e5*(x-x0)/x0,1,0,sig))

def lorentzian(x,x0,a,gam):
    return abs(a)*gam**2/(gam**2 + (x-x0)**2)


def linear(x,a,b):
    return a*x + b

def polynomial(x,*a):
    return sum([p*(x**i) for i,p in enumerate(a)])

class Slider(QWidget):
    def __init__(self, minimum, maximum, parent=None):
        super(Slider, self).__init__(parent=parent)
        self.verticalLayout = QVBoxLayout(self)
        self.label = QLabel(self)
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout = QHBoxLayout()
        spacerItem = QSpacerItem(0, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.slider = QSlider(self)
        self.slider.setOrientation(Qt.Vertical)
        self.horizontalLayout.addWidget(self.slider)
        spacerItem1 = QSpacerItem(0, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.resize(self.sizeHint())

        self.minimum = minimum
        self.maximum = maximum
        self.slider.valueChanged.connect(self.setLabelValue)
        self.x = None
        self.setLabelValue(self.slider.value())

    def setLabelValue(self, value):
        self.x = self.minimum + (float(value) / (self.slider.maximum() - self.slider.minimum())) * (
        self.maximum - self.minimum)
        self.label.setText("{0:.4g}".format(self.x))

class App(QtGui.QMainWindow):

    def __init__(self):
        super().__init__()
        
        self.title = "AstroFit"
        self.left = 20
        self.top = 20
        self.width = 1000
        self.height = 1000

        self.plt = []#list containing plot information for each added 1d plot
        self.lines = []#list containing infinite lines to allow for removal
        self.imv = []#list containing images
        self.names1d = []#list containing names1d of 1d plots to allow for removal
        self.names2d = []#list containing names of 2d plots to allow for removal
        self.lrs = []#list of linear regions for data selection (i.e. for fitting continua or EWs)
        self.regPlot = []#list to hold region plots
        self.yrange = [0,1]
        self.prior = []#list of continuum priors that can be used
        self.fill = None
        self.memory = [] #used to have access to unsmoothed data

        self.is1d = False
        self.is2d = False
        self.isTab = False
        self.is1dPos = False
        self.is2dPos = False

        self.flux = []
        self.err = []
        self.initUI()

        #must have data atributes that can be grabbed by the user upon request
        self.contfit = [[],[]] #A list of lists holding the continuum fit parameters and names for the fit
        self.conterr = []
        self.pdfs = {}
        self.cname = ''
        #TODO: the output to user should show the model and the parameter names
        self.ewfit = [] #gives equivalent width fitting values
        self.ewferr = []
        self.ew = [] #gives equivalent width
        self.ewerr = []
        self.ewBounds = None
        #TODO: should save to a file placed wherever the user desires
        #TODO: other fits? 2d likely, but this will be done later. what about 1d?

        #NOTE: self.arviz will hold the data used in arviz visualizations. also for saving data
        self.arviz = {}

        self.toolbar = QToolBar("Spec Tools")
        self.addToolBar(self.toolbar) 
        #Other tools
        self.sampler = qt.QComboBox(self)
        self.sampler.addItems([ 'Metropolis', 'NUTS', 'Slice', 'HamiltonianMC', 'BinaryMetropolis'])
        #Creating buttons
        cb.buttons(self)
        #Connecting buttons to functions
        conb.connect(self)
        #Placing buttons in toolbar
        Lay(self)


    def initUI(self):
        
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.resize(self.width,self.height)
        self.setWindowTitle(self.title)
        self.dtool = Dock("Toolbar",size=(1,1))
        self.dplot = Dock("plots",size=(500,300))
        self.dTable = Dock("Table",size=(500,300))
        self.regDock = Dock("Regions",size=(500,300))
        self.plot2d = Dock("2d image",size=(300,300))
        # TODO: widgets also have "close" function, or at least imageviews do
        self.regWin = pg.GraphicsWindow()
        self.area.addDock(self.dtool,'top')
        self.area.addDock(self.dplot,'bottom',self.dtool)
        self.area.addDock(self.dTable,'bottom',self.dplot)
        self.area.addDock(self.regDock,'below',self.dTable)
        self.area.addDock(self.plot2d,"right",self.dplot)
        self.Gwin1d = pg.GraphicsWindow()
        self.Gwin1d.resize(1000,600)
        self.table = pg.TableWidget()
        self.dplot.addWidget(self.Gwin1d)
        self.regDock.addWidget(self.regWin)
        self.dTable.addWidget(self.table)
        self.dtool.hideTitleBar()
        self.Gwin1d.setBackground('w')
        self.slide_smooth = Slider(0.1,10)

        self.show()

    def imageHoverEvent(self,event):
        """Show the position, pixel, and value under the mouse cursor.
        """
        '''
        if event.isExit():
            self.plot2d.setTitle("Not working")
            return
        '''
        #embed()
        point = event[0]
        data = self.imv.image
        vb = self.imv.getView().getViewBox()
        #embed()
        point = vb.mapSceneToView(event[0])
        i, j = point.y(), point.x()
        #embed()
        i = int(np.clip(i, 0, data.shape[0] - 1))
        j = int(np.clip(j, 0, data.shape[1] - 1))
        val = data[j, i]
        wl = self.imv.tVals[j,i]
        ppos = self.imv.mapToParent(point.toPoint())
        x, y = ppos.x(), ppos.y()
        self.plot2d.setTitle(r"pixel: (%d, %d)  value: %g  $wl (\AA)$: (%0.1f)" % (i, j, val,wl))

    def smooth_update(self):
        data = self.memory[0]
        wl = data[0].getData()[0] 
        flux = data[0].getData()[1]
        err = data[1].getData()[1]
        sig = self.slide_smooth.x
        sm_flux = gaussian_filter1d(flux,sig)
        self.plt[0].clear()
        self.plt[0].plot(wl,sm_flux,pen='b',stepMode=True)
        self.plt[0].plot(wl,err,pen='k',stepMode=True)
        if len(self.lrs) > 0:
            self.lrs[0].sigRegionChanged #TODO: doesn't resolve lost region problem
            self.updateLRplot()
            self.updateLR()

    def start_smoothing(self):
        if len(self.plt) > 0:
            self.dplot.addWidget(self.slide_smooth,row=0,col=1)
            self.slide_smooth.slider.valueChanged.connect(self.smooth_update)
            self.memory.append(self.plt[0].listDataItems())
    #TODO: generalize to multiplots

    def updateLR(self):
        #TODO: The y-range for each plot is only connected to the most recent region, why?
        for i in range(len(self.lrs)):
            if self.lrs[i].sigRegionChanged:
                xlen = self.lrs[i].getRegion()
                self.lrs[i].setRegion(xlen)
                abdata = self.regPlot[i].listDataItems()
                data = (abdata[0].getData()[0],abdata[0].getData()[1],abdata[1].getData()[1])
                mask = (data[0][:-1] > self.regPlot[i].getViewBox().viewRange()[0][0]) & (data[0][:-1] < self.regPlot[i].getViewBox().viewRange()[0][1])
                #mask was of different length then data, therefore had to adjust using slice
                if len(data[1][mask]) > 0: self.yrange = [np.min(data[1][mask]) - 0.8*np.min(data[1][mask]),1.5*np.max(data[1][mask])]

    def addRegionPlot(self,data):
        num = len(self.regPlot)
        if num < 2:
            self.regPlot.append(self.regWin.addPlot(title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux',stepMode=True)
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error',stepMode=True)
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)
        elif num >= 2 and num < 4:
            self.regPlot.append(self.regWin.addPlot(row=1,col=num-2,title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux',stepMode=True)
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error',stepMode=True)
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)
        elif num >= 4 and num < 6:
            self.regPlot.append(self.regWin.addPlot(row=2,col=num-4,title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux',stepMode=True)
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error',stepMode=True)
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)

    def addPlot(self,name):
        #TODO:need to figure out how to handle placement after removing plots
        #Likely better to change this to tabs instead of multiple plots in a 
        #single window
        num = len(self.plt)
        if num < 2:
            self.plt.append(self.Gwin1d.addPlot(title=name))
            n = len(self.plt)
            self.plt[n-1].setCursor(QtCore.Qt.CrossCursor)
        elif num >= 2 and num < 4:
            self.plt.append(self.Gwin1d.addPlot(row=1,col=num-2,title=name))
            n = len(self.plt)
            self.plt[n-1].setCursor(QtCore.Qt.CrossCursor)
        elif num >=4 and num < 6:
            self.plt.append(self.Gwin1d.addPlot(row=2,col=num-4,title=name))
            n = len(self.plt)
            self.plt[n-1].setCursor(QtCore.Qt.CrossCursor)
        else:
            qt.QMessageBox.about(self,"Too Many!","Too many plots, not adding")
            return -1

    def removePlot(self):
        choice, ok = qt.QInputDialog.getItem(self,"Removing 1d?","Which plot?:",self.names1d,0,False)
        if ok:
            loc = self.names1d.index(choice)
            self.Gwin1d.removeItem(self.plt[loc])
            self.plt.pop(loc)
            self.names1d.pop(loc)
            self.flux.pop(loc)
            self.err.pop(loc)
        else:
            qt.QMessageBox.about(self,"Done","Not Removing any plots")
   

    def removeRegPlot(self):
        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        choice, ok = qt.QInputDialog.getItem(self,"Removing Region","Which Region?:",np.array(check,dtype=str),0,False)
        ind = check.index(float(choice))
        if ok:
            self.regWin.removeItem(self.regPlot[ind])
            self.regPlot.pop(ind)
            self.plt[0].removeItem(self.lrs[ind])#TODO:need to figure out how to grab parent plot
            self.lrs.pop(ind)#The two should increase and be placed together in these lists
            #TODO: need to double check this works well
        else:
            qt.QMessageBox.about(self,"Done","Not Removing Region")

    def mouseMoveEvent(self,event):
        pos = event
        for plot in self.plt:
            vb = plot.getViewBox()
            if isinstance(event,QtCore.QPointF) and plot.sceneBoundingRect().contains(pos):#checks if variable is an instance of the given class
                point = vb.mapSceneToView(event)
                plot.setLabel("top",text='x={0:1f}, y={1:.2E}'.format(point.x(),point.y()))

    #grabs .fits file for 1d visualization
    def file1d(self):
        self.is1d = True
        self.openFileNameDialog(is1d=self.is1d)

    #grabs .fits file for 2d visualization (allows for multi-extension frames)
    def file2d(self):
        self.is2d = True
        self.openFileNameDialog(is2d=self.is2d)

    def fileflux(self):
        self.openFileNamesDialog(isfluxest=True)

    def set1dFluxCol(self):
        isFlux = True
        self.plotColoring(isFlux=isFlux)

    def set1dErrCol(self):
        isErr = True
        self.plotColoring(isErr=isErr)

    def fileTab(self):
        self.isTab = True
        self.openFileNameDialog(isTab = self.isTab)

    def coadd(self):
        self.openFileNamesDialog(iscoadd=True)
    
    def open_pyp(self):
        self.openFileNamesDialog(ispypeit=True)
    #clear
    def openFileNameDialog(self,is1d = False ,is2d = False,isTab = False):
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        fileName, _ = qt.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","","All Files (*);;Fits Files (*.fits);;Python files (*.py)",options = options)
        _,exten = os.path.splitext(fileName) #Used to place constraints on filetype, exten grabs file extension

        if exten == ".fits":
            if is1d:
                strs = fileName.split('/')
                check = self.addPlot(name=strs[len(strs)-1])
                if check != -1:
                    self.plot_1d_fits(fileName=fileName)
                    self.names1d.append(strs[len(strs) - 1])
            if is2d:
                self.show_2d_fits(fileName=fileName)
            if isTab:
                self.table_create(fileName=fileName)
        elif exten == ".txt":
            strs = fileName.split('/')
            check = self.addPlot(name=strs[len(strs)-1])
            if check != -1:
                self.plot_1d_txt(fileName=fileName)
                self.names1d.append(strs[len(strs) - 1])
        else:
            qt.QMessageBox.about(self,"Opening File","ERROR: (Program Name) only supports .fits files")
            #TODO: Need to come up with name for GUI application

    def openFileNamesDialog(self,iscoadd=False,ispypeit=False,isfluxest=False):
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        files, _ = qt.QFileDialog.getOpenFileNames(self, "QFileDialog.getOpenFileNames()", "","All Files (*);;Fits Files (*.fits);;Python files (*.py)", options=options)
        if files and ispypeit:
            self.pyp_coadd(files)
        if files and iscoadd:
            #embed()
            self.coadd_1d(files)
        if files and isfluxest:
            self.estimate_flux(files)

    def estimate_flux(self,files):
        ID = []
        specFlux = []
        photFlux = []
        ratio = []
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        #Used to grab catalog of photometry to compare with
        catName, _ = qt.QFileDialog.getOpenFileName(self,"Photometric Catalog","","All Files (*);;Fits Files (*.fits);;Python files (*.py)",options = options)
        #TODO: Be good to generalize this to other filters
        if catName.split('.')[-1] == 'fits': 
            cat = fits.open(catName)
            T = Table(cat[1].data)
        else: 
            cat = ascii.read(catName)
            T = Table(cat)
        for f in files:
            name = str(f.split('/')[-1].split('_')[0])
            arr = np.where(T['Names'] == int(name))
            check = T['Cluster'][arr[0]] == 'A1689'
            if arr[0].size == 0:
                #TODO: probably should have a note here
                continue
            
            elif arr[0].size > 1:
                corr = np.where(f.split('/')[-1].split('_')[1] == T['Cluster'][arr[0]])
                check = T['Cluster'][arr[0]][corr[0]] == 'A1689'
            
            if check:
                HST_dat = ascii.read('DATA/hst_wfc_475W.txt')
                hst_wl = HST_dat['col1']
                hst_tp = HST_dat['col2']
            else:
                HST_dat = ascii.read('DATA/hst_wfc_435W.txt')
                hst_wl = HST_dat['col1']
                hst_tp = HST_dat['col2']
            ID.append(name)
            data = fits.open(f)[1].data
            wl = data['wave']
            flux = data['flux']
            newflux = spectres(hst_wl,wl,flux,fill=0,verbose=False)
            photFlux.append(T['F475W/F435W'][arr[0]][0])
            if check: specFlux.append(1e-17*(4750**2/(3*1e18))*np.trapz(newflux*hst_tp,x=hst_wl)/np.trapz(hst_tp,x=hst_wl))
            else: specFlux.append(1e-17*(4350**2/(3*1e18))*np.trapz(newflux*hst_tp,x=hst_wl)/np.trapz(hst_tp,x=hst_wl))
            ratio.append(photFlux[-1]/specFlux[-1])
        newT = Table([ID,specFlux,photFlux,ratio],names=('ID','specFlux','photFlux','ratio'))
        newT.write('Slitloss.fits',format='fits')           
    
    def table_create(self,fileName=""):
        #TODO: should add in capabilities such as table modification
        # that can be written to a new fits file. Would also be nice
        # if users could do calculations directly from the table.
        self.table.clear()
        fits = fileName.find(".fits")
        dat = fileName.find(".dat")
        txt = fileName.find(".txt")
        if not(fits == -1):
            data = Table.read(fileName,format='fits')
            self.table.setColumnCount(len(data.colnames))
            self.table.setRowCount(len(data))
            self.table.setHorizontalHeaderLabels(data.colnames)
            for i in range(len(data)):
                self.table.setRow(i,data[i])
        if not(dat == -1):
            pass
        if not(txt == -1):
            pass

    #TODO: Consider Plotting BINNED data, can use stepMode = True, but requires len(x) == len(y) + 1
    def updatePlot(self,count,wl,flux,err):
        self.plt[count].clear()
        pen = pg.mkPen(color='b')
        self.plt[count].addLegend()
        wl = np.append(wl, np.mean(np.diff(wl)) + wl[-1])
        embed()
        self.flux.append(self.plt[count].plot(wl,flux,pen=pen,name='Flux',stepMode=True))
        self.err.append(self.plt[count].plot(wl,err,name='Error',stepMode=True))
    
    def headerDisplay(self,hdul):
        hdus = []
        for item in hdul:
            key = item.header.keys()
            strng = ""
            while True:
                try: k = key.send(None)
                except StopIteration: break
                strng += k + "= " + str(item.header[k]) + "\n"
            hdus.append(strng)
            qt.QMessageBox.about(self,"Showing {}".format(item),strng)

    def plot_1d_txt(self,fileName=""):
        count = len(self.plt) - 1
        if fileName:
            data = ascii.read(fileName)
            while True:
                x, ok = qt.QInputDialog.getItem(self,"data","Choose your x axis:",data.keys(),0,False)
                y, ok = qt.QInputDialog.getItem(self,"data","Choose your y axis:",data.keys(),0,False)
                err, ok = qt.QInputDialog.getItem(self,"data","Choose your err axis:",data.keys(),0,False)
                data[x][np.isnan(data[x])] = 0e40
                data[y][np.isnan(data[y])] = 0e40
                data[err][np.isnan(data[err])] = 0e40
                x = data[x]
                y = data[y]
                err = data[err]
                if (x is not None) and (y is not None):
                    self.updatePlot(count,x,y,err)
                else: qt.QMessageBox.about(self,"Showing {0},{1}".format(x,y),"One data set is not valid")
                Happy, ok = qt.QInputDialog.getItem(self,"Good","Happy?:",["True","False"],0,False)
                if Happy == "True":
                    break

    #TODO: add capability to move plots along y-axis and to have multiple 
    #plots on same set of axes. Need to make sure calcs are on base data
    #and not the manipulated set

    def plot_1d_fits(self,fileName=""):
        count = len(self.plt) - 1
        #TODO: add legend (should be done for fits as well). legend = pg.LegendItem, legend.setParentItem(self.plot_view)
        data = []
        if fileName:
            hdul = fits.open(fileName)
            for i in range(len(hdul)):
                try:
                    data.append(hdul[i].data)
                except RuntimeError:
                    print(i)
            #NOTE: can generate header using key = hdul[i].header.keys(), then key.send(None) generates the next key
            # can have loop: str += key.send(None) + '= ' + hdul[i].header[key.send(None)] + '\n'
            self.headerDisplay(hdul)
            #embed()
            while True:
                E, ok = qt.QInputDialog.getItem(self,"Which spectrum?","choose extension",[str(i) for i in np.arange(len(data))],0,False)
                E = int(E)
                x, ok = qt.QInputDialog.getItem(self,"data","Choose your x axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
                y, ok = qt.QInputDialog.getItem(self,"data","Choose your y axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
                err, ok = qt.QInputDialog.getItem(self,"data","Choose your err axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
                data[E][x][np.isnan(data[E][x])] = 0
                data[E][y][np.isnan(data[E][y])] = 0
                data[E][err][np.isnan(data[E][err])] = 0
                x = data[E][x]
                y = data[E][y]
                err = data[E][err]
                if (x is not None) and (y is not None): self.updatePlot(count,x,y,err)
                else: qt.QMessageBox.about(self,"Showing {0},{1}".format(x,y),"One data set is not valid")
                Happy, ok = qt.QInputDialog.getItem(self,"Good","Happy?:",["True","False"],0,False)
                if Happy == "True":
                    break
            #TODO: How to grab data arbitrarily?


            pg.SignalProxy(self.plt[count].scene().sigMouseMoved, rateLimit=60,slot=self.mouseMoveEvent)
            self.plt[count].scene().sigMouseMoved.connect(self.mouseMoveEvent)

    def Data_operations(self):
        """
        The goal of this function is to allow the user to perform operations
        on the data elements that are being plotted
        """
        choice, ok = qt.QInputDialog.getItem(self,"Which plot","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not adjusting","Chose no data, not performing any operations.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        data_names = {"wl":data[0].getData()[0],"flux":data[0].getData()[1],"err":data[1].getData()[1]}
        #data_names['wl'] = np.append(data_names['wl'],np.mean(np.diff(data_names['wl']))+data_names['wl'][-1])
        while True:
            name, ok = qt.QInputDialog.getItem(self,"data","Which data element?:",data_names.keys(),0,False)
            if not(ok):
                qt.QMessageBox.about(self,"No operation?","No data element was chosen")
            op, ok = qt.QInputDialog.getItem(self,"what Math?","What operation would you like to perform",["invert","square root","scalar multiplication"],0,False)
            if op == "invert":
                data_names[name] = 1/data_names[name]
            if op == "square root":
                data_names[name] = np.sqrt(data_names[name])
            if op == "scalar multiplication":
                scalar,_ = qt.QInputDialog.getDouble(self,"Scalar Multiplication","What value?",0,-1e200,1e200,10)
                data_names[name] = data_names[name]*scalar

            self.plt[dat_choice].clear()
            self.flux[dat_choice] = self.plt[dat_choice].plot(data_names['wl'],data_names['flux'],pen='b',stepMode=True)
            self.err[dat_choice] = self.plt[dat_choice].plot(data_names['wl'],data_names['err'],stepMode=True)

            Happy, ok = qt.QInputDialog.getItem(self,"Good","Finished?:",["True","False"],0,False)
            if Happy == "True":
                break

    def updateLRplot(self):
        for i in range(len(self.regPlot)):
            if self.lrs[i].sigRegionChanged:
                self.regPlot[i].setXRange(*self.lrs[i].getRegion(),padding=0)
                self.regPlot[i].setYRange(self.yrange[0],self.yrange[1],padding=0)
                plots = self.regPlot[i].listDataItems()
                plots[0].setPen('w')
                plots[1].setPen('r')
                #TODO: make these colors user defined. create some variable that 
                # the user can update, i.e. self.LRerrColor and self.LRfluxColor

    def whichPlot(self):
        choice, ok = qt.QInputDialog.getItem(self,"Adding region","Which plot?:",np.arange(len(self.plt),dtype=str),0,False)
        self.addLinearRegion(int(choice))

    def addLinearRegion(self,chosen):
        #NOTE:is there a means to extend the length of a list without changing its elements?
        xlen = self.plt[chosen].getViewBox().viewRange()[0]#first element is the bounds of the x-axis data
        #TODO: the logic here is not great as it doesn't delineate linear regions of different plots
        #Need to set linear regions as children of plots (in some sense)
        xran = np.sort(np.random.normal(np.mean(xlen),0.03*(xlen[1]-xlen[0]),size=(2)))
        #NOTE: xran is used to force the bounds of the linear region to be around the 
        # mean of the domain
        self.lrs.append(pg.LinearRegionItem(xran))
        count = len(self.lrs)
        self.plt[chosen].addItem(self.lrs[count-1])
        abdata = self.plt[chosen].listDataItems()
        data = (abdata[0].getData()[0],abdata[0].getData()[1],abdata[1].getData()[1])
        self.addRegionPlot(data)

        self.lrs[len(self.lrs)-1].sigRegionChanged.connect(self.updateLRplot)

    def nonParamEW(self):
        '''
        Module for Non-parameterized Equivalent width determination
        NOTE: the logic assumes you are finding an emission line
        '''

        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        if self.fill is not None: self.plt[dat_choice].removeItem(self.fill)

        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Do not currently support no region choice")
            return
        mask_choice = check.index(float(reg))
        lr = self.lrs[mask_choice].getRegion()


        if len(data) > 0:
            wl = data[0].getData()[0]
            flux = data[0].getData()[1]
            err = data[1].getData()[1]
            reply = qt.QMessageBox.question(self,'pickles','do you want a pickle?',qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)
            if reply == qt.QMessageBox.Yes:
                options = qt.QFileDialog.Options()
                options |= qt.QFileDialog.DontUseNativeDialog
                pick,ok = qt.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","","All Files (*);;pickles (*.pkl)",options = options) 
            if len(self.contfit[0]) > 0:
                items = ("{}".format(i) for i in range(len(self.contfit[0])))
                cont, ok = qt.QInputDialog.getItem(self,"Continuum","Which Continuum?:",self.pdfs.keys(),0,False)
                if ok:
                    #NOTE: this should grab the relevant Probability density functions
                    samples = self.pdfs[cont][0]
                else:
                    qt.QMessageBox.about(self,"Done","Exiting without fit")
                    return
                mask = (wl[:-1] > lr[0]) & (wl[:-1] < lr[1])
                finalflux = flux[mask]
                finalwl = wl[:-1][mask]
                finalerr = err[mask]
                ind = np.where(finalflux == np.max(finalflux))
                maxwl = finalwl[ind]
                cont_pdf = Pow(maxwl,samples['amp'],samples['alpha'],samples['unity'])
                left = qt.QInputDialog.getDouble(self,"left","What velocity left maximum?: ",0,False)
                right = qt.QInputDialog.getDouble(self,"right","What velocity right maximum?: ",0,False)
                newMask = (finalwl > maxwl*(1-(left[0]/3e5))) & (finalwl < maxwl*(1+(right[0]/3e5)))
                #NOTE: this gives consistent results, though they may not be the best choice
                finalflux = finalflux[newMask]
                finalwl = finalwl[newMask]
                finalerr = finalerr[newMask]
                #TODO:fillbetween?
                Limit_curve_top = pg.PlotCurveItem(finalwl,finalflux)
                Limit_curve_bottom = pg.PlotCurveItem(finalwl,finalerr)
                self.fill = pg.FillBetweenItem(Limit_curve_bottom,Limit_curve_top,brush=(0,100,0,150))
                self.plt[dat_choice].addItem(self.fill)


                #Start EW determination
                fflux = [flux - cont_pdf for flux in finalflux]#NOTE: first step in equivalent width determination
                EW_pdf = np.trapz(fflux,x=finalwl,axis=0)/cont_pdf
                EW = np.mean(EW_pdf)
                self.pdfs['npEW'] = [EW_pdf,pd.core.series.Series([EW],index=["EW"])]
                pdf_err = np.std(EW_pdf)
                err_list = np.array([(finalwl[i+1] - finalwl[i])/2 * np.sqrt(finalerr[i+1]**2 + finalerr[i]**2) for i in range(len(finalwl)-1)])
                EWerr = np.sqrt(np.sum(err_list**2) + pdf_err**2)
                qt.QMessageBox.about(self,"Measured","Non-Parameterized EW: {0} \xb1 {1}".format(EW,EWerr))
            elif reply == qt.QMessageBox.Yes:
                cont_params = pd.read_pickle(pick)
                #TODO: two scenarios (currently) 1: use slope m and intercept b to estimate continuum
                # 2: use centroid, cont_amp, and unity to estimate continuum
                # 1: continuum = cont_params['m']*lambda + cont_params['b']
                # 2: continuum = cont_params['cont_amp']*(lamb/cont_params['unity'])**cont_params['alpha']
                # Then, take (integral - continuum)/continuum to get the equivalent width
                mask = (wl[:-1] > lr[0]) & (wl[:-1] < lr[1])
                finalflux = flux[mask]
                finalwl = wl[:-1][mask]
                finalerr = err[mask]
                ind = np.where(finalflux == np.max(finalflux))
                z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
                maxwl = 1215.56845*(1+z)
                print(maxwl)
                if 'cont_amp' not in cont_params.varnames:
                    cont_pdf = Pow(maxwl,cont_params['amp'],cont_params['alpha'],cont_params['unity'])
                else:
                    cont_pdf = Pow(maxwl,cont_params['cont_amp'],cont_params['alpha'],cont_params['unity'])
                left = qt.QInputDialog.getDouble(self,"left","What velocity left maximum?: ",0,False)
                right = qt.QInputDialog.getDouble(self,"right","What velocity right maximum?: ",0,False)
                newMask = (finalwl > maxwl*(1-(left[0]/3e5))) & (finalwl < maxwl*(1+(right[0]/3e5)))
                #NOTE: this gives consistent results, though they may not be the best choice
                finalflux = finalflux[newMask]
                finalwl = finalwl[newMask]
                finalerr = finalerr[newMask]
                Limit_curve_top = pg.PlotCurveItem(finalwl,finalflux)
                Limit_curve_bottom = pg.PlotCurveItem(finalwl,finalerr)#TODO: how to make fillbetween step?
                self.fill = pg.FillBetweenItem(Limit_curve_bottom,Limit_curve_top,brush=(0,100,0,150))
                self.plt[dat_choice].addItem(self.fill)


                #Start EW determination
                fluxes = np.random.normal(finalflux,scale=finalerr,size=(1000,len(finalflux)))#creating 1000 spectra according to distribution from fluxerr
                #This allows the uncertainty in the flux measurements to affect the EW distribution
                EW_pdf = np.zeros((len(fluxes),len(cont_pdf)))
                for i,flux in enumerate(fluxes):#TODO: check logic of this. seems flawed
                    fflux = [f - cont_pdf for f in flux]
                    EW_pdf[i] = np.trapz(fflux,x=finalwl,axis=0)/cont_pdf
                EW_pdf = EW_pdf.flatten()
                EW = np.mean(EW_pdf)
                self.pdfs['npEW'] = [EW_pdf,pd.core.series.Series([EW],index=["EW"])]
                if 'm' in cont_params.varnames: flux_cont = linear(maxwl,cont_params['m'],cont_params['b'])
                else: flux_cont = cont_pdf
                actflux = np.zeros((len(fluxes),len(flux_cont)))
                #embed()
                for i,flux in enumerate(fluxes):
                    tempflux = []
                    for f in flux:
                        tempflux.append(f-flux_cont)
                    actflux[i] = np.trapz(tempflux,x=finalwl,axis=0)
                self.pdfs['Flux'] = [actflux.flatten(),pd.core.series.Series([np.mean(actflux.flatten())],index=['Flux'])]
                pdf_err = np.std(EW_pdf)
                err_list = np.array([(finalwl[i+1] - finalwl[i])/2 * np.sqrt(finalerr[i+1]**2 + finalerr[i]**2) for i in range(len(finalwl)-1)])
                EWerr = np.sqrt(np.sum(err_list**2) + pdf_err**2)
                qt.QMessageBox.about(self,"Measured","Non-Parameterized EW: {0} \xb1 {1}".format(EW,EWerr))
            else:
                qt.QMessageBox.about(self,"No continuum","No continuum available, not measuring")
        else:
            qt.QMessageBox.about(self,"No data","No data to measure")

    #TODO: create something for grabbing filters (a couple hst for now) and 
    #calculating the ratio between spectra and photometry to estimate 
    #slitloss corrections

    def LineCenter(self):
        '''
        Module for determining line center non-parametrically
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()

        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Do not currently support no region choice")
            return
        mask_choice = check.index(float(reg))
        lr = self.lrs[mask_choice].getRegion()


        if len(data) > 0:
            wl = data[0].getData()[0]
            flux = data[0].getData()[1]
            err = data[1].getData()[1]
           
            mask = (wl[:-1] > lr[0]) & (wl[:-1] < lr[1])
            finalflux = flux[mask]
            finalwl = wl[:-1][mask]
            finalerr = err[mask]

            #Start Line center determination
            Lin = np.trapz(finalflux*finalwl/finalerr**2,x=finalwl)/np.trapz(finalflux/finalerr**2,x=finalwl)
            Linerr = np.trapz(finalflux*finalwl**2/finalerr**2,x=finalwl)/np.trapz(finalflux/finalerr**2,x=finalwl)-Lin**2
            F_nu = ((np.sum(1/finalerr**2))**2 - np.sum(1/finalerr**4))/(np.sum(1/finalerr**2))**2
            Linerr = np.sqrt(Linerr/F_nu)
            infline = pg.InfiniteLine(Lin,pen=(100,50,200))
            self.plt[dat_choice].addItem(infline)
            qt.QMessageBox.about(self,"Measured","Line center: {0} \xb1 {1}".format(Lin,Linerr))
            HST_dat = ascii.read('DATA/hst_wfc_475W.txt')
            hst_wl = HST_dat['col1']
            hst_tp = HST_dat['col2']
            newflux = spectres(hst_wl,finalwl,finalflux,fill=0,verbose=False)
            qt.QMessageBox.about(self,"Integral","Flux: {}".format(np.trapz(newflux*hst_tp,x=hst_wl)/np.trapz(hst_tp,x=hst_wl)))
        else:
            qt.QMessageBox.about(self,"No data","No data to measure")

    def fitting_1d(self):
        #TODO: rename to parameterEW or something along these lines to distinguish from brute force integration method
        """
        This will take in user input on the location of some peak and then fit a gaussian to that peak.
        """
        #TODO: here we need to first ask for which plot/dataset, then which continuum, and then which region
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()

        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Do not currently support no region choice")
            return
        mask_choice = check.index(float(reg))
        lr = self.lrs[mask_choice].getRegion()


        if len(data) > 0:
            wl = data[0].getData()[0][:-1]
            flux = data[0].getData()[1]
            err = data[1].getData()[1]

            Ans = qt.QMessageBox.question(self,"Masking","Mask Emission and Absorption lines?",qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)

            if Ans==qt.QMessageBox.Yes:
                z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
                if ok:
                    wlSiII = 1260*(1+z)
                    wlOISiII = 1303*(1+z)
                    wlCII = 1334*(1+z)
                    wlSiIV = 1393*(1+z)
                    wlCIV = 1550*(1+z)
                    wlHeII = 1640*(1+z)#TODO: re-impliment using list of lines instead
                    wlOIII = 1665*(1+z)
                    wlCIII = 1908*(1+z)
                    wlOIV = 1343*(1+z)
                    c = 3e5 #km/s
                    v = c*(wl-wlCIV)/wlCIV
                    CIVMask = (v > -1000) & (v < 400)
                    HeIIMask = (c*(wl-wlHeII)/wlHeII > -400) & (c*(wl-wlHeII)/wlHeII < 400)
                    OIIIMask = (c*(wl-wlOIII)/wlOIII > -400) & (c*(wl-wlOIII)/wlOIII < 400)
                    CIIIMask = (c*(wl-wlCIII)/wlCIII > -400) & (c*(wl-wlCIII)/wlCIII < 400)
                    OIVMask = (c*(wl-wlOIV)/wlOIV > -400) & (c*(wl-wlOIV)/wlOIV < 400)
                    SiIIMask = (c*(wl-wlSiII)/wlSiII > -400) & (c*(wl-wlSiII)/wlSiII < 400)
                    OISiIIMask = (c*(wl-wlOISiII)/wlOISiII > -400) & (c*(wl-wlOISiII)/wlOISiII < 400)
                    CIIMask = (c*(wl-wlCII)/wlCII > -400) & (c*(wl-wlCII)/wlCII < 400)
                    SiIVMask = (c*(wl-wlSiIV)/wlSiIV > -400) & (c*(wl-wlSiIV)/wlSiIV < 400)
                    #add1 = (wl > 4720) & (wl < 4750)
                    #add2 = (wl > 4908) & (wl < 4933)
                    TotMask = CIVMask | HeIIMask | OIIIMask | CIIIMask | OIVMask | SiIIMask | OISiIIMask | CIIMask | SiIVMask#| add1 | add2
                    x = wl[TotMask]
                    y = flux[TotMask]
                    CIVscatter = pg.ScatterPlotItem(x=x,y=y,pen=pg.mkPen('g'),brush=pg.mkBrush('g'))
                    self.plt[dat_choice].addItem(CIVscatter)
                    wl = wl[~TotMask]#~ is bitwise NOT which is how we mask in python
                    flux = flux[~TotMask]
                    err = err[~TotMask]
            else:
                qt.QMessageBox.about(self,"No Mask","Continuing fit with all data in bounds")
           
            mask = (wl > lr[0]) & (wl < lr[1])
            altmask = (wl > lr[0]) & (wl < 4200)
            finalwl = wl[mask] 
            index = np.where(flux[mask] == np.max(flux[mask]))
            peakWl = finalwl[index]
            #peakWl = wl[altmask][index]
            peakFl = flux[mask][index]
            #peakFl = flux[altmask][index]
            #finalwl = wl[mask]
            peakWl, ok = qt.QInputDialog.getDouble(self,"Get Emission line location","Centroid:",3000.0,0.0,10000.0,10)

            zB = ((peakWl - 0.01*peakWl), (peakWl + 0.01*peakWl))
            #zB = (zB[0][0],zB[1][0])#necessary b/c zB is created as array of arrays and numpy fails with array inputs
            sigB = (0.01, 15)
            ampB = (0,4*peakFl)
            ampB = (ampB[0],ampB[1][0])#same as z
            #embed()
            items = ('Gauss','Voigt')
            func, ok = qt.QInputDialog.getItem(self,"Get function","Function: ",items,0,False)
            if func == 'Gauss':
                self.Fitter(Powpgauss,data,flux[mask],err[mask],finalwl,[ampB,zB,sigB,(0,np.max(flux[mask])),(-5,5),(lr[0],lr[1])],name='EW',plt_name=dat_choice)
            elif func == 'Voigt':
                gvamp = (0,np.max(flux[mask]))
                gtau = np.array(zB) + 200
                diffwl = np.diff(wl)
                index = index[0][0]
                f0 = flux[mask][index+8]
                wl0 = finalwl[index+8]
                f1 = flux[mask][index+50]
                wl1 = finalwl[index+50]
                slope = abs(f1-f0)/abs(wl1-wl0)
                self.Fitter(piecewise,data,flux[mask],err[mask],finalwl,[(lr[0]+100,lr[0]+200),(lr[0],lr[1]),zB,(0,np.max(flux[mask])),ampB,(-5,5),sigB,(-1000,1000),(slope-0.2*slope,slope+0.2*slope)],name='Voigt',plt_name=dat_choice)
        else:
            qt.QMessageBox.about(self,"No data on screen","Not fitting")
                        
    def Cont_Fit(self):
        """
        This will ask for a wavelength range from the user. A continuum will be fit
        to the flux data with consideration of the error spectrum. Should also consider setting up
        a list of masks to the data to use for continuum fitting
        """
        #TODO: must ask 2 questions. 1: which plot? 2:which linear region? (this should be limited to those within the plot)
        # 1 gives the data to be analyzed, 2 gives the masking region (multiple regions?)
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return

        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Do not currently support no region choice")
            return

        #geting data for fitting continuum and plotting chosen points
        mask_choice = check.index(float(reg))
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        lr = self.lrs[mask_choice].getRegion()
        mask = (data[0].getData()[0][:-1] > lr[0]) & (data[0].getData()[0][:-1] < lr[1])
        flux = data[0].getData()[1][mask]
        wl = data[0].getData()[0][:-1][mask]
        err = data[1].getData()[1][mask]
        
        Ans = qt.QMessageBox.question(self,"Masking","Mask Emission and Absorption lines?",qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)

        if Ans==qt.QMessageBox.Yes:
            z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
            if ok:
                wlSiII = 1260*(1+z)
                wlOISiII = 1303*(1+z)
                wlCII = 1334*(1+z)
                wlSiIV = 1393*(1+z)
                wlCIV = 1550*(1+z)
                wlHeII = 1640*(1+z)#TODO: re-impliment using list of lines instead
                wlOIII = 1665*(1+z)
                wlCIII = 1908*(1+z)
                c = 3e5 #km/s
                v = c*(wl-wlCIV)/wlCIV
                CIVMask = (v > -1000) & (v < 400)
                SiIIMask = (c*(wl-wlSiII)/wlSiII > -400) & (c*(wl-wlSiII)/wlSiII < 400)
                OISiIIMask = (c*(wl-wlOISiII)/wlOISiII > -400) & (c*(wl-wlOISiII)/wlOISiII < 400)
                CIIMask = (c*(wl-wlCII)/wlCII > -400) & (c*(wl-wlCII)/wlCII < 400)
                SiIVMask = (c*(wl-wlSiIV)/wlSiIV > -400) & (c*(wl-wlSiIV)/wlSiIV < 400)
                HeIIMask = (c*(wl-wlHeII)/wlHeII > -400) & (c*(wl-wlHeII)/wlHeII < 400)
                OIIIMask = (c*(wl-wlOIII)/wlOIII > -400) & (c*(wl-wlOIII)/wlOIII < 400)
                CIIIMask = (c*(wl-wlCIII)/wlCIII > -400) & (c*(wl-wlCIII)/wlCIII < 400)
                #add1 = (wl > 4819) & (wl < 4981)
                #add2 = (wl > 5162) & (wl < 5305)
                TotMask = CIVMask | HeIIMask | OIIIMask | CIIIMask | SiIIMask | OISiIIMask | CIIMask | SiIVMask#| add1 | add2
                x = wl[TotMask]
                y = flux[TotMask]
                CIVscatter = pg.ScatterPlotItem(x=x,y=y,pen=pg.mkPen('g'),brush=pg.mkBrush('g'))
                self.plt[dat_choice].addItem(CIVscatter)
                wl = wl[~TotMask]#~ is bitwise NOT which is how we mask in python
                flux = flux[~TotMask]
                err = err[~TotMask]
        else:
            qt.QMessageBox.about(self,"No Mask","Continuing fit with all data in bounds")
        items = ('Power Law','Linear','Polynomial')
        func, ok = qt.QInputDialog.getItem(self,"Get function","Function: ",items,0,False)
        if func == 'Power Law' and ok:
            #TODO: Force left bound as "centroid" in Pow?
            self.Fitter(Pow,data,flux,err,wl,[(0,np.max(flux)),(-3,3),(np.min(wl),np.max(wl))],name='continuum',plt_name=dat_choice,cname=func)
            self.cname = 'pl'
            self.contfit[1].append('pl')
        elif func == 'Linear' and ok:
            self.line = partial(linear,b=lr[0])
            self.Fitter(self.line,data,flux,err,wl,[(np.min(flux),np.max(flux)),(np.min(flux)/np.max(wl),np.max(flux)/np.min(wl))],name='continuum',plt_name=dat_choice,cname=func)
            #TODO: is the slope range reasonable?
            self.cname = 'l'
            self.contfit[1].append('l')
        elif func == 'Polynomial' and ok:
            order, check = qt.QInputDialog.getInt(self,"Polynomial","What Order?:",0,False)
            if check:
                self.cname = 'p'
                self.contfit[1].append('p')
                guess = []
                guess.append((np.min(flux),np.max(flux)))
                for i in range(order):
                    guess.append(((np.min(flux)/np.max(wl))**(i+1),(np.max(flux)/np.min(wl)**(i+1))))
                self.Fitter(polynomial,data,flux,err,wl,guess,name='continuum',plt_name=dat_choice,cname=func,order=order)
        else:#TODO: consider using else, this is kind of weird coding
            qt.QMessageBox.about(self,"Continuum","Not fitting continuum")
        
                
        return
        
    #TODO:Need to be able to mask certain data points. maybe instead of z, just ask for data in some bounds?
    def Fitter(self,func,data,flux,err,wl,bounds,name = '',plt_name=None,cname=None,order=None):
        '''
        helper function for fitting algorithms (cont_fit, and ew)
        '''
        #TODO: Threading for progress bar? Maybe just add in self.MCMC code block?
        #embed()
        basic_model = pm.Model()
        if np.max(wl) - np.min(wl) > 500:
            midwl = np.mean(wl)
            leftun = midwl - 250
            rghtun = midwl + 250
        else:
            leftun = np.min(wl)
            rghtun = np.max(wl)
        with basic_model:
            
            #Priors on model
            #NOTE: these are elements in vals in the order they appear in the code
            if name == "continuum":
                if cname == "Power Law":
                    amp = pm.TruncatedNormal("amp",mu=(bounds[0][0]+bounds[0][1])/2,sigma=0.8*(bounds[0][1]-bounds[0][0]),testval=(bounds[0][0]+bounds[0][1])/2,lower=0.0000001)
                    alpha = pm.TruncatedNormal("alpha",mu=(bounds[1][0]+bounds[1][1])/2,sigma=0.8*(bounds[1][1]-bounds[1][0]),testval=(bounds[1][0]+bounds[1][1])/2,lower=-5,upper=5)
                    unity = pm.TruncatedNormal("unity",mu=(bounds[2][0]+bounds[2][1])/2,sigma=0.8*(bounds[2][1]-bounds[2][0]),testval=(bounds[2][0]+bounds[2][1])/2,lower=leftun,upper=rghtun)
                    #step = pm.HamiltonianMC()
                    #Expected value
                    mu = func(wl.astype(np.float32),amp,alpha,unity)
                if cname == "Linear":
                    incpt = pm.Uniform("incpt",bounds[0][0],bounds[0][1])
                    slope = pm.Uniform("slope",bounds[1][0],bounds[1][1])

                    mu = func(wl.astype(np.float32),incpt,slope)
                if cname == "Polynomial":
                    a = []
                    for i in range(order):
                        a.append(pm.Uniform("a{}".format(i),bounds[i][0],bounds[i][1]))
                    mu = func(wl.astype(np.float32),*a)
            if name == "EW":
                #Priors
                #embed()
                amp = pm.TruncatedNormal("amp",mu=(bounds[0][0]+bounds[0][1])/2,sigma=0.8*(bounds[0][1] - bounds[0][0]),testval=bounds[0][1]/2,lower=0)
                centroid = pm.Normal("centroid",mu=(bounds[1][0]+bounds[1][1])/2,sigma=0.8*(bounds[1][1] - bounds[1][0]))
                sigma = pm.TruncatedNormal("sigma",mu=(bounds[2][0]+bounds[2][1])/2,sigma=0.4*(bounds[2][1] - bounds[2][0]),testval=(bounds[2][0]+bounds[2][1])/2,lower=0)
                cont_amp = pm.TruncatedNormal("cont_amp",mu=(bounds[3][0]+bounds[3][1])/2,sigma=0.8*(bounds[3][1] - bounds[3][0]),testval=(bounds[3][0]+bounds[3][1])/2,lower=0.000001)
                alpha = pm.TruncatedNormal("alpha",mu=(bounds[4][0]+bounds[4][1])/2,sigma=0.8*(bounds[4][1] - bounds[4][0]),testval=(bounds[4][0]+bounds[4][1])/2,lower=-5,upper=5)
                unity = pm.TruncatedNormal("unity",mu=(bounds[5][0]+bounds[5][1])/2,sigma=0.2*(bounds[5][1] - bounds[5][0]),testval=(bounds[5][0]+bounds[5][1])/2,lower=leftun,upper=rghtun)

                #embed()
                #TODO: consider using non-centered reparameterization, i.e. amp = mu + sigma*amp_0, where amp_0 ~ N(0,1)
                #use 540 as example case
                #step = pm.HamiltonianMC()

                mu = func(wl.astype(np.float32),amp,centroid,sigma,cont_amp,alpha,unity)
            if name == "Voigt":
                switch = pm.TruncatedNormal("switch",mu=(bounds[0][0]+bounds[0][1])/2,sigma=0.8*(bounds[0][1] - bounds[0][0]),testval=(bounds[0][0]+bounds[0][1])/2,lower=bounds[0][0])
                unity = pm.TruncatedNormal("unity",mu=(bounds[1][0]+bounds[1][1])/2,sigma=0.2*(bounds[1][1] - bounds[1][0]),testval=(bounds[1][0]+bounds[1][1])/2,lower=leftun,upper=rghtun)
                centroid = pm.TruncatedNormal("centroid",mu=(bounds[2][0]+bounds[2][1])/2,sigma=0.8*(bounds[2][1] - bounds[2][0]),lower=bounds[2][0]-10,upper=bounds[2][1]+10)
                cont_amp = pm.TruncatedNormal("cont_amp",mu=(bounds[3][0]+bounds[3][1])/2,sigma=0.8*(bounds[3][1] - bounds[3][0]),testval=(bounds[3][0]+bounds[3][1])/2,lower=0.000001)
                amp = pm.TruncatedNormal("amp",mu=(bounds[4][0]+bounds[4][1])/2,sigma=0.8*(bounds[4][1] - bounds[4][0]),testval=(bounds[4][0]+bounds[4][1])/2,lower=0)
                alpha = pm.TruncatedNormal("alpha",mu=(bounds[5][0]+bounds[5][1])/2,sigma=0.8*(bounds[5][1] - bounds[5][0]),testval=(bounds[5][0]+bounds[5][1])/2,lower=-5,upper=5)
                sigma = pm.TruncatedNormal("sigma",mu=(bounds[6][0]+bounds[6][1])/2,sigma=0.4*(bounds[6][1] - bounds[6][0]),testval=(bounds[6][0]+bounds[6][1])/2,lower=0)
                b = pm.Normal("b",mu=(bounds[7][0]+bounds[7][1])/2,sigma=0.4*(bounds[7][1] - bounds[7][0]),testval=(bounds[7][0]+bounds[7][1])/2)
                m = pm.TruncatedNormal("m",mu=(bounds[8][0]+bounds[8][1])/2,sigma=0.4*(bounds[8][1] - bounds[8][0]),testval=(bounds[8][0]+bounds[8][1])/2,lower=0.0)
                mu = func(wl,switch,unity,centroid,cont_amp,amp,alpha,sigma,b,m)
                #mu = pm.math.switch(wl <= switch,polygauss(wl,amp,centroid,sigma,b,m),Pow(wl,cont_amp,alpha,unity))
                
            '''    
            print(basic_model.test_point)
            print(basic_model.check_test_point())
            q0 = step._logp_dlogp_func.dict_to_array(basic_model.test_point)
            p0 = step.potential.random()
            start = step.integrator.compute_state(q0,p0)
            print(start.energy)
            logp,dlogp = step.integrator._logp_dlogp_func(q0)
            print(logp)
            print(dlogp)
            '''
            #embed()
            #Likelihood of sampling distribution
            Y_obs = pm.Normal("Y_obs",mu=mu,sigma=err.astype(np.float32),observed=flux.astype(np.float32))
            step = self.sampler.currentText()
            cores = multi.cpu_count()
            if step == 'Metropolis':
                trace = pm.sample(20000,tune=25000,cores=10,init='adapt_diag',step=pm.step_methods.Metropolis())#used for testing parameter space
            elif (cores >= 10) and (step == 'NUTS'):
                trace = pm.sample(10000,tune=5000,target_accept=0.8,cores=10,init='adapt_diag')
            elif step == 'NUTS':
                trace = pm.sample(10000,tune=5000,target_accept=0.8,cores=cores,init='adapt_diag')#TODO: might be overworking computers here
            else:
                qt.QMessageBox.about(self,"Not supported","sampling method {} is not yet supported. Exiting without fit.".format(step))
                return
            #NOTE: vals['mean'].keys() gives the parameter names
            if cname == "Power Law":
                vals = az.summary(trace,round_to=10,var_names=['amp','alpha','unity'])
            elif name == "EW":
                vals = az.summary(trace,round_to=10,var_names=['amp','centroid','sigma','cont_amp','alpha','unity'])
            elif name == "Voigt":
                vals = az.summary(trace,round_to=10,var_names=['switch','unity','centroid','cont_amp','amp','alpha','sigma','b','m'])
            samples = pm.trace_to_dataframe(trace,varnames=vals['mean'].keys())

        if name == 'continuum':
            self.contfit[0].append(vals['mean'])
            name = self.cname + self.names1d[plt_name]
            if name in self.pdfs.keys():
                name = name + ".1" #TODO: need to allow for more names, though probably want to limit to not eat up all the memory
            self.pdfs[name] = (samples,vals['mean'])
            pen = (50,200,50) # sets the color of the line, TODO: Should make this user adjustable
            cont = 1.0
        
        if name == 'EW':
            self.ewfit.append(vals['mean'])
            name += self.names1d[plt_name]
            if name in self.pdfs.keys():
                name += ".1"
            self.pdfs[name] = (samples,vals['mean'])
            ewPdf = np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])/Pow(vals['mean']['centroid'],samples['cont_amp'],samples['alpha'],samples['unity'])
            ewMeas = np.sqrt(2*np.pi)*vals['mean']['amp']*vals['mean']['sigma']/Pow(vals['mean']['centroid'],vals['mean']['cont_amp'],vals['mean']['alpha'],vals['mean']['unity'])
            self.pdfs['EWPDF'] = [ewPdf,pd.core.series.Series([ewMeas],index=["EW"])]
            self.pdfs['Flux'] = [np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])]
            pen = (0,100,0)
            cont = 1.0
            #This is used for GUI image such that we can see the fit 
        if name == "Voigt":
            self.ewfit.append(vals['mean'])
            name += self.names1d[plt_name]
            if name in self.pdfs.keys():
                name += ".1"
            self.pdfs[name] = (samples,vals['mean'])
            ewPdf = np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])/Pow(vals['mean']['centroid'],samples['cont_amp'],samples['alpha'],samples['unity'])
            ewMeas = np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])/Pow(vals['mean']['centroid'],vals['mean']['cont_amp'],vals['mean']['alpha'],vals['mean']['unity'])
            self.pdfs['EWPDF'] = [ewPdf,pd.core.series.Series([ewMeas],index=["EW"])]
            self.pdfs['Flux'] = [np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])]
            pen = (0,100,0)
            cont = 1.0
            #TODO: plotting of function failing, unsure why. Also EW not being plotted because name doesn't include "EW"
        self.arviz[name] = trace 
        self.plt[plt_name].plot(data[0].getData()[0],cont*func(data[0].getData()[0],*vals['mean']),pen=pen)
        del(basic_model)

    def S_N(self):
        choice, ok = qt.QInputDialog.getItem(self,"Which spectrum?","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not checking","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Do not currently support no region choice")
            return
        mask_choice = check.index(float(reg))
        lr = self.lrs[mask_choice].getRegion()
        if len(data) > 0:
            wl = data[0].getData()[0]
            flux = data[0].getData()[1]
            err = data[1].getData()[1]
            mask = (wl[:-1] > lr[0]) & (wl[:-1] < lr[1])
            finalflux = flux[mask]
            finalwl = wl[:-1][mask]
            finalerr = err[mask]
            S_N_perpix = np.mean(finalflux/finalerr)
            totflux = np.sum(finalflux)
            toterr = np.sqrt(np.sum(finalerr**2))
            qt.QMessageBox.about(self,"S/N ratio","S/N per wavelength bin: {0}. S/N total {1}".format(S_N_perpix,totflux/toterr))

    
    def pickle_plot(self):
        basic_model = pm.Model()
        with basic_model:
            choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
            if not(ok):
                qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
                return
            options = qt.QFileDialog.Options()
            options |= qt.QFileDialog.DontUseNativeDialog
            pick,ok = qt.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","","pickles (*.pkl);;All Files (*)",options = options) 
            dat_choice = self.names1d.index(choice)
            data = self.plt[dat_choice].listDataItems()
            cont_params = pd.read_pickle(pick)
            vals = az.summary(cont_params,round_to=10)
            curve, ok = qt.QInputDialog.getItem(self,"what kind of pickle?","Choose curve:",['gauss','voigt','power','linear','polynomial'],0,False)
            if curve == 'gauss':
                func = Powpgauss
            elif curve == 'voigt':
                #TODO: this doesn't really work for the voigt profile
                func = piece
            elif curve == 'power':
                func = Pow
            elif curve == 'linear':
                func = linear
            elif curve == 'polynomial':
                func = polynomial
            pen = (0,100,0)
            embed()
            self.plt[dat_choice].plot(data[0].getData()[0],func(data[0].getData()[0],*vals['mean']),pen=pen)
        
        
    #TODO: Should consider allowing user to adjust parameterization
    # but should make this a visual thing like a slider bar (It would be helpful to output the sigma
    # compared with the original choice for sigma in the normal distribution)
    #TODO: I also want to show the corner plots for the parameterization
    # and allow the user to hover over the space and select a point
    # this point would then be used to plot the fit to the spectrum.
    # This would mostly be an educational thing, though it may be useful
    # for teasing out problem parameter spaces as well.
    def save_spectrum(self):
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        wl = data[0].getData()[0][:-1]
        flux = data[0].getData()[1]
        err = data[1].getData()[1]

        new_data = np.rec.array([flux,err,wl],
                    formats='float32,float32,float32',names='flux,sig,wave')
        bintable = fits.BinTableHDU(new_data)
        name, ok = qt.QInputDialog.getText(self,"Filename","Name the fits file:")
        dir_ = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtGui.QFileDialog.ShowDirsOnly)
        if len(dir_) == 0:
            return
        cwd = os.getcwd()
        os.chdir(dir_)        
        if not(ok):
            name = "basic"
        bintable.writeto(name + ".fits")
        os.chdir(cwd)  


    def save_data(self):
        '''
        routine to save parameter fits to .fits file
        '''
        basic_model = pm.Model()
        with basic_model:
            choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
            if not(ok):
                return
            dat_to_save = self.arviz[choice]
            name, ok = qt.QInputDialog.getText(self,"Filename","Name the fits file:")
            if not(ok):
                return
            dir_ = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtGui.QFileDialog.ShowDirsOnly)
            if len(dir_) == 0:
                return
            vals = az.summary(self.arviz[choice],round_to=5,stat_funcs=[conf_low,conf_high])
            cwd = os.getcwd()
            os.chdir(dir_)
            if os.path.exists(name+'.pkl'):
                os.remove(name+'.pkl')
            if os.path.exists(name+'.txt'):
                os.remove(name+'.txt')
            with open(name+'.pkl','wb') as output:
                pickle.dump(dat_to_save,output,pickle.HIGHEST_PROTOCOL)
            with open(name+'.txt','a') as output:
                for i in range(len(vals['mean'])):
                    output.write("{0}   {1}   {2}\n".format(vals['mean'][i],vals['conf_low'][i]-vals['mean'][i],vals['conf_high'][i]-vals['mean'][i]))
            os.chdir(cwd)
        del(basic_model)
                
    def arviz_density(self):
        """
        show density plot
        """
        basic_model = pm.Model()
        with basic_model:
            choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
            if not(ok):
                return
            az.plot_density(self.arviz[choice])
            plt.show()
        del(basic_model)

    def arviz_autocorrelation(self):
        '''
        show autocorrelation plot
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_autocorr(self.arviz[choice])
        plt.show()

    def arviz_energy(self):
        '''
        show energy plot
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_energy(self.arviz[choice])
        plt.show()

    def arviz_forest(self):
        '''
        show forest plot
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_forest(self.arviz[choice])
        plt.show()

    def arviz_joint(self):
        '''
        essentially a corner plot with marginals=True
        '''
        #TODO: Not able to select non-parameterized pdf if fits performed
        basic_model = pm.Model()
        with basic_model:
            if len(self.arviz.keys()) > 0:
                choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
                if not(ok):
                    return
                vals = az.summary(self.arviz[choice],round_to=5,stat_funcs=[conf_low,conf_high])#var_names?
                az.rcParams['plot.max_subplots'] = 80
                axes = az.plot_pair(self.arviz[choice],kind=['hexbin',"kde"],kde_kwargs={"fill_last":False},marginals=True,point_estimate="mean",marginal_kwargs={"quantiles":[0.16,0.5,0.84]},figsize=(len(vals['mean']),len(vals['mean'])))
                for i in range(len(axes)):
                    axes[i][i].set_title(r'${0}={1}_{{{2}}}^{{{3}}}$'.format(vals['mean'].keys()[i],vals['mean'][i],vals['conf_low'][i]-vals['mean'][i],vals['conf_high'][i]-vals['mean'][i]))
                if (choice.find('EW') != -1) or (choice.find('Voigt') != -1):    
                    fig2 = plt.figure()
                    ax2 = plt.axes()
                    y,x,_ = ax2.hist(self.pdfs['EWPDF'][0],bins=100,density=True)
                    #embed()
                    up = conf_high(self.pdfs['EWPDF'][0])
                    down = conf_low(self.pdfs['EWPDF'][0])
                    m = np.mean(self.pdfs['EWPDF'][0])
                    ax2.vlines(m,0,np.max(y),color='k')
                    ax2.text(np.mean(x) - np.std(x),0.8*np.max(y),r'${0}^{{{1}}}_{{{2}}}$'.format(m,up-m,down-m))
                    ax2.set_xlabel('EW')
                    fig3 = plt.figure()
                    ax3 = plt.axes()
                    fy,fx,_ = ax3.hist(self.pdfs['Flux'][0],bins=100,density=True)
                    fup = conf_high(self.pdfs['Flux'][0])
                    fdown = conf_low(self.pdfs['Flux'][0])
                    fm = np.mean(self.pdfs['Flux'][0])
                    ax3.vlines(fm,0,np.max(fy),color='k')
                    ax3.text(np.mean(fx)-np.std(fx),0.8*np.max(fy),r'${0}^{{{1}}}_{{{2}}}$'.format(fm,fup-fm,fdown-fm))
            else:
                choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.pdfs.keys(),0,False)
                fig1 = plt.figure()
                ax1 = plt.axes()
                pandadat = pd.DataFrame(self.pdfs[choice][0])#TODO: dataframes are much faster, good idea to
                #port all math into pandas dataframes
                y,x,_ = ax1.hist(pandadat,bins=100,density=True)
                up = pandadat.quantile(q=0.841,interpolation='linear')[0]
                down = pandadat.quantile(q=0.159,interpolation='linear')[0]
                m = pandadat.quantile(q=0.5,interpolation='linear')[0]
                ax1.vlines(m,0,np.max(y),color='k')
                ax1.text(np.mean(x) - np.std(x),0.8*np.max(y),r'${0}^{{{1}}}_{{{2}}}$'.format(m,up-m,down-m),fontsize=25)
            plt.show()
        del(basic_model)

    def arviz_parallel(self):
        '''
        show parallel plot
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_parallel(self.arviz[choice])
        plt.show()

    def arviz_trace(self):
        '''
        show outcome of run
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_trace(self.arviz[choice])
        plt.show()

    def arviz_plot_posterior(self):
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_posterior(self.arviz[choice])
        plt.show()

    def Flux_to_Lum(self):

        choice, ok = qt.QInputDialog.getItem(self,"Which plot","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not converting","Chose no data, not converting.")
            return
        z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
        if not(ok):
            qt.QMessageBox.about(self,"No redshift","can't convert without redshift")
            return
        
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        wl = data[0].getData()[0]
        measurements = data[0].getData()[1]
        measurement_errs = data[1].getData()[1]

        H0 = 70.0
        Om0 = 0.3
        Ob0 = 0.05
        Tcmb0 = 2.725
    

        cosmo = cosmology.FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0, Tcmb0=Tcmb0)



        flux_to_lum = lambda flux, lum_dist: np.multiply(np.multiply(flux, np.square(lum_dist)), 4.*np.pi)  ## bolometric fluxes and luminosities
        lum_to_flux = lambda lum, lum_dist:  np.divide(np.divide(lum, np.square(lum_dist)), 4.*np.pi)

        cm_in_Mpc = 3.086e24

        obj_lum_dist_Mpc = cosmo.luminosity_distance(z).value
        obj_lum_dist_cm  = obj_lum_dist_Mpc * cm_in_Mpc


        ## The (1+redshift) factor accounts for the fact that flux and luminosity are densities per unit WAVELENGTH
        ## The (1+redshift) factor would not be applied to bolometric luminosities and fluxes
        ## The (1+redshift) factor would be in the denominator / numerator respectively if densities were per unit frequency
        ## See Hogg+2002 K-correction paper

        self.plt[dat_choice].clear()
        converted_meas = flux_to_lum(measurements, obj_lum_dist_cm)
        conv_meas_errs = flux_to_lum(measurement_errs, obj_lum_dist_cm)
        converted_meas = np.multiply(converted_meas, 1.+ z)
        conv_meas_errs = np.multiply(conv_meas_errs, 1.+ z)
        self.plt[dat_choice].plot(wl,converted_meas,pen='b',stepMode=True)
        self.plt[dat_choice].plot(wl,conv_meas_errs,stepMode=True)


        return converted_meas, conv_meas_errs


    def plotColoring(self,isFlux = False, isErr=False):
        choice, ok = qt.QInputDialog.getItem(self,"Which to color?","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Coloring","Chose no data, not Coloring.")
            return
        dat_choice = self.names1d.index(choice)
        color = qt.QColorDialog.getColor()
        if not(self.flux==None) and isFlux:
            self.flux[dat_choice].setPen(color)
        if not(self.err==None) and isErr:
            self.err[dat_choice].setPen(color)

    def upperLimit(self):
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
        Lya = 1215.56845*(1+z)#Caroll and Ostlie second edition equation 5.7
        width, ok = qt.QInputDialog.getDouble(self,"Get width",r"$\Delta w:$",2.0,0.0,20.0,10)
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()
        wl = data[0].getData()[0][:-1]
        err = data[1].getData()[1]
        diffwl = np.diff(wl)
        diffwl = np.append(diffwl,2.18)
        mask = (wl >= Lya - width/2) & (wl <= Lya + width/2)
        Flux5sig = 5*np.sqrt(np.sum(err[mask]**2*diffwl[mask]**2))
        qt.QMessageBox.about(self,"Upper Limit",r"The 5$\sigma$ upper limit is {}".format(Flux5sig))

    #TODO: broken, got to deal with len(y) + 1 = len(x) details for step plot
    def artificial_gauss(self):
        '''
        Module for adding gaussians to existing data.
        intended for sensitivity measurements
        ''' 
        choice, ok = qt.QInputDialog.getItem(self,"Which to measure?","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Adding","Chose no data, not adding.")
            return
        dat_choice = self.names1d.index(choice)
        data = self.plt[dat_choice].listDataItems()

        check = [self.regPlot[i].getViewBox().viewRange()[0][0] for i in np.arange(len(self.regPlot))]
        reg, ok = qt.QInputDialog.getItem(self,"Choose a region","Which plot?:",np.array(check,str),0,False)#TODO: allow for multiple regions? or only abs masking?
        if not(ok):
            qt.QMessageBox.about(self,"Not Adding","Do not currently support no region choice")
            return

        mask_choice = check.index(float(reg))
        lr = self.lrs[mask_choice].getRegion()
        mask = (data[0].getData()[0][:-1] > lr[0]) & (data[0].getData()[0][:-1] < lr[1])
        origFlux = data[0].getData()[1]
        origWl = data[0].getData()[0]
        flux = data[0].getData()[1][mask]
        wl = np.array(data[0].getData()[0][mask])
        err = data[1].getData()[1][mask]

        sig, ok = qt.QInputDialog.getDouble(self,"Gaussian Width", "What is the width of the gaussians?:",0.0,0,50.0,5)
        if not(ok):
            qt.QMessageBox.about(self,"Not Adding","Do not currently support no sigma choice")
        FWHM = 2*np.sqrt(2*np.log(2))*sig
        #pix_count = FWHM/2.18 #from QA of pypeit, only for Chris
        #sig_data = np.sqrt(pix_count)*err
        sig_data = []
        sig_mean = []
        diffwl = np.diff(wl)
        diffwl = np.append(diffwl, 2.18)#from QA of pypeit, only for Chris
        for w in wl:
            masker = (wl > w -FWHM/2) & (wl < w + FWHM/2)
            sig_data.append(np.sqrt(np.sum(err[masker]**2*diffwl[masker]**2)))
            sig_mean.append(np.mean(flux[masker]))
        sig_data = np.array(sig_data)
        sig_mean = np.array(sig_mean)
        
        #pen1sig = pg.mkPen((139,0,0),style=QtCore.Qt.DashLine)#darkred
        #pen2sig = pg.mkPen((128,128,0),style=QtCore.Qt.DashLine)#olive
        #pen3sig = pg.mkPen((0,100,0),style=QtCore.Qt.DashLine)#darkgreen

        done = "False"
        run = 0
        Flux = []
        success = []
        while not(done == "True"):
            y = origFlux[mask]
            data[0].setData(origWl,origFlux,pen='b')
            #if run == 0: 
            #    self.plt[dat_choice].plot(wl,sig_data+sig_mean,pen=pen1sig,name='1 sigma')
            #    self.plt[dat_choice].plot(wl,2*sig_data+sig_mean,pen=pen2sig,name='2 sigma')
            #    self.plt[dat_choice].plot(wl,3*sig_data+sig_mean,pen=pen3sig,name='3 sigma')
            #    self.plt[dat_choice].addLegend()
            pause, ok = qt.QInputDialog.getInt(self,'pausing','Data Reset',0,0,0,0)
            run += 1
            #x_0 = np.random.uniform(lr[0],lr[1])
            #i = (np.abs(wl-x_0)).argmin()#rough location for c
            i = np.random.choice(wl.shape[0],1)
            x_0 = wl[i]
            fl = np.random.uniform(0.1,2.0)
            Flux.append(fl)
            A = fl/(np.sqrt(2*np.pi)*sig)
            gauss = gauss0(wl,A,x_0,np.sqrt(2)*sig)
            y += gauss
            data[0].setData(wl,y)
            win,ok = qt.QInputDialog.getInt(self,"Flux: {0}, 1-sig: {1}, 2-sig: {2}, 3-sig: {3}, 5-sig: {4}".format(fl,sig_data[i],2*sig_data[i],3*sig_data[i],5*sig_data[i]),"did you find it?(1=yes,0=no): ",0,0,1,1)
            success.append(win)
            done,ok = qt.QInputDialog.getItem(self,"Total runs: {}".format(run),"Done?:",["True","False"],0,False)
        data[0].setData(origWl,origFlux,pen='b')
        bool_list = list(map(bool,success))
        Flux = np.array(Flux)
        qt.QMessageBox.about(self,"Finished","You identified {0} out of {1} gaussians with mininum found {2}".format(np.sum(success),len(success),np.min(Flux[bool_list])))
        #EW = sqrt(2pi)*A*sigma/Cont(x_0)
    
    def Remove_lines(self):
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        dat_choice = self.names1d.index(choice)
        for L in self.lines:
            self.plt[dat_choice].removeItem(L)
        self.lines = []
        
    def Show_lines(self):
        '''
        Method for showing all emission/absorption lines
        based on redshift given by user
        '''
        choice, ok = qt.QInputDialog.getItem(self,"Which to fit","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Fitting","Chose no data, not fitting.")
            return
        z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
        if not(ok):
            qt.QMessageBox.about(self,"Not plotting","Chose no data, not plotting.")
            return
        dat_choice = self.names1d.index(choice)
        Lyman = 1216
        oLy = (z+1)*np.array(1216)
        Neblines = [1265,1309,1533,1661,1666,1909] #SiII*1265,1309,1533;OIII]1661,1666;CIII]1909
        Lowion = [1260.4221,1302.1685,1304.3702,1334.5323,1526.70698,1608.45085,1670.7886] #SiII1260;OI1302;SiII1304;CII1334;SiII1526;FeII1608;AlII1670
        Highion = [1393.76018,1402.77291,1548.2049,1550.77845] #SiIV1393,1402;CIV1548,1550   
        Steidel04 = [1334.5323,1854.71829,1862.79113,2344.2129601,2374.4603294,2382.7641781,2586.6495659,2600.1724835,2796.3542699,2803.5314853]# CII1334;AlIII1854,1862;FeII2344,2374,2382;FeII2586,2600;MgII2796,2803
        oHighion = (z+1)*np.array(Highion)
        oLowion = (z+1)*np.array(Lowion)
        oNeblines = (z+1)*np.array(Neblines)
        oSteidel04 = (z+1)*np.array(Steidel04)
        #TODO: add method for removing these lines, also want to lock the movability of these lines to each other
        #Actually it would be interesting to lock all the nebular lines seperately from high/low and vice versa
        Lyline = pg.InfiniteLine(angle=90,movable=True,hoverPen='g',pen=(0,120,0),label=str(int(Lyman)),markers='o',name='High ionization lines',labelOpts={'rotateAxis':(1,0),'color':'k'})
        self.lines.append(Lyline)
        self.plt[dat_choice].addItem(Lyline)
        Lyline.setPos(oLy)
        for count,H in enumerate(oHighion):
            Iline = pg.InfiniteLine(angle=90,movable=True,hoverPen='g',pen=(0,120,60),label=str(int(Highion[count])),markers='o',name='High ionization lines',labelOpts={'rotateAxis':(1,0),'color':'k'})
            self.plt[dat_choice].addItem(Iline)
            self.lines.append(Iline)
            Iline.setPos(H)
        for count,L in enumerate(oLowion):
            Iline = pg.InfiniteLine(angle=90,movable=True,hoverPen='g',pen=(120,0,60),label=str(int(Lowion[count])),markers='o',name='Low-Ionization Absorption lines',labelOpts={'rotateAxis':(1,0),'color':'k'})
            self.plt[dat_choice].addItem(Iline)
            self.lines.append(Iline)
            Iline.setPos(L)
        for count,N in enumerate(oNeblines):
            Iline = pg.InfiniteLine(angle=90,movable=True,hoverPen='g',pen=(120,60,0),label=str(int(Neblines[count])),markers='o',name='Nebular lines',labelOpts={'rotateAxis':(1,0),'color':'k'})
            self.plt[dat_choice].addItem(Iline)
            self.lines.append(Iline)
            Iline.setPos(N)
        for count,S in enumerate(oSteidel04):
            Iline = pg.InfiniteLine(angle=90,movable=True,hoverPen='g',pen=(0,60,120),label=str(int(Steidel04[count])),markers='o',name='Steidel 2004',labelOpts={'rotateAxis':(1,0),'color':'k'})
            self.plt[dat_choice].addItem(Iline)
            self.lines.append(Iline)
            Iline.setPos(S)
        self.plt[dat_choice].addLegend()
    #TODO: generalize, currently constructed for pypeit 1.3.3
    def coadd_1d(self,files):
        hdul = []
        a = []
        for f in files:
            
            hdul.append(fits.open(f))
            if f.split('/')[-1] == "1212_M114":
                a.append(1.28)
            elif f.split('/')[-1] == "1212_M1":
                a.append(1.12)
            else:
                a.append(1)
        resdata = []
        new_grid = np.arange(3000,6000,2.18)
        count = 0
        print(a)
        for h in hdul:
            data = h[1].data
            try: resdata.append(spectres(new_grid,data['wave'],a[count]*data['flux'],a[count]*1/np.sqrt(data['ivar'])))
            except KeyError:  resdata.append(spectres(new_grid,data['wave'],a[count]*data['flux'],a[count]*data['sig']))
            else: pass
            count += 1

        ivar_new = np.array([i[1] for i in resdata])
        ivar_new = 1/np.square(ivar_new)
        flux_new = np.array([i[0] for i in resdata])
        new_data = np.rec.array([np.sum(flux_new*ivar_new,axis=0)/np.sum(ivar_new,axis=0),np.sqrt(1/np.sum(ivar_new,axis=0)),new_grid],
                    formats='float32,float32,float32',names='flux,sig,wave')
        bintable = fits.BinTableHDU(new_data)
        name, ok = qt.QInputDialog.getText(self,"Filename","Name the fits file:")
        dir_ = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtGui.QFileDialog.ShowDirsOnly)
        if len(dir_) == 0:
            return
        cwd = os.getcwd()
        os.chdir(dir_)        
        if not(ok):
            name = "basic"
        bintable.writeto(name + ".fits")
        os.chdir(cwd)        



    '''
    NOTE: This begins the 2D image utilities
    '''
    #TODO: likely "easy" to setup w/ photutils
    #amazing examples of analysis tools with pyqtgraph
    #launch with python -m pyqtgraph.examples
    #image analysis example gives a great tool

    #NOTE: for analysis of 2D spectra there exists a fitting algorithm with astropy (astropy.modeling.(fitting,models))
    # it seems to be good for analysis of flat frams
    # can look at ccdproc for coadding
    # image_view widget has a close function that may be preferable to closing the whole dock
    # image_view has a getROI function that returns the ROI plotWidget, this could be 
    # saved or transferred to the 1dPlotting for further analysis (fitting, wl determination, etc.)
    #TODO: 2d coadding is a feature I want to add in the near future. for pypeit I want to 
    # grab the bitmask and turn it into a list of 1's and zeros to use on any image. 
    # then we can just multiply element by element with the desired matrix to get our image
    # (this presumes that bits that are "on" are not desired and should be given a value of 0)

    #TODO: Pypeit coadding that is useful does: (Science - skymodel)*np.sqrt(ivarmodel)*(mask == 0)
    #TODO: We can extract individual slits using pypeit's spec2dobj Class
    # i.e. spec2DObj = spec2dobj.Spec2DObj.from_file(filename)
    # slitmask = spec2DObj.slits.slit_img(flexure=spec2DObj.sci_spat_flexure)
    # then the slit is spec2DObj.sciimg*(slitmask==slit_ID)
    # where the ID's are given by filename as extension name DET##-SLITS
    # grab the data and then they are given by ['spat_id']
    def pyp_coadd(self,files):
        multi = []
        choice, ok = qt.QInputDialog.getItem(self,"Which Detector?","choose a data set:",['1','2'],0,False)
        choice = int(choice)
        ID = []
        Intrument = []
        err = []
        for f in files:
            hdul = fits.open(f)
            hdr = hdul[0].header
            nums = np.array([1,3,5,9,8])
            if (hdul[1].header['EXTNAME'].find('DET01') != -1) and (choice == 2):
                nums += 12
            multi.append((hdul[nums[0]].data-hdul[nums[1]].data)*np.sqrt(hdul[nums[2]].data)*(hdul[nums[3]].data==0))
            err.append(hdul[nums[2]].data)#inverse variance
            #embed()
            try:
                ID.append(hdr['FRAMEID'])
                Intrument.append(hdr['INSTRUME'])
            except:
                pass
        #TODO: how to check if mosfire? shift using scipy.ndimage.shift(image,shift,cval=0.0). 
        #How to check if B image? INSTRUME: gives intrument name. FRAMEID: gives whether A or B
        for e,I in enumerate(Intrument):
            if I == "MOSFIRE":
                if ID[e] == "B":
                    multi[e] = sc.ndimage.shift(multi[e],np.array([0,-14]))
                    err[e] = sc.ndimage.shift(err[e],np.array([0,-14]))
                else:
                    pass

        Intrument = np.array(Intrument)
        multi = np.array(multi)
        err = np.array(err)
        if np.any(Intrument == "MOSFIRE"):
            data = np.median(multi,axis=0)
            merr = np.sqrt((np.sum(err*multi**2,axis=0)/np.sum(err,axis=0) - (np.sum(err*multi,axis=0)/np.sum(err,axis=0))**2)/(len(multi)-1))*np.sqrt(np.pi*(2*len(multi)-1)/(4*len(multi)-4))
        else: 
            '''       
            data = np.average(multi,axis=0,weights=err)#divide by zero, but not a problem if np.sum(x*w)/np.sum(w)?
            merr = np.sqrt((np.sum(err*multi**2,axis=0)/np.sum(err,axis=0) - (np.sum(err*multi,axis=0)/np.sum(err,axis=0))**2)/(len(multi)-1))
            '''
            data = np.median(multi,axis=0)
            if len(multi) > 1: merr = np.sqrt((np.sum(err*multi**2,axis=0)/np.sum(err,axis=0) - (np.sum(err*multi,axis=0)/np.sum(err,axis=0))**2)/(len(multi)-1))*np.sqrt(np.pi*(2*len(multi)-1)/(4*len(multi)-4))
            else: merr = err
        merr[np.isnan(merr)] = 1000
        self.imv = pg.ImageView(view=pg.PlotItem())
        self.plot2d.addWidget(self.imv)
        if data.ndim != 2:
            qt.QMessageBox.about(self,"Error","Data has {0} but images require 2".format(data.ndim))
            self.imv.close()
            return
        self.imv.setImage(data,levels=(-10,10),xvals=hdul[nums[4]].data)#xvals sets the wavelengths, later grabbed as self.imv.tVals
        self.imv.setCursor(QtCore.Qt.CrossCursor)
        self.isigprox = pg.SignalProxy(self.imv.scene.sigMouseMoved,rateLimit=60,slot=self.imageHoverEvent)
        self.err2d = merr

    def save_coadd(self):
        if self.imv:
            data = self.imv.image
            try:
                wl = self.imv.tVals
                tabhdu = fits.ImageHDU(wl)
            except:
                pass
            name, ok = qt.QInputDialog.getText(self,"Filename","Name the fits file:")
            if not(ok):
                return
            dir_ = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtGui.QFileDialog.ShowDirsOnly)
            if len(dir_) == 0:
                return
            cwd = os.getcwd()
            os.chdir(dir_)
            img = fits.ImageHDU(data)
            prim = fits.PrimaryHDU()
            primhead = prim.header
            primhead['EXTEND'] = True
            new_hdul = fits.HDUList([prim,img,tabhdu])
            new_hdul.writeto(name,overwrite=True)
            os.chdir(cwd)

    def Extract1d(self):
        data = self.roi.getArrayRegion(self.imv.image,self.imv.imageItem,axes=(0,1))
        wlarr = self.roi.getArrayRegion(self.imv.tVals,self.imv.imageItem,axes=(0,1))
        err = self.roi.getArrayRegion(self.err2d,self.imv.imageItem,axes=(0,1))
        p,pcov = curve_fit(gauss0,np.arange(len(data[0])),np.sum(data,axis=0),(np.max(np.sum(data,axis=0)),np.where(np.sum(data,axis=0)==np.max(np.sum(data,axis=0)))[0][0],1))
        lower = int(p[1] - 3*p[2])
        upper = int(p[1]+3*p[2])
        spec = data[:,lower:upper]
        err1d = err[:,lower:upper]
        flux = np.sum(spec,axis=1)
        ferr = np.sqrt(np.sum(err1d**2,axis=1))
        wl = wlarr[:,int(np.floor(p[1]))]
        new_data = np.rec.array([flux,ferr,wl],
                    formats='float32,float32,float32',names='flux,sig,wave')
        bintable = fits.BinTableHDU(new_data)
        dir_ = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtGui.QFileDialog.ShowDirsOnly)
        if len(dir_) == 0:
            return
        cwd = os.getcwd()
        os.chdir(dir_) 
        bintable.writeto("spec1d.fits",overwrite=True)
        os.chdir(cwd)

    def makeROI(self):
        self.roi = pg.RectROI([10,50],20,pen='r')
        self.imv.addItem(self.roi)

    def remove2D(self):
        self.imv.close()

    def show_2d_fits(self,fileName = ""):
        try: hdul = fits.open(fileName)
        except: qt.QMessageBox.about(self,"Error","Something went wrong opening {0}".format(fileName))

        header = hdul[0].header
        s = []
        #TODO: need to generalize this. EXT0 might apply to a few cases (certainly pypeit),
        # but not all (i.e. student generated fits files)
        for i in header:
            if i.find('EXT0') != -1:
                s.append(i)
        if len(s) == 0:
            s = np.arange(len(hdul))
            choice, ok = qt.QInputDialog.getItem(self,"Which extension?","choose a data set:",s.astype(str),0,False)
            choice = int(choice)
        else:
            choice, ok = qt.QInputDialog.getItem(self,"Which extension?","choose a data set:",[header[i] for i in s],0,False)        
        data = hdul[choice].data
        self.imv = pg.ImageView(view=pg.PlotItem())
        self.plot2d.addWidget(self.imv)
        if data.ndim != 2:
            qt.QMessageBox.about(self,"Error","Data has {0} dimensions but images require 2 dimensions".format(data.ndim))
            self.imv.close()
            return
        self.imv.setImage(data,levels=(-10,10))
        self.imv.setCursor(QtCore.Qt.CrossCursor)
        self.isigprox = pg.SignalProxy(self.imv.scene.sigMouseMoved,rateLimit=60,slot=self.imageHoverEvent)
    
            
if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon("Images/Logo_final.png"))#This will display the logo on macs
    #TODO: Create better logo
    ex = App()
    sys.exit(app.exec_())