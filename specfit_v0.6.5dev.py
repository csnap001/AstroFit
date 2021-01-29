"""
Created Thursday 8th, 2019 by Christopher Snapp-Kolas

TODO: Rename the file so that it better describes its use.
TODO: add "clear screen" option for clearing displays
TODO: Add threading of functions to allow continued functionality of other QtGui objects
TODO: Add progress bar for fits: self.progress = QtGui.ProgressBar(self) -> self.progress.setValue(**some increasing number**)
This is a difficult task, as pymc3 holds all this info internally and doesn't appear to allow for access
TODO: 2 primary functions need to be built. 1: need non-parameterized equivalent width determination
2: Need line center determination. use paper that fred sent (find the first moment)
generally speaking there is more to be done than this, but this will be sufficient for now
TODO: Re-organize storage of data. There should be meta-data fxns that feed to data analysis fxns
Should write out and organize my thoughts on this. 
TODO: Should save data to files as well (.dat?), currently using pickle, works for me but not generally what people
will want
TODO: Add operations to data elements (mostly to account for inverse variance stuff)

TODO: Add in Stellar population synthesis model plotting from Hoki (which is a wrapper to BPASS)
This likely needs more work with data containers to help out with the workflow of the code

Module for GUI spectroscopic fitting environment based on pymc3
and astropy. (Possibly, desired) This module will also have basic image arithmatic capabilities.

Meta data:

Class Functions:

"""

import PyQt5.QtWidgets as qt
from pyqtgraph.Qt import QtCore,QtGui
import pyqtgraph as pg
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.modeling import models
from pyqtgraph.dockarea import *
import numpy as np
import sys
import create_buttons as cb
import connect_buttons as conb
import Initialize_View_Widgets as IVW
from Layouts import Lay
from Get_Pix_Position import Pos, off
import time
from My_Packages import FuncDef as fxn
import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
from functools import partial
import os
import corner
import scipy as sc
from IPython import embed
from astropy import cosmology
import pickle

def conf_low(data,half=0.34):
    mean = np.mean(data)
    sd = np.std(data)
    kde = sc.stats.gaussian_kde(data,0.05)
    i = 0
    conf = 0
    while conf < half:
        i += 1
        conf = kde.integrate_box_1d(mean-i*0.01*sd,mean)
    return mean - i*0.01*sd

def conf_high(data,half=0.34):
    mean = np.mean(data)
    sd = np.std(data)
    kde = sc.stats.gaussian_kde(data,0.05)
    i = 0
    conf = 0
    while conf < half:
        i += 1
        conf = kde.integrate_box_1d(mean,mean+i*0.01*sd)
    return mean + i*0.01*sd

class App(QtGui.QMainWindow):

    def __init__(self):
        super().__init__()

        self.title = "SpecFit"
        self.left = 20
        self.top = 20
        self.width = 1000
        self.height = 1000

        self.plt = []#list containing plot information for each added 1d plot
        self.names1d = []#list containing names1d of 1d plots to allow for removal
        self.lrs = []#list of linear regions for data selection (i.e. for fitting continua or EWs)
        self.regPlot = []#list to hold region plots
        self.yrange = [0,1]
        self.prior = []#list of continuum priors that can be used
        self.fill = None

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


    def initUI(self):
        
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.resize(self.width,self.height)
        self.setWindowTitle(self.title)
        self.dtool = Dock("Toolbar",size=(1,1))
        self.dplot = Dock("plots",size=(500,300))
        self.dTable = Dock("Table",size=(500,300))
        self.regDock = Dock("Regions",size=(500,300))
        self.regWin = pg.GraphicsWindow()
        self.area.addDock(self.dtool,'top')
        self.area.addDock(self.dplot,'bottom',self.dtool)
        self.area.addDock(self.regDock,'below',self.dplot)
        self.area.addDock(self.dTable,'right',self.dplot)
        self.Gwin1d = pg.GraphicsWindow()
        self.Gwin1d.resize(1000,600)
        self.table = pg.TableWidget()
        self.dplot.addWidget(self.Gwin1d)
        self.regDock.addWidget(self.regWin)
        self.dTable.addWidget(self.table)
        self.dtool.hideTitleBar()
        self.Gwin1d.setBackground('w')
        #Creating buttons
        cb.buttons(self)

        #Connecting functions to buttons, that is setting up the action a button does
        conb.connect(self)

        #Create 1d, 2d, buttons, and table widgets
        '''
        pg.setConfigOption('background','w')
        pg.setConfigOption('foreground','k')
        IVW.Views(self)
        #self.plot_view.setMouseTracking(True)
        #pg.SignalProxy(self.plot_view.scene().sigMouseMoved, rateLimit=60,slot=self.mouseMoveEvent)
        #self.plot_view.scene().sigMouseMoved.connect(self.mouseMoveEvent)
        '''
        Lay(self)
        self.show()

    #TODO: need to be able to update each region not just the most recent
    def updateLR(self):
        #TODO: The y-range for each plot is only connected to the most recent region, why?
        for i in range(len(self.lrs)):
            if self.lrs[i].sigRegionChanged:
                xlen = self.lrs[i].getRegion()
                self.lrs[i].setRegion(xlen)
                abdata = self.regPlot[i].listDataItems()
                data = (abdata[0].getData()[0],abdata[0].getData()[1],abdata[1].getData()[1])
                mask = (data[0] > self.regPlot[i].getViewBox().viewRange()[0][0]) & (data[0] < self.regPlot[i].getViewBox().viewRange()[0][1])
                if len(data[1][mask]) > 0: self.yrange = [np.min(data[1][mask]) - 0.8*np.min(data[1][mask]),1.5*np.max(data[1][mask])]

    def addRegionPlot(self,data):
        num = len(self.regPlot)
        if num < 2:
            self.regPlot.append(self.regWin.addPlot(title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux')
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error')
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)
        elif num >= 2 and num < 4:
            self.regPlot.append(self.regWin.addPlot(row=1,col=num-2,title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux')
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error')
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)
        elif num >= 4 and num < 6:
            self.regPlot.append(self.regWin.addPlot(row=2,col=num-4,title="Region data"))
            self.regPlot[len(self.regPlot)-1].addLegend()
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[1],color='w',name='Flux')
            self.regPlot[len(self.regPlot)-1].plot(data[0],data[2],color='red',name='Error')
            self.regPlot[len(self.regPlot)-1].sigXRangeChanged.connect(self.updateLR)

    def addPlot(self,name):
        #TODO:need to figure out how to handle placement after removing plots
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

    def set1dFluxCol(self):
        isFlux = True
        self.plotColoring(isFlux=isFlux)

    def set1dErrCol(self):
        isErr = True
        self.plotColoring(isErr=isErr)

    def fileTab(self):
        self.isTab = True
        self.openFileNameDialog(isTab = self.isTab)

    #clear
    def openFileNameDialog(self,is1d = False ,is2d = False,isTab = False):
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        fileName, _ = qt.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","","All Files (*);;Fits Files (*.fits);;Python files (*.py)",options = options)
        name,exten = os.path.splitext(fileName) #Used to place constraints on filetype, exten grabs file extension

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
        else:
            qt.QMessageBox.about(self,"Opening File","ERROR: (Program Name) only supports .fits files")
            #TODO: Need to come up with name for GUI application

    def openFileNamesDialog(self):
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        files, _ = qt.QFileDialog.getOpenFileNames(self, "QFileDialog.getOpenFileNames()", "","All Files (*);;Fits Files (*.fits);;Python files (*.py)", options=options)
        if files:
            print(files)
    
    def table_create(self,fileName=""):
        
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

    #NOTE: Each fits file extension (i.e. hdul[1]) as an NAXIS memeber that states the 
    # number of data axes. So in bill's data the header has NAXIS = 1 and then in the 
    # image files we have NAXIS = 2. We can use this to implement restrictions
    #TODO: Consider Plotting BINNED data, can use stepMode = True, but requires len(x) == len(y) + 1
    def updatePlot(self,count,wl,flux,err):
        self.plt[count].clear()
        pen = pg.mkPen(color='b')
        self.plt[count].addLegend()
        self.flux.append(self.plt[count].plot(wl,flux,pen=pen,name='Flux'))
        self.err.append(self.plt[count].plot(wl,err,name='Error'))
    
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
                E, ok = qt.QInputDialog.getItem(self,"Which spectrum?","choose extension",[str(i) for i in np.arange(len(data)+1)],0,False)
                E = int(E)
                x, ok = qt.QInputDialog.getItem(self,"data","Choose your x axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
                y, ok = qt.QInputDialog.getItem(self,"data","Choose your y axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
                err, ok = qt.QInputDialog.getItem(self,"data","Choose your err axis:",[data[E].columns[i].name for i in range(len(data[E].columns))],0,False)
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

        while True:
            name, ok = qt.QInputDialog.getItem(self,"data","Which data element?:",data_names.keys(),0,False)
            if not(ok):
                qt.QMessageBox.about(self,"No operation?","No data element was chosen")
            op, ok = qt.QInputDialog.getItem(self,"what Math?","What operation would you like to perform",["invert","square root"],0,False)
            if op == "invert":
                data_names[name] = 1/data_names[name]
            if op == "square root":
                data_names[name] = np.sqrt(data_names[name])

            self.plt[dat_choice].clear()
            self.flux = self.plt[dat_choice].plot(data_names['wl'],data_names['flux'],pen='b')
            self.err = self.plt[dat_choice].plot(data_names['wl'],data_names['err'])

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
           
            if len(self.contfit[0]) > 0:
                items = ("{}".format(i) for i in range(len(self.contfit[0])))
                cont, ok = qt.QInputDialog.getItem(self,"Continuum","Which Continuum?:",self.pdfs.keys(),0,False)
                if ok:
                    #NOTE: this should grab the relevant Probability density functions
                    samples = self.pdfs[cont][0]
                else:
                    qt.QMessageBox.about(self,"Done","Exiting without fit")
                    return
                mask = (wl > lr[0]) & (wl < lr[1])
                finalflux = flux[mask]
                finalwl = wl[mask]
                finalerr = err[mask]
                ind = np.where(finalflux == np.max(finalflux))
                maxwl = finalwl[ind]
                cont_pdf = fxn.Pow(maxwl,samples['amp'],samples['alpha'],samples['unity'])
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
            else:
                qt.QMessageBox.about(self,"No continuum","No continuum available, not measuring")
        else:
            qt.QMessageBox.about(self,"No data","No data to measure")

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
           
            #TODO: in order to get full uncertainty of given quantity, need to define
            #prior probability as the distribution of the continuum 
            # (which has its own prior of course)
            mask = (wl > lr[0]) & (wl < lr[1])
            finalflux = flux[mask]
            finalwl = wl[mask]
            finalerr = err[mask]

            #Start Line center determination
            Lin = np.trapz(finalflux*finalwl/finalerr**2,x=finalwl)/np.trapz(finalflux/finalerr**2,x=finalwl)
            Linerr = np.trapz(finalflux*finalwl**2/finalerr**2,x=finalwl)/np.trapz(finalflux/finalerr**2,x=finalwl)-Lin**2
            F_nu = ((np.sum(1/finalerr**2))**2 - np.sum(1/finalerr**4))/(np.sum(1/finalerr**2))**2
            Linerr = np.sqrt(Linerr/F_nu)
            infline = pg.InfiniteLine(Lin,pen=(100,50,200))
            self.plt[dat_choice].addItem(infline)
            qt.QMessageBox.about(self,"Measured","Line center: {0} \xb1 {1}".format(Lin,Linerr))
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
            wl = data[0].getData()[0]
            flux = data[0].getData()[1]
            err = data[1].getData()[1]

            Ans = qt.QMessageBox.question(self,"Masking","Mask Emission and Absorption lines?",qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)

            if Ans==qt.QMessageBox.Yes:
                z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
                if ok:
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
                    add1 = (wl > 4720) & (wl < 4750)
                    add2 = (wl > 4908) & (wl < 4933)
                    TotMask = CIVMask | HeIIMask | OIIIMask | CIIIMask | OIVMask | add1 | add2
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
            #altmask = (wl > lr[0]) & (wl < 4000)
            finalwl = wl[mask] 
            index = np.where(flux[mask] == np.max(flux[mask]))
            peakWl = finalwl[index]
            #peakWl = wl[altmask][index]
            peakFl = flux[mask][index]
            #peakFl = flux[altmask][index]
            
            zB = ((peakWl - 0.1*peakWl), (peakWl + 0.1*peakWl))
            zB = (zB[0][0],zB[1][0])#necessary b/c zB is created as array of arrays and numpy fails with array inputs
            sigB = (0.01, 15)
            ampB = (0,4*peakFl)
            ampB = (ampB[0],ampB[1][0])#same as z
            #embed()
            self.Fitter(fxn.Powpgauss,data,flux[mask],err[mask],finalwl,[ampB,zB,sigB,(0,np.max(flux[mask])),(-5,5),(lr[0],lr[1])],name='EW',plt_name=dat_choice)

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
        mask = (data[0].getData()[0] > lr[0]) & (data[0].getData()[0] < lr[1])
        flux = data[0].getData()[1][mask]
        wl = data[0].getData()[0][mask]
        err = data[1].getData()[1][mask]
        
        Ans = qt.QMessageBox.question(self,"Masking","Mask Emission and Absorption lines?",qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)

        if Ans==qt.QMessageBox.Yes:
            z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
            if ok:
                wlCIV = 1550*(1+z)
                wlHeII = 1640*(1+z)#TODO: re-impliment using list of lines instead
                wlOIII = 1665*(1+z)
                wlCIII = 1908*(1+z)
                c = 3e5 #km/s
                v = c*(wl-wlCIV)/wlCIV
                CIVMask = (v > -1000) & (v < 400)
                HeIIMask = (c*(wl-wlHeII)/wlHeII > -400) & (c*(wl-wlHeII)/wlHeII < 400)
                OIIIMask = (c*(wl-wlOIII)/wlOIII > -400) & (c*(wl-wlOIII)/wlOIII < 400)
                CIIIMask = (c*(wl-wlCIII)/wlCIII > -400) & (c*(wl-wlCIII)/wlCIII < 400)
                #add1 = (wl > 4819) & (wl < 4981)
                #add2 = (wl > 5162) & (wl < 5305)
                TotMask = CIVMask | HeIIMask | OIIIMask | CIIIMask #| add1 | add2
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
            self.Fitter(fxn.Pow,data,flux,err,wl,[(0,np.max(flux)),(-3,3),(np.min(wl),np.max(wl))],name='continuum',plt_name=dat_choice,cname=func)
            self.cname = 'pl'
            self.contfit[1].append('pl')
        elif func == 'Linear' and ok:
            self.line = partial(fxn.linear,b=lr[0])
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
                self.Fitter(fxn.polynomial,data,flux,err,wl,guess,name='continuum',plt_name=dat_choice,cname=func,order=order)
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
                    step = pm.HamiltonianMC()
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
                
                amp = pm.Normal("amp",mu=(bounds[0][0]+bounds[0][1])/2,sigma=0.8*(bounds[0][1] - bounds[0][0]),testval=bounds[0][1]/2)
                centroid = pm.Normal("centroid",mu=(bounds[1][0]+bounds[1][1])/2,sigma=0.8*(bounds[1][1] - bounds[1][0]))
                sigma = pm.TruncatedNormal("sigma",mu=(bounds[2][0]+bounds[2][1])/2,sigma=0.4*(bounds[2][1] - bounds[2][0]),testval=(bounds[2][0]+bounds[2][1])/2,lower=0)
                cont_amp = pm.TruncatedNormal("cont_amp",mu=(bounds[3][0]+bounds[3][1])/2,sigma=0.8*(bounds[3][1] - bounds[3][0]),testval=(bounds[3][0]+bounds[3][1])/2,lower=0.0000001)
                alpha = pm.TruncatedNormal("alpha",mu=(bounds[4][0]+bounds[4][1])/2,sigma=0.8*(bounds[4][1] - bounds[4][0]),testval=(bounds[4][0]+bounds[4][1])/2,lower=-5,upper=5)
                unity = pm.TruncatedNormal("unity",mu=(bounds[5][0]+bounds[5][1])/2,sigma=0.2*(bounds[5][1] - bounds[5][0]),testval=(bounds[5][0]+bounds[5][1])/2,lower=leftun,upper=rghtun)

                #embed()
                #TODO: consider using non-centered reparameterization, i.e. amp = mu + sigma*amp_0, where amp_0 ~ N(0,1)
                #use 540 as example case
                step = pm.HamiltonianMC()

                mu = func(wl.astype(np.float32),amp,centroid,sigma,cont_amp,alpha,unity)
            print(basic_model.test_point)
            print(basic_model.check_test_point())
            q0 = step._logp_dlogp_func.dict_to_array(basic_model.test_point)
            p0 = step.potential.random()
            start = step.integrator.compute_state(q0,p0)
            print(start.energy)
            logp,dlogp = step.integrator._logp_dlogp_func(q0)
            print(logp)
            print(dlogp)
            #embed()
            #Likelihood of sampling distribution
            Y_obs = pm.Normal("Y_obs",mu=mu,sigma=err.astype(np.float32),observed=flux.astype(np.float32))
            
            #trace = pm.sample(20000,tune=5000,cores=6,init='adapt_diag',step=pm.step_methods.Metropolis())#used for testing parameter space
            trace = pm.sample(10000,tune=5000,target_accept=0.80,cores=10)
            #NOTE: vals['mean'].keys() gives the parameter names
            if cname == "Power Law":
                vals = az.summary(trace,round_to=10,var_names=['amp','alpha','unity'])
            elif name == "EW":
                vals = az.summary(trace,round_to=10,var_names=['amp','centroid','sigma','cont_amp','alpha','unity'])
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
            ewPdf = np.sqrt(2*np.pi)*np.array(samples['amp'])*np.array(samples['sigma'])/fxn.Pow(vals['mean']['centroid'],samples['cont_amp'],samples['alpha'],samples['unity'])
            ewMeas = np.sqrt(2*np.pi)*vals['mean']['amp']*vals['mean']['sigma']/fxn.Pow(vals['mean']['centroid'],vals['mean']['cont_amp'],vals['mean']['alpha'],vals['mean']['unity'])
            self.pdfs['EWPDF'] = [ewPdf,pd.core.series.Series([ewMeas],index=["EW"])]
            pen = (0,100,0)
            cont = 1.0
            #This is used for GUI image such that we can see the fit 
        self.arviz[name] = trace 
        self.plt[plt_name].plot(data[0].getData()[0],cont*func(data[0].getData()[0],*vals['mean']),pen=pen)
        del(basic_model)
        
    #TODO: Should consider allowing user to adjust parameterization
    #but should make this a visual thing like a slider bar (It would be helpful to output the sigma
    # compared with the original choice for sigma in the normal distribution)
    def save_data(self):
        '''
        routine to save parameter fits to .fits file
        '''
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
        cwd = os.getcwd()
        os.chdir(dir_)
        with open(name,'wb') as output:
            pickle.dump(dat_to_save,output,pickle.HIGHEST_PROTOCOL)
        os.chdir(cwd)
        #TODO: only saves to specGUI. need to figure out how to choose save location
                
    def arviz_density(self):
        """
        show density plot
        """
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        az.plot_density(self.arviz[choice])
        plt.show()

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
        choice, ok = qt.QInputDialog.getItem(self,"Which run?","choose a data set:",self.arviz.keys(),0,False)
        if not(ok):
            return
        embed()
        vals = az.summary(self.arviz[choice],round_to=5,stat_funcs=[conf_low,conf_high])#var_names?
        axes = az.plot_pair(self.arviz[choice],kind=['hexbin',"kde"],kde_kwargs={"fill_last":False},marginals=True,point_estimate="mean",marginal_kwargs={"quantiles":[0.16,0.5,0.84]})
        for i in range(len(axes)):
            axes[i][i].set_title(r'${0}={1}_{{{2}}}^{{{3}}}$'.format(vals['mean'].keys()[i],vals['mean'][i],vals['conf_low'][i]-vals['mean'][i],vals['conf_high'][i]-vals['mean'][i]))
        plt.show()

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
        self.plt[dat_choice].plot(wl,converted_meas,pen='b')
        self.plt[dat_choice].plot(wl,conv_meas_errs)


        return converted_meas, conv_meas_errs

    def showPDFS(self):
        '''
        Helper function for displaying corner plot of fit parameters.
        This will display the results of a given fit.
        TODO: self.pdfs should hold all non-deleted fits (that is, it should have all continua fits)
        '''
        func, ok = qt.QInputDialog.getItem(self,"Get Fit","Which Fit?: ",self.pdfs.keys(),0,False)

        if ok:
            labels = self.pdfs[func][1].keys()
            corner.corner(self.pdfs[func][0],quantiles=[0.16,0.5,0.84],show_titles=True,
                        labels=labels,verbose=True,plot_contours=True,title_fmt=".3E",truths=self.pdfs[func][1],
                        levels=(0.68,)) #must be (values, # of parameters), (i.e. (365,4) corresponds to a fit with four parameters)
            plt.show()
                
    def plotColoring(self,isFlux = False, isErr=False):
        choice, ok = qt.QInputDialog.getItem(self,"Which to color?","Choose data set:",self.names1d,0,False)
        if not(ok):
            qt.QMessageBox.about(self,"Not Coloring","Chose no data, not Coloring.")
            return
        dat_choice = self.names1d.index(choice)
        color = qt.QColorDialog.getColor()
        if not(len(self.flux) == 0) and isFlux:
            self.flux[dat_choice].setPen(color)
        if not(len(self.err) == 0) and isErr:
            self.err[dat_choice].setPen(color)

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
        mask = (data[0].getData()[0] > lr[0]) & (data[0].getData()[0] < lr[1])
        origFlux = data[0].getData()[1]
        origWl = data[0].getData()[0]
        flux = data[0].getData()[1][mask]
        wl = data[0].getData()[0][mask]
        err = data[1].getData()[1][mask]

        sig, ok = qt.QInputDialog.getDouble(self,"Gaussian Width", "What is the width of the gaussians?:",0.0,0,50.0,5)
        if not(ok):
            qt.QMessageBox.about(self,"Not Adding","Do not currently support no sigma choice")

        
        done = "False"
        run = 0
        EW = []
        success = []
        FWHM = 2*np.sqrt(2*np.log(2))*sig
        while not(done == "True"):
            y = origFlux[mask]
            data[0].setData(origWl,origFlux,pen='b')
            pause, ok = qt.QInputDialog.getInt(self,'pausing','Data Reset',0,0,0,0)
            run += 1
            x_0 = np.random.uniform(lr[0],lr[1])
            ew = np.random.uniform(1,20)
            EW.append(ew)
            mask_x = (wl > x_0 - FWHM/2) & (wl < x_0 + FWHM/2)
            F_cont = np.mean(y[mask_x])
            A = ew*F_cont/(np.sqrt(2*np.pi)*sig)
            gauss = fxn.gauss0(wl,A,x_0,np.sqrt(2)*sig)
            y += gauss
            data[0].setData(wl,y)
            win,ok = qt.QInputDialog.getInt(self,"EW: {}".format(ew),"did you find it?(1=yes,0=no): ",0,0,1,1)
            success.append(win)
            done,ok = qt.QInputDialog.getItem(self,"Total runs: {}".format(run),"Done?:",["True","False"],0,False)
        data[0].setData(origWl,origFlux,pen='b')
        bool_list = list(map(bool,success))
        EW = np.array(EW)
        qt.QMessageBox.about(self,"Finished","You identified {0} out of {1} gaussians with mininum found {2}".format(np.sum(success),len(success),np.min(EW[bool_list])))
        #EW = sqrt(2pi)*A*sigma/Cont(x_0)
    '''
    NOTE: This begins the 2D image utilities
    currently these are broken and require setup
    Not particularly concerned with this currently
    '''
    #TODO: likely "easy" to setup w/ photutils
    def combine_img_ext(self):
        """
        This will take the images, remove the overscan section, and then combine to 
        have a single image which will then be displayed. prefereably it will show
        the previous extensions on a seperate screen for comparison. This will require
        the use of astropy and numpy for image mathematics (try to avoid explicit for loops
        as these will eat up memory)
        TODO: Look further into pyqtgraph examples. they are very powerful.
        """
        pass

    def imgAdd(self):
        pass
    
    def imgSub(self):
        pass

    def imgMult(self):
        pass

    def imgDiv(self):
        pass

    def get1d_from2d(self):
        #likely will call plot_1d_fits. pg.affineSlice might be helpful here
        pass

    def get_wvlngth(self):
        #for wavelength solution, use pixmap. likely need get1d_from2d
        pass

    #TODO: this is currently broken, need method for adding frames to view one by one for multiple frame analysis
    # See pyqtgraph examples for help on this. Can easily accomplish this with current structure, but not interesting now
    # using docking
    #NOTE: for analysis of 2D spectra there exists a fitting algorithm with astropy (astropy.modeling.(fitting,models))
    # it seems to be good for analysis of flat frams
    def show_2d_fits(self,fileName = ""):
        text = ""
        self.view = self.img_view.addViewBox()
        while not(text == "Y" or text == "N"):
            text, _ = qt.QInputDialog.getText(self,"Get text","Multiframe? (Y,N)",qt.QLineEdit.Normal,"")
            #TODO: How to add to graphicsLayoutWidget?
            if text == 'Y':
                if fileName:
                    isFits = fileName.find('.fits')
                    data = []
                    if not(isFits == -1):
                        f = fits.open(fileName)
                        for i in range(len(f)-1):
                            data.append(f[i+1].data)
                            img = pg.ImageView(view=pg.PlotItem())
                            print(img)
                            img.setImage(data[i])
                            imgt = pg.ImageItem()
                            lay = qt.QVBoxLayout()
                            lay.addWidget(img)
                            imgt.setImage(data[i])
                            img3 = pg.makeARGB(data[i])
                            img4 = pg.makeQImage(data[i])
                            self.img_view.addItem(qt.QPushButton('hello',self))
                        f.close()
            elif text == 'N':
                #This is single frame so add single img
                if fileName:
                    isFits = fileName.find('.fits')
                    self.img_view.clear()
                    if not(isFits == -1):
                        f = fits.open(fileName)
                        img = f[1].data
                        self.img_view.setImage(img)
                        f.close()
            else:
                msg = qt.QMessageBox()
                msg.setText('Please choose Y or N')
                told = msg.exec_()
    
            
if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())