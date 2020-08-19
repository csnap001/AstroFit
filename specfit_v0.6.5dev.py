"""
Created Thursday 8th, 2019 by Christopher Snapp-Kolas

TODO: Rename the file so that it better describes its use.
TODO: Should there be a toolbar?
TODO: comment all code for easier readability
TODO: add "clear screen" option for clearing displays
TODO: Add threading of functions to allow continued functionality of other QtGui objects
TODO: Add progress bar for fits: self.progress = QtGui.ProgressBar(self) -> self.progress.setValue(**some increasing number**)
TODO: Rework Gui using Pyqtgraph examples in site-packages within python, These seem more powerful

Module for GUI spectroscopic fitting environment based on specutils
and astropy. (Possibly, desired) This module will also have basic image arithmatic capabilities.

Meta data:

Class Functions:

"""

import PyQt5.QtWidgets as qt
from pyqtgraph.Qt import QtCore,QtGui
import pyqtgraph as pg
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
from My_Packages import MCMC_fine_tuning as mcmc#TODO:Consider changing to pymultinest
from My_Packages import FuncDef as fxn
from My_Packages import confidence as conf
#import matplotlib.pyplot as plt
from functools import partial
import os
import corner
from IPython import embed


class App(QtGui.QMainWindow):
    posSignal = QtCore.pyqtSignal(list) #this must be defined as a class variable
    XSignal = QtCore.pyqtSignal(tuple)
    #outside of any function in order to work
    def __init__(self):
        super().__init__()
        '''
        setting images for cursors
        '''
        Lpm = QtGui.QPixmap('Cursor_imgs/Left_under.png')
        Lpm = Lpm.scaled(30,30) #can use this to scale
        #TODO: set position along left/right vertical line
        self.Lcursor = QtGui.QCursor(Lpm) #TODO: set cursor position to corner
        Rpm = QtGui.QPixmap('Cursor_imgs/Right_under.png')
        Rpm = Rpm.scaled(30,30)
        self.Rcursor = QtGui.QCursor(Rpm)

        self.title = "SpecFit"
        self.left = 20
        self.top = 20
        self.width = 1000
        self.height = 1000

        self.plt = []#list containing plot information for each added 1d plot
        self.names1d = []#list containing names1d of 1d plots to allow for removal

        self.is1d = False
        self.is2d = False
        self.isTab = False
        self.is1dPos = False
        self.is2dPos = False

        self.flux = None
        self.err = None
        self.isMask = False
        self.Mask = []

        self.xpos = None
        self.ypos = None
        self.Lbound = None
        self.Rbound = None
        self.ximg = []
        self.yimg = []
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
        self.fit = []

        #connects signals to functions allowing for the .emit() to cause an event
        self.posSignal.connect(self.getPos)
        self.XSignal.connect(self.getX)
        

    def initUI(self):
        self.area = DockArea()
        self.setCentralWidget(self.area)
        self.resize(self.width,self.height)
        self.setWindowTitle(self.title)
        self.dtool = Dock("Toolbar",size=(1,1))
        self.dplot = Dock("plots",size=(500,300))
        self.dTable = Dock("Table",size=(500,300))
        self.area.addDock(self.dtool,'top')
        self.area.addDock(self.dplot,'bottom',self.dtool)
        self.area.addDock(self.dTable,'right',self.dplot)
        self.Gwin1d = pg.GraphicsWindow()
        self.table = pg.TableWidget()
        #plt = self.Gwin1d.addPlot(title="plotting")#NOTE: not fixed
        #plt.addLegend()#NOTE:addLegend must come before plot
        #plt.plot([1,2,3],[1,2,3],pen=(255,0,0),name='red')
        #lr = pg.LinearRegionItem([1,2])#get values using lr.getRegion()
        #print(plt.Title)
        #argument should be based off of data
        #TODO:add linear region button (have to ensure it goes to desired plot)
        #TODO:add "add plot" button to allow for displaying multiple plots at a time
        #plt.addItem(lr)
        self.dplot.addWidget(self.Gwin1d)
        self.dTable.addWidget(self.table)
        self.dtool.hideTitleBar()
        
        self.fitProgress = QtGui.QProgressBar(self)
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

    
    def addPlot(self,name):
        #TODO:need to figure out how to handle placement after removing plots
        num = len(self.plt)
        if num < 2:
            self.plt.append(self.Gwin1d.addPlot(title=name))
        elif num >= 2 and num < 4:
            self.plt.append(self.Gwin1d.addPlot(row=1,col=num-2,title=name))
        elif num >=4 and num < 6:
            self.plt.append(self.Gwin1d.addPlot(row=2,col=num-4,title=name))
        else:
            qt.QMessageBox.about(self,"Too Many!","Too many plots, not adding")
            return -1

    def removePlot(self):
        choice, ok = qt.QInputDialog.getItem(self,"Removing 1d?","Which plot?:",self.names1d,0,False)
        loc = self.names1d.index(choice)
        #embed()
        if ok:
            self.Gwin1d.removeItem(self.plt[loc])
            self.plt.pop(loc)
        else:
            qt.QMessageBox.about(self,"Done","Not Removing any plots")

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
                check = self.addPlot(name=fileName)
                if check != -1:
                    self.plot_1d_fits(fileName=fileName)
                    self.names1d.append(fileName)
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
        
        self.tabWid.clear()
        fits = fileName.find(".fits")
        dat = fileName.find(".dat")
        txt = fileName.find(".txt")
        if not(fits == -1):
            data = Table.read(fileName,format='fits')
            self.tabWid.setColumnCount(len(data.colnames))
            self.tabWid.setRowCount(len(data))
        
            self.tabWid.setHorizontalHeaderLabels(data.colnames)
            for i in range(len(data)):
                for j in range(len(data.colnames)):
                    self.tabWid.setItem(i,j,qt.QTableWidgetItem(str(data[i][j]))) #argument must be string
        if not(dat == -1):
            pass
        if not(txt == -1):
            pass

    #TODO: currently this is specifically desigened for bill's data
    # it would be better to simply plot the data and then allow the user 
    # to rearrange by selecting which datum will go along the x and y axes.
    # could also have option for plotting multiple data elements, such
    # as 'raw' and 'error' (with the ability to have both at the same time)
    # this requires knowledge about astropy headers.
    #TODO: Should also check data dimensions to ensure that the input file is not
    # an image file. Must force one dimensional spectrum data
    #NOTE: Each fits file extension (i.e. hdul[1]) as an NAXIS memeber that states the 
    # number of data axes. So in bill's data the header has NAXIS = 1 and then in the 
    # image files we have NAXIS = 2. We can use this to implement restrictions
    #TODO: Consider Plotting BINNED data, since this more accurately represents the physical 
    # measurement. must consider how to do this. Ask Brian. Use .step
    #NOTE: Wavelength data can be given by some start wavelength, delta wavelength, and number of pixels
    # Therefore, it is necessary to grab this information and construct the wavelength bins
    def plot_1d_fits(self,fileName=""):
        count = len(self.plt) - 1
        #TODO: add legend (should be done for fits as well). legend = pg.LegendItem, legend.setParentItem(self.plot_view)
        #TODO: add title as filename that has been opened and plotted
        if fileName:
            hdul = fits.open(fileName)
            stwl = hdul[0].header['CRVAL1']
            print(stwl)
            step = hdul[0].header['CDELT1']
            flux = hdul[1].data
            print(flux[0])
            wl = [stwl + step*i for i in range(len(flux))]
            err = hdul[2].data
            pen = pg.mkPen(color='b')
            self.flux = self.plt[count].plot(wl,flux,pen=pen)
            self.err = self.plt[count].plot(wl,err)
            pg.SignalProxy(self.plt[count].scene().sigMouseMoved, rateLimit=60,slot=self.mouseMoveEvent)
            self.plt[count].scene().sigMouseMoved.connect(self.mouseMoveEvent)
            
            
       

    def get1d_from2d(self):
        #likely will call plot_1d_fits. pg.affineSlice might be helpful here
        pass

    def get_wvlngth(self):
        #for wavelength solution, use pixmap. likely need get1d_from2d
        pass

    #TODO: this is currently broken, need method for adding frames to view one by one for multiple frame analysis
    # See pyqtgraph examples for help on this
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
    #TODO: Error results seem unphysical. Method for grabbing instrument uncertainties?
    def fitting_1d(self):
        #calling specutils or some MCMC routine (perhaps emcee)
        #TODO: rename to parameterEW or something along these lines to distinguish from brute force integration method
        """
        This will take in user input on the location of some peak and then fit a gaussian to that peak.
        """
        #TODO: Need to propogate errors from continuum fit

        if self.ewBut.isChecked():
            self.plot_view.setCursor(self.Lcursor)
            self.is1dPos = True
            Pos(self,self.plot_view)
        if not(self.ewBut.isChecked()):
            if(not(self.xpos == None)):
                #TODO: getting continuum choice from user. How to name these?
                # also need to remove all continua once a new plot is loaded
                data = self.plot_view.listDataItems()
                if len(data) > 0:
                    data = [data[0].getData(),data[1].getData()]
                    wl = data[0][0]
                    flux = data[0][1]
                    err = data[1][1]
                   
                    if len(self.contfit[0]) > 0:
                        #these assume that contfit is NOT an array of arrays, which it is if multiple continua exist
                        items = ("{}".format(i) for i in range(len(self.contfit[0])))
                        cont, ok = qt.QInputDialog.getItem(self,"Continuum","Which Continuum?:",items,0,False)
                        self.fit = self.contfit[0][int(cont)]
                        self.cname = self.contfit[1][int(cont)]
                        if self.cname == 'pl' and ok:
                            continuum = fxn.Pow(wl,*self.fit)
                        elif self.cname == 'l' and ok:
                            continuum = self.line(wl,*self.fit)
                        elif self.cname == 'p' and ok:
                            continuum = fxn.polynomial(wl,*self.fit)
                        else:
                            qt.QMessageBox.about(self,"Done","Exiting without fit")
                            return

                        normflux = flux/continuum
                        mask = (wl > self.ewBounds[0]) & (wl < self.ewBounds[1])
                        finalflux = normflux[mask]
                        finalwl = wl[mask]
                        finalerr = err[mask]/continuum[mask]
                        #TODO: Need to propogate errors more appropriately. 
                        # perhaps have code block that uses all values of MCMC simualtion
                        # and grab the sd from there? 
                        index = np.where(finalflux == np.max(finalflux))
                        peakWl = finalwl[index]
                        peakFl = finalflux[index]

                        

                        centB = (self.xpos - 0.05*self.xpos,self.xpos + 0.05*self.xpos)
                        self.plot_view.plot([centB[0],centB[0]],[0,self.ypos],pen = 'c')
                        self.plot_view.plot([centB[1],centB[1]],[0,self.ypos],pen='m')
                        zcent = np.round((peakWl - 1215.67)/1215.67,6) #specific to lyman alpha
                        #TODO: Need to add choice for which line is being fit, or an arbitrary position
                        #TODO: is np.round necessary?
                        zB = ((zcent - 0.2*zcent), (zcent + 0.2*zcent))#TODO:Is this too constraining?
                        zB = (zB[0][0],zB[1][0])#necessary b/c zB is created as array of arrays and numpy fails with array inputs
                        sigB = (0.001*(self.ewBounds[1] - self.ewBounds[0]), 0.2*(self.ewBounds[1]-self.ewBounds[0]))#TODO: Probably not best choice for bounds
                        ampB = (peakFl - 0.2*peakFl,peakFl + 0.2*peakFl) #TODO: is peak too tightly constrained?
                        ampB = (ampB[0][0],ampB[1][0])#same as zB

                        self.Fitter(fxn.Lyalph,data,finalflux,finalerr,finalwl,[ampB,zB,sigB],name='EW')
                    else:
                        qt.QMessageBox.about(self,"Warning!","No continua available! First fit continuum")
                else:
                    qt.QMessageBox.about(self,"No data on screen","Not fitting")
                        
        
    #TODO: Need to track and remove plots, perhaps on user request
    #TODO: Need to propagate errors to any other fits (How to do this?)
    def Cont_Fit(self):
        """
        This will ask for a wavelength range from the user. A continuum will be fit
        to the flux data with consideration of the error spectrum. Should also consider setting up
        a list of masks to the data to use for continuum fitting
        """
        #TODO: Program crashes if no points are chosen, need to have a handler for this that informs the user
        self.isMask = True
        #self.Mcount, ok = qt.QInputDialog.getInt(self,"","How many masks?")
        if self.contBut.isChecked():
            self.plot_view.setCursor(self.Lcursor)#Setting cursor image
            Pos(self,self.plot_view)
        if not(self.contBut.isChecked()):
            self.isMask = False
            self.count = 0
            self.plot_view.setCursor(QtCore.Qt.CursorShape(2))
            self.plot_view.scene().sigMouseClicked.disconnect()
            data = self.plot_view.listDataItems()
            #geting data for fitting continuum and plotting chosen points
            if(len(data) > 0):
                boolMask = [len(self.Mask[i]) == 2 for i in range(len(self.Mask))]
                if (self.Mask) and np.all(boolMask):
                    data = [data[0].getData(),data[1].getData()]
                    trange = self.plot_view.viewRange()
                    yrange = trange[1]
                    holdToRemove = []
                    fmask = []
                    for i in range(len(self.Mask)):
                        line_data = [[self.Mask[i][0],self.Mask[i][0]],[yrange[0],yrange[1]]]
                        rine_data = [[self.Mask[i][1],self.Mask[i][1]],[yrange[0],yrange[1]]]
                        b = pg.PlotDataItem(x=line_data[0],y=line_data[1],pen='b')
                        r = pg.PlotDataItem(x=rine_data[0],y=rine_data[1],pen='r')#symbol is another argument, best for scatterplotitem
                        self.plot_view.addItem(b)
                        self.plot_view.addItem(r)
                        holdToRemove.append(b) #used to remove these lines if you do not perform the fit
                        holdToRemove.append(r)
                        fmask.append((data[0][0] > self.Mask[i][0]) & (data[0][0] < self.Mask[i][1]))
                        
                else:
                    qt.QMessageBox.about(self,"No choice", "No bounds chosen or Not an even number of bounds: #bounds = {}, not fitting".format(boolMask))
                    self.Mask = []
                    #TODO: make sure that all plotted bounds are deleted in this process
                    return
                
            else:
                qt.QMessageBox.about(self,"No Data","No data to fit to!")
                self.Mask = []
                return
            if(len(fmask)>1):fmask = [fmask[i]|fmask[i+1] for i in range(len(fmask) - 1)]
            print(fmask)
            flux = data[0][1][fmask]
            wl = data[0][0][fmask]
            err = data[1][1][fmask]#TODO: Check data type of fmask
            #Seems like all that needs to be done is to change fmask -> tuple(fmask)?
            '''
            FutureWarning: Using a non-tuple sequence for 
            multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. 
            In the future this will be interpreted as an array index, `arr[np.array(seq)]`, 
            which will result either in an error or a different result.
            err = data[1][1][fmask]
            '''
   
            Ans = qt.QMessageBox.question(self,"Masking","Mask Emission and Absorption lines?",qt.QMessageBox.Yes|qt.QMessageBox.No,qt.QMessageBox.No)


            if Ans==qt.QMessageBox.Yes:
                z, ok = qt.QInputDialog.getDouble(self,"Get redshift","z:",2.0,0.0,10.0,10)
                if ok:
                    wlCIV = 1550*(1+z)
                    wlHeII = 1640*(1+z)#TODO: re-impliment using list of lines instead
                    c = 3e5 #km/s
                    v = c*(wl-wlCIV)/wlCIV
                    CIVMask = (v > -1000) & (v < 400)
                    HeIIMask = (c*(wl-wlHeII)/wlHeII > -400) & (c*(wl-wlHeII)/wlHeII < 400)
                    TotMask = CIVMask | HeIIMask
                    x = wl[TotMask]
                    y = flux[TotMask]
                    CIVscatter = pg.ScatterPlotItem(x=x,y=y,pen=pg.mkPen('g'),brush=pg.mkBrush('g'))
                    self.plot_view.addItem(CIVscatter)
                    wl = wl[~TotMask]#~ is bitwise NOT which is how we mask in python
                    flux = flux[~TotMask]
                    err = err[~TotMask]
            else:
                qt.QMessageBox.about(self,"No Mask","Continuing fit with all data in bounds")

            items = ('Power Law','Linear','Polynomial')
            func, ok = qt.QInputDialog.getItem(self,"Get function","Function: ",items,0,False)
            if func == 'Power Law' and ok:
                #TODO: Force left bound as "centroid" in Pow?
                self.Fitter(fxn.Pow,data,flux,err,wl,[(np.min(flux),np.max(flux)),(-3,1),(np.min(wl),np.max(wl))],name='continuum')
                self.cname = 'pl'
                self.contfit[1].append('pl')
            if func == 'Linear' and ok:
                self.line = partial(fxn.linear,b=self.Mask[0][0])
                self.Fitter(self.line,data,flux,err,wl,[(np.min(flux),np.max(flux)),(np.min(flux)/np.max(wl),np.max(flux)/np.min(wl))],name='continuum')
                #TODO: is the slope range reasonable?
                self.cname = 'l'
                self.contfit[1].append('l')
            if func == 'Polynomial' and ok:
                order, check = qt.QInputDialog.getInt(self,"Polynomial","What Order?:",0,False)
                if check:
                    self.cname = 'p'
                    self.contfit[1].append('p')
                    guess = []
                    guess.append((np.min(flux),np.max(flux)))
                    for i in range(order):
                        guess.append(((np.min(flux)/np.max(wl))**(i+1),(np.max(flux)/np.min(wl)**(i+1))))

                    self.Fitter(fxn.polynomial,data,flux,err,wl,guess,name='continuum')
            if not(ok):#TODO: consider using else, this is kind of weird coding
                qt.QMessageBox.about(self,"Continuum","Not fitting continuum")
                self.Mask = []
                for i in holdToRemove:
                    self.plot_view.removeItem(i)
                
                
                
            self.Mask = []
        return
        
    
    def Fitter(self,func,data,flux,err,wl,bounds,name = ''):
        '''
        helper function for fitting algorithms (cont_fit, and ew)
        '''
        #TODO: Threading for progress bar? Maybe just add in self.MCMC code block?
        self.fitProgress.setValue(50)
        mymc = mcmc.fit(func,wl,flux,err, 7000,*bounds) #was 1000
        hists = []
        ps = []
        #perc = []
        test = []
  

        self.fitProgress.setValue(75)
        for i in range(len(mymc[0])):
            #perc.append(np.percentile(mymc[0][i],[16,50,84]))
            hists.append(np.histogram(mymc[0][i],bins = 250)) 
            test.append(np.where(hists[i][0] == np.max(hists[i][0])))#TODO: Is this the best fit? should devise an algorithm that determines best fit given cross correlations
            '''
            if(np.max(hists[i][0]) >= 8*np.mean(hists[i][0])):#For distributions that are essentially "Noise" it seems best to simply take the mean value
                ps.append(hists[i][1][test[i]])
            else:
            '''
            ps.append(np.percentile(mymc[0][i],50))#WARNING: just taking mean value, may not always be best


        hists = np.array(hists)
        value = [conf.ConfInt(hists[i][0],hists[i][1][1:],0.68) for i in range(len(hists))]
        print('68 interval: {}'.format(value))
        if name == 'continuum':
            self.contfit[0].append(ps)
            self.conterr.append(value)
            if name in self.pdfs.keys():
                name = self.cname + name #adding on extra name to differentiate continuum, won't work indefinitly
            self.pdfs[name] = (mymc[0],ps)
            pen = (100,90,0) # sets the color of the line, TODO: Should make this user adjustable
            cont = 1.0
        if name == 'EW':
            self.ewfit.append(ps)
            self.ewferr.append(value)
            self.pdfs[name] = (mymc[0],ps)
            ewPdf = np.sqrt(2*np.pi)*np.array(mymc[0][0])*np.array(mymc[0][2])
            ewMeas = np.sqrt(2*np.pi)*ps[0]*ps[2]
            self.pdfs['EWPDF'] = [ewPdf,[ewMeas]]
            pen = (0,100,0)
            print("result of EW fit {}".format(self.ewfit))

            #This is used for GUI image such that we can see the fit
            #TODO: the whole self.fit train of thought is very kludgy, need to come up with cleaner method
            if self.cname == 'pl':
                cont = fxn.Pow(data[0][0],*self.fit)
            elif self.cname == 'l':
                cont = self.line(data[0][0],*self.fit)
            elif self.cname == 'p':
                cont = fxn.polynomial(data[0][0],*self.fit)
        self.fitProgress.setValue(100)    
        self.plot_view.plot(data[0][0],cont*func(data[0][0],*ps),pen=pen)

    def showPDFS(self):
        '''
        Helper function for displaying corner plot of fit parameters.
        This will display the results of a given fit.
        TODO: Need to make the display more informative, right now it is unclear
        what function we are looking at and which parameter goes where.
        TODO: self.pdfs should hold all non-deleted fits (that is, it should have all continua fits)
        '''
        func, ok = qt.QInputDialog.getItem(self,"Get Fit","Which Fit?: ",self.pdfs.keys(),0,False)

        if ok:
            labels = ['a{}'.format(i) for i in range(len(self.pdfs[func][0]))]
            corner.corner(np.transpose(self.pdfs[func][0]),bins=250,quantiles=[0.16,0.5,0.84],show_titles=True,
                        labels=labels,verbose=True,plot_contours=True,title_fmt=".2E",truths=self.pdfs[func][1],
                        levels=(0.68,)) #must be (values, # of parameters), (i.e. (365,4) corresponds to a fit with four parameters)
            #plt.show()
                
    def plotColoring(self,isFlux = False, isErr=False):
        color = qt.QColorDialog.getColor()
        if not(self.flux == None) and isFlux:
            self.flux.setPen(color)
        if not(self.err==None) and isErr:
            self.err.setPen(color)
    
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

    count = 0    
    def onClick(self,event):
        pos = event.scenePos()
        vb = self.plot_view.getViewBox()
        if self.is1dPos:
            
            
            if self.count == 0:
                self.plot_view.setCursor(self.Rcursor)
                self.Lbound = vb.mapSceneToView(pos).x()
                
                
            if self.count == 1:
                self.plot_view.setCursor(QtCore.Qt.CursorShape(2))
                self.Rbound = vb.mapSceneToView(pos).x()
                
            if self.count == 2:
                self.plot_view.scene().sigMouseClicked.disconnect()
                self.ewBut.toggle()
                self.is1dPos = False #this is necessary b/c otherwise this piece of
                #code will be called again when pressing other buttons since it will
                #still be true
                x = vb.mapSceneToView(pos).x()#grabs x position relative to axes
                y = vb.mapSceneToView(pos).y()#grabs y position relative to axes
                self.posSignal.emit([x,y,self.Lbound,self.Rbound])
                self.count = 0
                return
            self.count += 1
        if self.is2dPos:
            pass
            #do something with 2d
        if self.isMask:#Continuum bounding
            #TODO: if multiple left and right bounds are chosen and then
            # the fitting is canceled, only the most recent plot is removed
            # need to remove all masks that were chosen. TODO additionally, need
            # a check for whether the left bound is actually to the left
            # of the right bound
            # TODO: Update fitting region choice by using LinearRegionItem([])?
            if self.count == 0:
                self.plot_view.setCursor(self.Rcursor)
                self.xpos = vb.mapSceneToView(pos).x()
            if self.count == 1:
                self.plot_view.setCursor(self.Lcursor)
                self.ypos = vb.mapSceneToView(pos).x() #TODO: is ypos named correctly?
                self.XSignal.emit((self.xpos,self.ypos))
                self.count = 0
                return
            self.count += 1
           

        
    @QtCore.pyqtSlot(list)
    def getPos(self,value):
       
        self.xpos = value[0]
        self.ypos = value[1]
        print(self.xpos,self.ypos)
        self.ewBounds = (value[2],value[3])
        self.fitting_1d()
        

    @QtCore.pyqtSlot(tuple)
    def getX(self,value):
        self.Mask.append(value)
        
        


            



if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())