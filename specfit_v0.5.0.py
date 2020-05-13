"""
Created Thursday 8th, 2019 by Christopher Snapp-Kolas

TODO: Rename the file so that it better describes its use.
TODO: Should there be a toolbar?
TODO: comment all code for easier readability
TODO: add "clear screen" option for clearing displays
TODO: Add threading of functions to allow continued functionality of other gui objects

Module for GUI spectroscopic fitting environment based on specutils
and astropy. (Possibly, desired) This module will also have basic image arithmatic capabilities.

Meta data:

Class Functions:

"""

import PyQt5.QtWidgets as qt
import PyQt5.QtCore as core
import PyQt5.QtGui as gui
import pyqtgraph as pg
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.modeling import models
import numpy as np
import sys
import create_buttons as cb
import connect_buttons as conb
import Initialize_View_Widgets as IVW
from Layouts import Lay
from Get_Pix_Position import Pos, off
import time
from My_Packages import MCMC_fine_tuning as mcmc
from My_Packages import FuncDef as fxn
from My_Packages import confidence as conf
import matplotlib.pyplot as plt
from functools import partial
import os
import corner


class App(qt.QWidget):
    posSignal = core.pyqtSignal(list) #this must be defined as a class variable
    XSignal = core.pyqtSignal(tuple)
    #outside of any function in order to work
    def __init__(self):
        super().__init__()
        '''
        setting images for cursors
        '''
        Lpm = gui.QPixmap('Cursor_imgs/Left_under.png')
        Lpm = Lpm.scaled(30,30) #can use this to scale
        #TODO: set position along left/right vertical line
        self.Lcursor = gui.QCursor(Lpm) #TODO: set cursor position to corner
        Rpm = gui.QPixmap('Cursor_imgs/Right_under.png')
        Rpm = Rpm.scaled(30,30)
        self.Rcursor = gui.QCursor(Rpm)

        self.title = "SpecFit"
        self.left = 20
        self.top = 20
        self.width = 1000
        self.height = 1000

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
        self.contfit = [] #A list of lists holding the continuum fit parameters
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

        #connects signals to functions allowing for the .emit() to cause an event
        self.posSignal.connect(self.getPos)
        self.XSignal.connect(self.getX)
        

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top,self.width,self.height)

        
        #Creating buttons
        cb.buttons(self)

        #Connecting functions to buttons, that is setting up the action a button does
        conb.connect(self)

        #Create 1d, 2d, buttons, and table widgets
        pg.setConfigOption('background','w')
        pg.setConfigOption('foreground','k')
        IVW.Views(self)

        Lay(self)

        self.show()


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


    def openFileNameDialog(self,is1d = False ,is2d = False,isTab = False):
        options = qt.QFileDialog.Options()
        options |= qt.QFileDialog.DontUseNativeDialog
        fileName, _ = qt.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()","","All Files (*);;Fits Files (*.fits);;Python files (*.py)",options = options)
        name,exten = os.path.splitext(fileName) #Used to place constraints on filetype, exten grabs file extension

        if exten == ".fits":
            if is1d:
                self.plot_1d_fits(fileName=fileName)
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
    def plot_1d_fits(self,fileName=""):
        self.plot_view.clear()
        if fileName:
            hdul = fits.open(fileName)
            stwl = hdul[0].header['CRVAL1']
            step = hdul[0].header['CDELT1']
            flux = hdul[1].data
            wl = [stwl + step*i for i in range(len(flux))]
            err = hdul[2].data
            pen = pg.mkPen(color='k')
            self.flux = self.plot_view.plot(wl,flux,pen=pen)
            self.err = self.plot_view.plot(wl,err)
            
            
       

    def get1d_from2d(self):
        #likely will call plot_1d_fits. pg.affineSlice might be helpful here
        pass

    def get_wvlngth(self):
        #for wavelength solution, use pixmap. likely need get1d_from2d
        pass

    #TODO: this is currently broken, need method for adding frames to view one by one for multiple frame analysis
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

    def fitting_1d(self):
        #calling specutils or some MCMC routine (perhaps emcee)
        #TODO: rename to parameterEW or something along these lines to distinguish from brut force integration method
        """
        This will take in user input on the location of some peak and then fit a gaussian to that peak.
        We should use specutils here.
        """

        if self.ewBut.isChecked():
            self.plot_view.setCursor(self.Lcursor)
            self.is1dPos = True
            Pos(self,self.plot_view)
        if not(self.ewBut.isChecked()):
            if(not(self.xpos == None)):
                #TODO: grab continuum of choice from user or just use last continuum?
                # currrently not doing either. crashes if multiple continua exist
                data = self.plot_view.listDataItems()
                if len(data) > 0:
                    data = [data[0].getData(),data[1].getData()]
                    wl = data[0][0]
                    flux = data[0][1]
                    err = data[1][1]
                   
                    if len(self.contfit) > 0:
                        #these assume that contfit is NOT an array of arrays, which it is if multiple continua exist
                        if self.cname == 'pl':
                            continuum = fxn.Pow(wl,*self.contfit[0])
                        elif self.cname == 'l':
                            continuum = self.line(wl,*self.contfit[0])
                        elif self.cname == 'p':
                            continuum = fxn.polynomial(wl,*self.contfit[0])
                        normflux = flux/continuum
                        mask = (wl > self.ewBounds[0]) & (wl < self.ewBounds[1])
                        finalflux = normflux[mask]
                        finalwl = wl[mask]
                        finalerr = err[mask]
                        plt.plot(finalwl,finalflux)
                        plt.show()
                        #contc = np.median(continuum[mask])
                        index = np.where(finalflux == np.max(finalflux))
                        peakWl = finalwl[index]
                        peakFl = finalflux[index]
                        #halfPeak = peakFl/2
                        

                        centB = (self.xpos - 0.05*self.xpos,self.xpos + 0.05*self.xpos)
                        self.plot_view.plot([centB[0],centB[0]],[0,self.ypos],pen = 'c')
                        self.plot_view.plot([centB[1],centB[1]],[0,self.ypos],pen='m')
                        zcent = np.round((peakWl - 1215.67)/1215.67,6) #specific to lyman alpha
                        zB = ((zcent - 0.01*zcent), (zcent + 0.01*zcent))
                        #zB = (2.664020,2.664024)
                        zB = (zB[0][0],zB[1][0])#necessary b/c zB is created as array of arrays and numpy fails with array inputs
                        sigB = (0.005*(self.ewBounds[1] - self.ewBounds[0]), 0.1*(self.ewBounds[1]-self.ewBounds[0]))
                        ampB = (peakFl - 0.01*peakFl,peakFl + 0.01*peakFl)
                        ampB = (ampB[0][0],ampB[1][0])#same as zB

                        self.Fitter(fxn.Lyalph,data,finalflux,finalerr,finalwl,[ampB,zB,sigB],name='EW')
                    else:
                        qt.QMessageBox.about(self,"Warning!","No continua available! First fit continuum")
                else:
                    qt.QMessageBox.about(self,"No data on screen","Not fitting")
                        
        
    #TODO: Need to track and remove plots, perhaps on user request
    def Cont_Fit(self):
        """
        This will call specutils and ask for a wavelength range from the user. A continuum will be fit
        to the flux data with consideration of the error spectrum. Should also consider setting up
        a list of masks to the data to use for continuum fitting
        """
        self.isMask = True
        #self.Mcount, ok = qt.QInputDialog.getInt(self,"","How many masks?")
        if self.contBut.isChecked():
            self.plot_view.setCursor(self.Lcursor)#Setting cursor image
            Pos(self,self.plot_view)
        if not(self.contBut.isChecked()):
            self.isMask = False
            self.plot_view.setCursor(core.Qt.CursorShape(2))
            self.plot_view.scene().sigMouseClicked.disconnect()
            data = self.plot_view.listDataItems()
            #geting data for fitting continuum and plotting chosen points
            if(len(data) > 0):
                data = [data[0].getData(),data[1].getData()]
                trange = self.plot_view.viewRange()
                yrange = trange[1]
                holdToRemove = []
                for i in range(len(self.Mask)):
                    line_data = [[self.Mask[i][0],self.Mask[i][0]],[yrange[0],yrange[1]]]
                    rine_data = [[self.Mask[i][1],self.Mask[i][1]],[yrange[0],yrange[1]]]
                    b = pg.PlotDataItem(x=line_data[0],y=line_data[1],pen='b')
                    r = pg.PlotDataItem(x=rine_data[0],y=rine_data[1],pen='r')#symbol is another argument, best for scatterplotitem
                    self.plot_view.addItem(b)
                    self.plot_view.addItem(r)
                    holdToRemove.append(b)
                    holdToRemove.append(r)
                    fmask = (data[0][0] > self.Mask[i][0]) & (data[0][0] < self.Mask[i][1])
                    #TODO: this does not account for multiple masks         
            flux = data[0][1][fmask]
            wl = data[0][0][fmask]
            err = data[1][1][fmask]

            items = ('Power Law','Linear','Polynomial')
            func, ok = qt.QInputDialog.getItem(self,"Get function","Function: ",items,0,False)
            if func == 'Power Law' and ok:
                self.Fitter(fxn.Pow,data,flux,err,wl,[(np.min(flux),np.max(flux)),(-3,1),(np.min(wl),np.max(wl))],name='continuum')
                self.cname = 'pl'
            if func == 'Linear' and ok:
                self.line = partial(fxn.linear,b=self.Mask[0][0])
                self.Fitter(self.line,data,flux,err,wl,[(np.min(flux),np.max(flux)),(np.min(flux)/np.max(wl),np.max(flux)/np.min(wl))],name='continuum')
                #TODO: is the slope range reasonable?
                self.cname = 'l'
            if func == 'Polynomial' and ok:
                order, check = qt.QInputDialog.getInt(self,"Polynomial","What Order?:",0,False)
                if check:
                    self.cname = 'p'
                    guess = []
                    guess.append((np.min(flux),np.max(flux)))
                    for i in range(order):
                        guess.append(((np.min(flux)/np.max(wl))**(i+1),(np.max(flux)/np.min(wl)**(i+1))))

                    self.Fitter(fxn.polynomial,data,flux,err,wl,guess,name='continuum')
            if not(ok):
                qt.QMessageBox.about(self,"Continuum","Not fitting continuum")
                for i in holdToRemove:
                    self.plot_view.removeItem(i)
                
                
                
            self.Mask = []
        return
        
    
    def Fitter(self,func,data,flux,err,wl,bounds,name = ''):
        '''
        helper function for fitting algorithms (cont_fit, and ew)
        '''
        mymc = mcmc.fit(func,wl,flux,err, 1500,*bounds)
        #TODO:Need to optimize mcmc (use emcee? don't like interface, but perhaps Tim and Rem can clarify)
        # (Using MCMC_fine_tuning.py) maybe rewrite MCMC_fine_tuning in c++ and then import the compiled file
        hists = []
        params = []
        ps = []
        #perc = []
        test = []
  

   
        for i in range(len(mymc[0])):
            #perc.append(np.percentile(mymc[0][i],[16,50,84]))
            hists.append(np.histogram(mymc[0][i],bins = 250)) 
            params.append(np.max(hists[i][0]))
            test.append(np.where(hists[i][0] == np.max(hists[i][0])))#TODO: Is this the best fit? should devise an algorithm that determines best fit given cross correlations
            ps.append(hists[i][1][test[i]])
            plt.figure(i)
            plt.plot(hists[i][1][1:],hists[i][0])

        plt.show()
        hists = np.array(hists)
        value = [conf.ConfInt(hists[i][0],hists[i][1][1:],0.68) for i in range(len(hists))]
        if name == 'continuum':
            self.contfit.append(ps)
            self.conterr.append(value)
            self.pdfs[name] = mymc[0]
            pen = (100,90,0)
            cont = 1.0
        if name == 'EW':
            self.ewfit.append(ps)
            self.ewferr.append(value)
            self.pdfs[name] = mymc[0]
            pen = (0,100,0)
            print("result of EW fit {}".format(self.ewfit))

            #This is used for GUI image such that we can see the fit
            if self.cname == 'pl':
                cont = fxn.Pow(data[0][0],*self.contfit[0])
            elif self.cname == 'l':
                cont = self.line(data[0][0],*self.contfit[0])
            elif self.cname == 'p':
                cont = fxn.polynomial(data[0][0],*self.contfit[0])
            
        self.plot_view.plot(data[0][0],cont*func(data[0][0],*ps),pen=pen)

    def showPDFS(self):
        '''
        Helper function for displaying corner plot of fit parameters.
        This will display the results of a given fit.
        TODO: need to add on a pdf display of the EW pdf derived from the EW fit
        '''
        func, ok = qt.QInputDialog.getItem(self,"Get Fit","Which Fit?: ",self.pdfs.keys(),0,False)
        if ok:
            for i in range(len(self.pdfs[func])):
                if np.frexp(self.pdfs[func][i][0])[1] < -15:#weird function that gives (x1,x2) s.t. n = x1 * 2**x2, note that 2^-15 ~ 10^-5
                    self.pdfs[func][i] = np.log10(self.pdfs[func][i])
            labels = ['a{}'.format(i) for i in range(len(self.pdfs[func]))]
            corner.corner(np.transpose(self.pdfs[func]),bins=100,quantiles=[0.16,0.5,0.84],show_titles=True,
                        labels=labels,verbose=True,plot_contours=True,title_fmt=".3f") #must be (values, # of parameters), (i.e. (365,4) corresponds to a fit with four parameters)
            plt.show()
                
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
                self.plot_view.setCursor(core.Qt.CursorShape(2))
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
        if self.isMask:
            #TODO: if multiple left and right bounds are chosen and then
            # the fitting is canceled, only the most recent plot is removed
            # need to remove all masks that were chosen. TODO additionally, need
            # a check for whether the left bound is actually to the left
            # of the right bound
            if self.count == 0:
                self.plot_view.setCursor(self.Rcursor)
                self.xpos = vb.mapSceneToView(pos).x()
            if self.count == 1:
                self.plot_view.setCursor(self.Lcursor)
                self.ypos = vb.mapSceneToView(pos).x()
                self.XSignal.emit((self.xpos,self.ypos))
                self.count = 0
                return
            self.count += 1
           

        
    @core.pyqtSlot(list)
    def getPos(self,value):
       
        self.xpos = value[0]
        self.ypos = value[1]
        print(self.xpos,self.ypos)
        self.ewBounds = (value[2],value[3])
        self.fitting_1d()
        

    @core.pyqtSlot(tuple)
    def getX(self,value):
        self.Mask.append(value)
        
        


            



if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())