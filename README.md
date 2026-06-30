# Astrofit
 The purpose of this code is to serve as a one stop for all things astronomy research.
 Any code that is used in the community is meant to be accessible through this hub.
 Current primary functionality for this code is to view and analyze 1D and 2D fits
 files. 

# 1D viewing and fitting
 The primary analysis codes for 1D come through PyMC, which is a code that uses MCMC
 fitting algorithms to find the best fit to data for a given model and to have uncertainties
 on those fit values. In this coding environment this is primarily done to fit continua, 
 emission lines, and absorption lines in 1D spectra. Though in principle this could be 
 expanded. Various quality assurance (QA) plots can also be generated to assess how 
 successful the fits were. Future improvements to the code will also include an
 interactive plot that will overlay the fit for a given parameter set from the fitting
 algorithm. This can be used either as a learning tool for understanding MCMC and 
 parameter spaces in general, or as a means of optimizing fits. This might be used to
 improve the first guess for fit parameters in the future.

 Additionally, spectral energy distribution (SED) fitting using either Hoki or BPASS will
 be implimented in the future. This will likely fall under the 1D fitting and viewing 
 algorithms as the end result of such analysis is a spectrum.

# 2D viewing and fitting
 Currently this code allows viewing and minimal analysis of 2D fits files. There is some
 ability to extract 1D structions from regions of interest (ROIs), but little analysis
 capabilities beyond this at this time. future iterations of the code will include 
 segmentation mapping capabilities. Ideally this will include the ability to change the 
 sigma values to show which pixels are observed for a given sigma value. This will be
 primarily an educational tool that will hopefully help researchers visualize there data.

# Tables
 Currently this code also supports fits tables, but with zero functionality beyond viewing data
 in the tables.
