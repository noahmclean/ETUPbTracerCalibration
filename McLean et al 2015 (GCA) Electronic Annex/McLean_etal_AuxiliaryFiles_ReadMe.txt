Auxiliary Material for

Evaluating Uncertainties in the Calibration of Isotopic Reference Materials and Multi-Element Isotopic Tracers (EARTHTIME Tracer Calibration Part II).

Noah McLean, Daniel Condon, Blair Schoene, sam Bowring

(Massachusetts Institute of Technology)

Geochemistry, Geophysics, Geosystems, Submitted September 2013

Introduction

These auxilliary files are MATLAB scripts that perform the algorithms described in this manuscript.  The MATLAB workspaces (.mat files) consisting of parsed data are available from the authors for the EARTHTIME tracer calibration data.  The raw data is provided in Condon et al. Metrology and Traceability of U-Pb Isotope Dilution Geochronology (EARTHTIME Tracer Calibration Part I) auxilliary materials. 

text01_AmelinTarantola_MLE.m evaluates the Pb standard inter-calibration between NBS 981, NBS 982, and Puratronic Pb, using data from Amelin and Davis (2006).

text02_BlankIC_LinearRegression.m regresses a 4-dimensional line through tracer-blank mixing data using the linear regression algorithm of McLean (2013) and overdispersion calculations of Vermeesch (2011).  This script is run iteratively with the next code.

text03_BlankIC_ifyouknowtracer.m evaluates the multivariate overdispersion in the Pb Blank IC dataset, calculating a Pb Blank IC and uncertainty, along with covariance matrices.

text04_TarantolaUIC_Mean.m inverts for the uranium IC (233/235 and 238/235) of the tracer using critical mixtures with U500 and CRM112a and measurements of the tracer uranium alone.  This code evaluates the maximum likelihood estimate of the U IC.

text05_TarantolaUIC.m performs Monte Carlo simulations of the U IC calculation in the previous code, with new simulations changing values of parameters with systematic uncertainties (see text).

text06_inverseGravTrac.m populates variables and parameters used for inverting the tracer-gravimetric solution data to calculate the 235U/205Pb and 202Pb/205Pb of the tracers.

text07_inverseGravTrac_Mean.m calculates the maximum likelihood estimates of the model parameters for this inverse problem.

text08_inverseGravTracMC.m calculates the Monte Carlo simulations of the parameters with systematic uncertainties in the gravimetric solution - tracer inversion problem.

text09_inverseGravTrac_MCfiles.m performs inversion of the gravimetric-tracer mixture data for these Monte Carlo simulations, saving the results in a tab-delimited text file for later parsing and interpretation.