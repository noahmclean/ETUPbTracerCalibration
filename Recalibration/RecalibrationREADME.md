# EARTHTIME ET(2)535 Tracer Recalibration

The MATLAB code and datasets in this folder are provided to reproduce the ET(2)535 tracer calibration results from the raw data provided in the sister papers Condon et al. (2015) and McLean et al. (2015), both published in *Geochimica et Cosmochimica Acta*.  

The measured Pb and U isotopic data used for the tracer calibration was provided as electronic supplements to the Condon et al. (2015) article and are provided in a separate folder of this repository.  The MATLAB codes used to transform the measured data into the EARTHTIME tracer isotopic composition are provided as electronic supplements to the McLean et al. (2015) article and are also provided in a separate folder of this repository.

This folder contains the the same data used in the original EARTHTIME tracer calibration exercise, provided here as MATLAB workspace files (.mat). The original codes were written to work with these MATLAB data already in the workspace.  The MATLAB codes in this folder have been renamed and edited only to indicate location of important input parameters that might change when re-calibrating the tracer. This might happen, for instance, when a new isotopic composition for Pb and/or U isotopic standards is published. 

## Pb Reference Material Intercalibration

The original EARTHTIME tracer calibration used published data from Amelin and Davis (2006), who mixed the Pb reference materials NBS 981, NBS 982, and Puratronic Pb with their own 202Pb-205Pb tracer. That data, reproduced in Condon et al. (2015), is provided here in the MATLAB workspace file AmelinTarantolaData.mat. The "Tarantola" reference is to the author of the algorithm used to invert the data. 

### Code descriptions:

`AmelinTarantola_Mean_recalibration.m` inverts all the measured data for the 202Pb/205Pb of the Amelin and Davis tracer, as well as the inter-calibrated isotopic compositions of NBS 981, NBS 982, and Puratronic Pb. The ICs of all parameters are traceable to the assumed 208Pb/206Pb ratio of NBS 981. To change this ratio from its original value of 2.1681 (Catanzaro et al., 1968; NIST certificate, 1991), first load the raw data in `AmelinTarantolaData.mat`, then edit line 217 of the code (the variable `r68g_981`) and run the code.  The results appear in the vector `um`.  Elements 12 to 20 of this vector are the maximum likelihood estimates of [204/206_981 207/206_981, 208/206_981, 204/206_982, 207/206_982, 208/206_982, 204/206_Pur, 207/206_Pur, 208/206_Pur]. 

`AmelinTarantola_MC_recalibration.m` performs the same inversion, but propagates uncertainties via the mixed linear and Monte Carlo uncertainty propagation algorithm described in McLean et al. (2015). Monte Carlo trials are output in the file `'MCintercal.txt`.


## References

Amelin, Y., & Davis, W. J. (2006). Isotopic analysis of lead in sub-nanogram quantities by TIMS using a 202Pb–205Pb spike. Journal of Analytical Atomic Spectrometry, 21(10), 1053-1061.

Condon, D. J., Schoene, B., McLean, N. M., Bowring, S. A., & Parrish, R. R. (2015). Metrology and traceability of U–Pb isotope dilution geochronology (EARTHTIME Tracer Calibration Part I). Geochimica et Cosmochimica Acta, 164, 464-480.

McLean, N. M., Condon, D. J., Schoene, B., & Bowring, S. A. (2015). Evaluating uncertainties in the calibration of isotopic reference materials and multi-element isotopic tracers (EARTHTIME Tracer Calibration Part II). Geochimica et Cosmochimica Acta, 164, 481-501.