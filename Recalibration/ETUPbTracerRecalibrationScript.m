%% Perform the ET (2)535 tracer calibration using datasets and codes from
% Condon et al. (2015) and McLean et al. (2015), reproducing the original
% results or updating the assumed values.

% 1. load Amelin and Davis (2006) 202-205 spiked 981/982/Puratronic Pb data
load AmelinTarantolaData.mat

% 2. Run inverse problem to solve for ICs of 981, 982, and Puratronic Pb
%
% RECALIBRATION: Change the 208/206 of 981 on line 217, all Pb ICs are 
% calibration to this one ratio (see Condon et al., 2015 for details)
AmelinTarantola_Mean_recalibration
figCondonTable1 = uifigure;
uiCondonTable1 = uitable(figCondonTable1, ...
    "Data", buildTable("Condon Table 1", um));

%3. Perform Monte Carlo uncertainty propagation on the inverse problem, 
% to give the 208/206 of 981 the right expected value and uncertainty
%
% RECALIBRATION: Change the IC of NBS 981 and its uncertainty on lines
% 217 and 218
% AmelinTarantola_MC_recalibration

