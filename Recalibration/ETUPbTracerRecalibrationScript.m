%% Perform the ET (2)535 tracer calibration using datasets and codes from
% Condon et al. (2015) and McLean et al. (2015), reproducing the original
% results or updating the assumed values.

% 0. Create figures and tabs to hold output from McLean and Condon papers
figCondon = uifigure("Name", "Condon", "Position", [100, 100, 600, 400]);
gridCondon = uigridlayout(figCondon, ...
    "ColumnWidth", {'1x'}, "RowHeight", {'1x'});
tabGroupCondon = uitabgroup(gridCondon);
figMcLean = uifigure("Name", "McLean", "Position", [100, 700, 600, 400]);
gridMcLean = uigridlayout(figMcLean, ...
    "ColumnWidth", {'1x'}, "RowHeight", {'1x'});
tabGroupMcLean = uitabgroup(gridMcLean);

% 1. load Amelin and Davis (2006) 202-205 spiked 981/982/Puratronic Pb data
load AmelinTarantolaData.mat

% 2. Run inverse problem to solve for ICs of 981, 982, and Puratronic Pb
%
% RECALIBRATION: Change the 208/206 of 981 on line 217, all Pb ICs are 
% calibration to this one ratio (see Condon et al., 2015 for details)
AmelinTarantola_Mean_recalibration
umMaxLik = um; % save off maximum likelihood model vector

tabCondonTable1 = uitab(tabGroupCondon, "Title", "Table 1");
table1Condon = buildTable("Condon Table 1", umMaxLik);
uitable(tabCondonTable1, "Data", table1Condon, ...
    "Units", "normalized", "Position", [0, 0, 1 1]);

% %3a. Perform Monte Carlo uncertainty propagation on the inverse problem, 
% % to give the 208/206 of 981 the right expected value and uncertainty
% %
% % RECALIBRATION: Change the IC of NBS 981 and its uncertainty on lines
% % 217 and 218. To do more MC trials, increase nM on line 155.
% rng("shuffle") 
% AmelinTarantola_MC_recalibration

% 3b. Parse MC intercalibration results, calculate measured and systematic
% uncertainty contributions. Note: change input filenames if needed.
parseMCintercal

tabMcLeanTable3 = uitab(tabGroupMcLean, "Title", "Table 3");
%create inputs for Table 3 of McLean
table3McLean = buildTable("McLean Table 3", ...
    umMaxLik(12:20), twosigmaM(12:20), twosigmaT(12:20));
uitable(tabMcLeanTable3, "Data", table3McLean, ...
    "Units", "normalized", "Position", [0, 0, 1 1]);
