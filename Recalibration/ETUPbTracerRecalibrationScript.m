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

tabMcLeanTable4 = uitab(tabGroupMcLean, "Title", "Table 4");
%create inputs for Table 3 of McLean
table4McLean = buildTable("McLean Table 4", rhotot(12:20, 12:20));
uitable(tabMcLeanTable4, "Data", table4McLean, ...
    "Units", "normalized", "Position", [0, 0, 1 1]);


% 4. Blank and tracer "minor Pb isotope" IC using linear regression with
% overdispersion term to account for Pb blank IC variability

% data from Condon et al supplement "ts05_Tracer-Blank Mixing Line.xls"
load Tracer_Blank_Workspace.mat

% must covert alphaPb measurements from old measurements to new ones
% needs original and new 208/206 ratio of 981.
r86_981_old = 2.1681;
r86_981_new = 2.1658;
ET535LoadingBlanks(:,1) = ...
    0.5*((1+2*ET535LoadingBlanks(:,1))*r86_981_new/r86_981_old - 1);
ET2535LoadingBlanks(:,1) = ...
    0.5*((1+2*ET2535LoadingBlanks(:,1))*r86_981_new/r86_981_old - 1);

% initialize ET535 IC for overdispersion calculation
conc205t = 9.884*10^-12;
r45t = 9.000000000000000e-05;
r65t = 6e-4;%3.887398722782950e-04;
r75t = 6e-4;%2.960685875948440e-04;
r85t = 9e-4;%7.443761622554650e-04;
data = ET535LoadingBlanks';

% iteratively solve for overdispersion and linear regression
for i_loop = 1:5
    BlankIC_ifyouknowtracer_recalibration

    trIC45 = 0.00009;   %starting dist from y-int to tracer
    BlankIC_LinearRegression_recalibration

    % re-run
    r45t = trIC45; % to use if doing blank/tracer IC from 1st principles
    r65t = trIC65;
    r75t = trIC75;
    r85t = trIC85;
end % for i_loop

% display ET535 IC data in output uitables

tabCondonTable3 = uitab(tabGroupCondon, "Title", "Table 3");
table3Condon = buildTable("Condon Table 3 ET535Pb", ...
    trIC_ET535, covtrbl_ET535);
uitable(tabCondonTable3, "Data", table3Condon, ...
    "Units", "normalized", "Position", [0, 0, 1 1]);

% make plot of Blank IC Variability, McLean et al. Figure 6
BlankIC_ifyouknowtracer_recalibration
BlankIC_ifyouknowtracer_recalibration_plots
BlankIC_LinearRegression_recalibration
BlankIC_LinearRegression_recalibration_plots

