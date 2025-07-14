%% Create concsM_981 concsM_982 and concsM_Pur variables for
% ET tracer recalibration. concsM are impurity concentrations in ppb.

% Get GDMS data
% Define the path and filename
data_folder = 'Condon et al 2015 (GCA) Electronic Annex';
file_name = 'ts04_Measurements_GDMS.xlsx';

% Construct the full relative path
relative_path_to_file = fullfile('..', data_folder, file_name);

% Load the data
GDMS_data_981 = readmatrix(relative_path_to_file, ...
    "Sheet", "raw data", "Range","C3:D80");