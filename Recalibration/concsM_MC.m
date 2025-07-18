function [concsM_981, concsM_982, concsM_Pur] = concsM_MC(nM)
%% Create concsM_981 concsM_982 and concsM_Pur variables for
% ET tracer recalibration. concsM are impurity concentrations in ppb.

% set number of Monte Carlo samples (comment to do this externally)
%nM = 1e6;

% Get GDMS data
% Define the path and filename
data_folder = 'Condon et al 2015 (GCA) Electronic Annex';
file_name = 'ts04_Measurements_GDMS_981_982_Pur.xlsx';

% Construct the full relative path
relative_path_to_file = fullfile('..', data_folder, file_name);

% ranges for data from Condon et al. supplement xlsx file
%  ranges for     [NBS 981,  NBS 982,  Puratronic]
variable_ranges = ["C3:D80", "E3:F80", "I3:J80"];

% for each of the three reference materials:
for iRM = 1:length(variable_ranges)

    % Load the data
    GDMS_data = readmatrix(relative_path_to_file, ...
        "Sheet", "raw data", "Range", variable_ranges(iRM));

    % exclude missing element measurements
    GDMS_data = GDMS_data(~isnan(GDMS_data(:,2)), :);
    nElements = size(GDMS_data, 1);

    % distinguish elements reported as "< x" from those reported as "x"
    isLessThan = ~isnan(GDMS_data(:,1));
    nLessThan = sum(isLessThan);
    nTriangle = nElements - nLessThan;

    % filter data by distribution
    GDMS_lessThan = GDMS_data(isLessThan, 2);
    GDMS_triangle = GDMS_data(~isLessThan, 2);

    % create uniformly distributed synthetic data
    lb = zeros(nLessThan, nM);
    ub = repmat(GDMS_lessThan, 1, nM);
    MC_GDMS_Uniform = unifrnd(lb, ub);

    % create triangularly distributed synthetic data
    % for a triangular distribution, a <= c <= b (see Wikipedia)
    a = repmat(GDMS_triangle*0.5, 1, nM);
    b = repmat(GDMS_triangle*1.5, 1, nM);
    c = repmat(GDMS_triangle, 1, nM);
    U = rand(nTriangle, nM);
    Fc = (c - a) ./ (b - a);
    use_first_branch = U < Fc;
    use_second_branch = ~use_first_branch;
    first_branch = a + sqrt(U .* (b-a) .* (c-a));
    second_branch = b - sqrt((1-U) .* (b-a) .* (b-c));

    MC_GDMS_Triangle = first_branch .* use_first_branch + ...
        second_branch .* use_second_branch;

    MC_impurity = sum(MC_GDMS_Uniform) + sum(MC_GDMS_Triangle);

    switch iRM
        case 1
            concsM_981 = MC_impurity;
        case 2
            concsM_982 = MC_impurity;
        case 3
            concsM_Pur = MC_impurity;
    end

end % for iRM = 1:nRMs


end % function