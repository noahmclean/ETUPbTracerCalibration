function tableOut = buildTable(name, varargin)
%FORMATOUTPUTS Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    name string
end

arguments (Input, Repeating)
    varargin
end

arguments (Output)
    tableOut table
end


switch name

    %%%%%%%%%%%%%%%%%%%%%%%
    case "Condon Table 1"
    
    % input arguments: um
    um = cell2mat(varargin);
    
    rowNames = [
        "Pb source metal"
        "U source metal"
        "Pb metal weight"
        "U metal weight"
        "Pb metal purity"
        "U metal purity"
        "204Pb/206Pb"
        "207Pb/206Pb"
        "208Pb/206Pb"
        "238U/235U"
        "238U/206Pb"];
    JMM = [
        "Puratronic"
        "CRM 115"
        "0.89463"
        "5.15157"
        "0.9999890"
        "0.99977"
        compose("%1.7f", um(18))
        compose("%1.6f", um(19))
        compose("%1.5f", um(20))
        "491.548"
        ""];
    RPGSC = [
        "NBS 982"
        "CRM 112a"
        "0.5175"
        "0.25564"
        "0.9999767"
        "0.99975"
        compose("%1.7f", um(15))
        compose("%1.6f", um(16))
        compose("%1.5f", um(17))
        "137.841"
        ""];
    ET = [
        "NBS 981"
        "CRM 112a"
        "0.31973"
        "5.34448"
        "0.9999986"
        "0.99975"
        compose("%1.7f", um(12))
        compose("%1.6f", um(13))
        compose("%1.4f", um(14))        
        "137.841"
        ""];
    tableOut = table(JMM, RPGSC, ET);
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = ["JMM", "RP/GSC", "ET"];
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 3"

    % input arguments: umMaxLik, twosigmaM, twosigmaT
    umMaxLik = cell2mat(varargin(1));
    twosigmaM = cell2mat(varargin(2));
    twosigmaT = cell2mat(varargin(3));
    
    % Define row names and data for McLean Table 3
    rowNames = [
        "981 204Pb/206Pb"
        "981 207Pb/206Pb"
        "981 208Pb/206Pb"
        "982 204Pb/206Pb"
        "982 207Pb/206Pb"
        "982 208Pb/206Pb"
        "Pur. 204Pb/206Pb"
        "Pur. 207Pb/206Pb"
        "Pur. 208Pb/206Pb"];
    nBlocks = [160; ""; ""; 160; ""; ""; 36; ""; ""];
    umMaxColumn = formatColumn(umMaxLik, [7 6 4 7 6 6 7 6 5]);
    twosigmaMColumn = formatColumn(twosigmaM, [7 6 4 7 6 6 7 6 5]);
    twosigmaMColumn(3) = "-";
    twosigmaTColumn = formatColumn(twosigmaT, [6 5 4 6 5 5 6 5 5]);
    
    tableOut = table(umMaxColumn, twosigmaMColumn, twosigmaTColumn, nBlocks);
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = ...
        ["Wtd. Mean", "±2σᵃ", "±2σᵇ", "n (blocks)"];
    

    %%%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 4"
    
    rhotot = cell2mat(varargin);
    rowNames = [
        "981 204Pb/206Pb"
        "981 207Pb/206Pb"
        "981 208Pb/206Pb"
        "982 204Pb/206Pb"
        "982 207Pb/206Pb"
        "982 208Pb/206Pb"
        "Pur. 204Pb/206Pb"
        "Pur. 207Pb/206Pb"
        "Pur. 208Pb/206Pb"];
    colNames = rowNames;
    tableOut = table('Size', [9, 9], 'VariableTypes', repelem("string", 9));
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = colNames;

    for irow = 1:9
        for jcol = 1:9

            if irow == jcol % ones on the diagonal
                tableOut{irow, jcol} = "1";
            elseif irow > jcol % lower triangle
                tableOut{irow, jcol} = ...
                    compose("%1.3f", rhotot(irow, jcol));
            else, tableOut{irow, jcol} = "";
            end % if rowrho == colrho

        end % for colrho
    end % for rowrho

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "Condon Table 3 ET535Pb"

    trIC = cell2mat(varargin(1));
    covtrbl = cell2mat(varargin(2));
    trIC_1s = sqrt(diag(covtrbl(1:4, 1:4)));

    tracerICstrings = [
        "204Pb/205Pb"
        "206Pb/205Pb"
        "207Pb/205Pb"
        "208Pb/205Pb"
        "233U/235U"
        "238U/235U"
        "conc. 205Pb (mol/g)"
        "conc. 235U (mol/g)"
        ];
    rowNames = [
        "ET535 " + tracerICstrings
        "ET2535 202Pb/205Pb"
        "ET2535 " + tracerICstrings];
    colNames = ["Value", "±1σ (abs)"];
    tableOut = table('Size', [17, 2], 'VariableTypes', repelem("string", 2));
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = colNames;
    
    tableOut{1, 1} = compose("%1.5f", trIC(1));
    tableOut{2, 1} = compose("%1.6f", trIC(2));
    tableOut{3, 1} = compose("%1.6f", trIC(3));
    tableOut{4, 1} = compose("%1.6f", trIC(4));
    tableOut{1, 2} = compose("%1.5f", trIC_1s(1));
    tableOut{2, 2} = compose("%1.6f", trIC_1s(2));
    tableOut{3, 2} = compose("%1.6f", trIC_1s(3));
    tableOut{4, 2} = compose("%1.6f", trIC_1s(4));
    

    %%%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 6"
    
    covtrbl = cell2mat(varargin);
    rowNames = [
        "ET535 204Pb/205Pb"
        "ET535 206Pb/205Pb"
        "ET535 207Pb/205Pb"
        "ET535 208Pb/205Pb"
        "blank 206Pb/204Pb"
        "blank 207Pb/204Pb"
        "blank 208Pb/204Pb"];
    colNames = rowNames;
    tableOut = table('Size', [7, 7], 'VariableTypes', repelem("string", 7));
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = colNames;
    
    for irow = 1:7
        for jcol = 1:7

            if irow == jcol % ones on the diagonal
                tableOut{irow, jcol} = "1";
            elseif irow > jcol % lower triangle
                rhoij = covtrbl(irow,jcol)/...
                        sqrt(covtrbl(irow,irow)*covtrbl(jcol,jcol));
                tableOut{irow, jcol} = ...
                    compose("%1.3f", rhoij);
            else, tableOut{irow, jcol} = "";
            end % if rowrho == colrho

        end % for colrho
    end % for rowrho


    %%%%%%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 5 ET535Pb"
    
    trIC_ET535 = cell2mat(varargin(1));
    covtrbl_ET535 = cell2mat(varargin(2));
    mu = cell2mat(varargin(3));    
    xs = cell2mat(varargin(4)); 

    rowNames = [
        "ET535 204Pb/205Pb"
        "ET535 206Pb/205Pb"
        "ET535 207Pb/205Pb"
        "ET535 208Pb/205Pb"];
    colNames = ["ET535 value", "ET535 ±2σ", "ET2535 value",  "ET2535 ±2σ", " ", ...
        "Pb blank value", "Pb blank ±2σ"];

    tableOut = table('Size', [4, 7], 'VariableTypes', repelem("string", 7));
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = colNames;

    tableOut{:,1} = formatColumn(trIC_ET535, [6, 5, 5, 5]);
    tableOut{:,2} = formatColumn(2*sqrt(diag(covtrbl_ET535(1:4,1:4))), ...
                                 [6, 5, 5, 5]);
    tableOut{:,5} = ["206/204"; "207/204"; "208/204"; " "];
    tableOut{:,6} = [formatColumn(mu, [2, 2, 2]); " "];
    tableOut{:,7} = [formatColumn(2*sqrt(diag(xs)), [2, 2, 2]); " "];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "Condon Table 3 ET2535Pb"

    % add ET2535 Pb IC to existing Condon Table 3 table
    tableOut = cell2table(varargin(1));
    tableOut = tableOut.('Var1'){1,1}; % not easy to get back from cell!
    trIC_ET2535 = cell2mat(varargin(2));
    covtrbl = cell2mat(varargin(3));

    tableOut{10:13,1} = formatColumn(trIC_ET2535, [6, 5, 5, 5]);
    tableOut{10:13,2} = formatColumn(sqrt(diag(covtrbl(1:4,1:4))), ...
                        [6, 5, 5, 5]);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 5 ET2535Pb"
    
    % add ET2535 Pb IC to existing McLean Table 5 table
    tableOut = cell2table(varargin(1));
    tableOut = tableOut.('Var1'){1,1}; % not easy to get back from cell!
    trIC_ET2535 = cell2mat(varargin(2));
    covtrbl = cell2mat(varargin(3));

    tableOut{:,3} = formatColumn(trIC_ET2535, [6, 5, 5, 5]);
    tableOut{:,4} = formatColumn(2*sqrt(diag(covtrbl(1:4,1:4))), ...
                    [6, 5, 5, 5]);


    %%%%%%%%%%%%%%%%%%%%%
    case "McLean Table 8"

    um = cell2mat(varargin(1));
    covmeas = cell2mat(varargin(2));
    covtot = cell2mat(varargin(3));

    rowNames = [
        "233U/235U"
        "238U/235U"];
    colNames = ["MLE", "±2σᵃ", "±2σᵇ", "ρᶜ"];
    
    tableOut = table('Size', [2, 4], 'VariableTypes', repelem("string", 4));
    tableOut.Properties.RowNames = rowNames;
    tableOut.Properties.VariableNames = colNames;

    tableOut{:,1} = formatColumn(um(1:2), [6, 8]);
    tableOut{:,2} = formatColumn(2*sqrt(diag(covmeas(1:2,1:2))), [6, 8]);
    tableOut{:,3} = formatColumn(2*sqrt(diag(covtot(1:2,1:2))), [5, 8]);
    tableOut{1,4} = formatColumn(covtot(1,2)/sqrt(covtot(1,1)*covtot(2,2)), 3);
    tableOut{2,4} = "";


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    case "Condon Table 3 U IC"

    % add U IC to Condon et al. Table 3
    tableOut = cell2table(varargin(1));
    tableOut = tableOut.('Var1'){1,1}; % not easy to get back from cell!
    um = cell2mat(varargin(2));
    covtot = cell2mat(varargin(3));

    tableOut{5:6,1} = formatColumn(um(1:2), [6, 8]);
    tableOut{14:15,1} = formatColumn(um(1:2), [6, 8]);
    tableOut{5,2} = compose("%1.6f", sqrt(covtot(1,1)));
    tableOut{6,2} = compose("%1.2e", sqrt(covtot(2,2)));
    tableOut{14,2} = compose("%1.6f", sqrt(covtot(1,1)));
    tableOut{15,2} = compose("%1.1e", sqrt(covtot(2,2)));

end % switch name


%% Helper functions

function stringVec = formatColumn(columnVec, numDigits)
    
    nRows = length(columnVec);
    stringVec = strings(nRows, 1);
    for iRow = 1:nRows
        formatString = "%9." + num2str(numDigits(iRow)) + "f";
        stringVec(iRow) = compose(formatString, columnVec(iRow));
    end % for iRow

end % formatColumn


end % function buildTable