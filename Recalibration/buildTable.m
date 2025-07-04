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

end


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