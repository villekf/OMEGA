function [output] = loadGATEData(filename,variableName, GATE_char)
%LOADGATEDATA Loads GATE raw data
%   Loads saved raw list-mode GATE data, by first trying the specified file
%   format. If that is not found, the other two are checked as well. If no
%   file is found, an error is thrown.

output = [];
for gg = 1 : 3
    if gg == 1
        try
            output = loadStructFromFile([filename, char(GATE_char{gg}) '.mat'], variableName);
            status = true;
        catch
            warning(['Specified file format (' char(GATE_char{gg}) ') could not be found, trying other file formats']);
            status = false;
        end
    else
        try
            output = loadStructFromFile([filename, char(GATE_char{gg}) '.mat'], variableName);
            status = true;
        catch
            status = false;
        end
    end
    if status
        disp(['Loaded ' char(GATE_char{gg}) ' data']);
        break
    elseif gg == 3 && ~status
        error('No measurement data found')
    end
end
end

