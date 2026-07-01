function [T,cols] = palm_parquetread(filename)
% Read a Parquet file.
%
% [T,cols] = palm_parquetread('data.parquet')
%
% T    : Array containing the values (without column headers)
% cols : Cell array containing column names
%
% This function uses DuckDB as backend.
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Feb/2026
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Check that file exists
if ~exist(filename,'file')
    error('File not found: %s',filename);
end

% Create temporary CSV file (data only, no headers)
csvfile = [tempname() '.csv'];
cmd     = sprintf(['duckdb -c "COPY (SELECT * FROM read_parquet(''%s'')) ' ...
                   'TO ''%s'' (FORMAT CSV, HEADER FALSE, DELIMITER '','');"'], ...
                  filename,csvfile);
[status,output] = system(cmd);
if status ~= 0
    error('DuckDB failed while reading Parquet file %s:\n%s',filename,output);
end

% Get column names from Parquet file
cmd = sprintf('duckdb -csv -c "DESCRIBE (SELECT * FROM read_parquet(''%s''));"',filename);
[status, cols] = system(cmd);
if status ~= 0
    cols = {};
else
    lines = regexp(strtrim(cols),'\r?\n','split'); % keep first field
    lines(cellfun(@isempty,lines)) = [];
    cols = cellfun(@(s)strtrim(strtok(s,',')),lines,'UniformOutput',false);
    cols = cols(2:end); % drop the header
end

% Load the data (assumed numeric only)
T = load(csvfile);

% Clean up
if exist(csvfile,'file')
    delete(csvfile);
end