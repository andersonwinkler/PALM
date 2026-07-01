function palm_parquetwrite(filename,T,cols)
% Write a numeric array and its column names to a Parquet file.
%
% palm_parquetwrite('data.parquet', T, cols)
%
% filename : Name of the Parquet file to create
% T        : Array containing the values (no column headers)
% cols     : Column names (cell array)
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

% Check inputs
if ~isnumeric(T) || ~isreal(T)
    error('Data must be a real numeric array.');
end
if ~iscell(cols)
    error('Column names  must be a cell array.');
end
if numel(cols) ~= size(T,2)
    error('Number of column names (%d) does not match the number of data columns (%d).', ...
          numel(cols), size(T,2));
end

% Create a temporary CSV file: a header row with the names, then the data
csvfile = [tempname() '.csv'];
fid     = fopen(csvfile, 'w');
if fid == -1
    error('Could not open temporary file for writing: %s', csvfile);
end
hdr = cellfun(@(c) ['"' strrep(c,'"','""') '"'],cols,'UniformOutput',false);
fprintf(fid,'%s\n',strjoin(hdr,','));
fclose(fid);
dlmwrite(csvfile,T,'-append','delimiter',',','precision',17); %#ok<DLMWT>

% Convert the CSV into Parquet, forcing every column to DOUBLE
typelist = strjoin(repmat({'''DOUBLE'''}, 1, numel(cols)), ',');
cmd      = sprintf(['duckdb -c "COPY (SELECT * FROM read_csv(''%s'', ' ...
                    'header = true, delim = '','', types = [%s])) ' ...
                    'TO ''%s'' (FORMAT PARQUET);"'], ...
                   csvfile,typelist,filename);
[status,output] = system(cmd);

% Clean up the temporary file (whether or not DuckDB succeeded)
if exist(csvfile,'file')
    delete(csvfile);
end
if status ~= 0
    error('DuckDB failed while writing Parquet file %s:\n%s', filename, output);
end