function [M,PPH] = palm_vestread(filename)
% Read an FSL "vest" file, e.g. design matrix,
% t-contrasts or f-contrasts.
%
% M = vestread(filename);
%
% filename : File name.
% M        : Matrix.
% PPH      : Peak-to-peak heights.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012
% http://brainder.org

% Read the whole file at once
fid = fopen(filename,'r');
tmp = textscan(fid,'%s');
fclose(fid);

% Get the number of columns
nW = tmp{1}(find(strcmp(tmp{1},'/NumWaves'))+1);
nW = str2double(nW);

% Get the number of rows
nP = tmp{1}(find(strcmp(tmp{1},'/NumPoints') | strcmp(tmp{1},'/NumContrasts'))+1);
nP = str2double(nP);

% Get the peak-to-peak heights
pos = find(strcmp(tmp{1},'/PPheights'));
PPH = tmp{1}(pos+1:pos+nW);
PPH = str2double(PPH)'; % if there is no PPH, this returns NaN

% Reshape to a matrix
M = str2double(tmp{1}(find(strcmp(tmp{1},'/Matrix'))+1:end));
M = reshape(M,[nW,nP])';
