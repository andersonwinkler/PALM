function Mset = palm_msetread(fname)
% Read a .mset file, i.e. a file containing
% multiple matrices (2D arrays)
%
% Mset = palm_msetread(filename)
%
% filename : File name.
% Mset     : Cell array with the matrices.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2015
% http://brainder.org

% Read the input file
fid = fopen(fname,'r');
T = textscan(fid,'%s');
fclose(fid);
T = T{1};

% Read each section that begins with '/Matrix'
% and put as an element in a cell array
midx = find(strcmpi(T,'Matrix'));
Mset = cell(size(midx));
cnt = 1;
for m = midx',
    nR = str2double(T{m+1});
    nC = str2double(T{m+2});
    nV = nR*nC;
    mat = str2double(T(m+3:m+3+nV-1));
    Mset{cnt} = reshape(mat,nC,nR)';
    cnt = cnt + 1;
end
