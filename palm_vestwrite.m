function palm_vestwrite(filename,M)
% Write an FSL "vest" file, e.g. design matrix,
% t-contrasts or f-contrasts.
%
% vestwrite(filename,M);
%
% filename : File name.
% M        : Matrix.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Number of waves, points and peak-to-peak heights
[nP,nW] = size(M);
PPH = max(M,[],1) - min(M,[],1);

% Formatting string
fstr = horzcat('%0.6e',repmat('\t%0.6e',1,nW-1),'\n');

% Write to the disk
fid = fopen(filename,'w');
fprintf(fid,'/NumWaves\t%d\n',nW);
fprintf(fid,'/NumPoints\t%d\n',nP);
fprintf(fid,horzcat('/PPHeights\t',fstr),PPH);
fprintf(fid,'\n/Matrix\n');
fprintf(fid,fstr,M');
fclose(fid);
