function Y4d = palm_conv2to4(Y2d,siz)
% Convert a 2D array (t,x*y*z) into a 4D dataset (x,y,z,t)
% 
% Usage:
% Y4d = conv2to4(Y2d,size3d);
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012
% http://brainder.org

% tmp = reshape(Y2d,[size(Y2d,1) siz]);
% Y4d = permute(tmp,[2 3 4 1]);

Y4d = reshape(Y2d',[siz(:)' size(Y2d,1)]);
