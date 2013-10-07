function Y2d = palm_conv4to2(Y4d)
% Convert a 4D dataset (x,y,z,t) into a 2D array (t,x*y*z)
% that can be used in a multivariate GLM
% 
% Usage:
% Y2d = conv4to2(Y4d);
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012
% http://brainder.org

% tmp = permute(Y4d,[4 1 2 3]);
% siz = size(tmp);
% Y2d = reshape(tmp,[size(tmp,1) prod(siz(2:end))]);

Y2d = reshape(Y4d,numel(Y4d)/size(Y4d,4),size(Y4d,4))';

