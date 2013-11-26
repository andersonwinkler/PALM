function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: double.m 1143 2008-02-07 19:33:33Z spm 

%
% niftilib $Id: double.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%


out = double(numeric(fa));

