function o = horzcat(varargin)
% Horizontal concatenation of file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: horzcat.m 1143 2008-02-07 19:33:33Z spm 

%
% niftilib $Id: horzcat.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%


o    = cat(2,varargin{:});
return;

