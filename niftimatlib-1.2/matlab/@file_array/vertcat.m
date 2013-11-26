function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: vertcat.m 1143 2008-02-07 19:33:33Z spm 

%
% niftilib $Id: vertcat.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%



o = cat(1,varargin{:});
return;
