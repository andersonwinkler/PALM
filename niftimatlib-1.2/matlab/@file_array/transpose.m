function varargout = transpose(varargin)
% Transposing is not allowed.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: transpose.m 1143 2008-02-07 19:33:33Z spm 

%
% niftilib $Id: transpose.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%



error('file_array objects can not be transposed.');
