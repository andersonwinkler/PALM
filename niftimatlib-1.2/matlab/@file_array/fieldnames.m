function t = fieldnames(obj)
% Fieldnames of a file-array object
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: fieldnames.m 1143 2008-02-07 19:33:33Z spm 

%
% niftilib $Id: fieldnames.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%


t = {...
    'fname'
    'dim'
    'dtype'
    'offset'
    'scl_slope'
    'scl_inter'
};
