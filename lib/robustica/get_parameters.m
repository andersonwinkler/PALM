function argin = get_parameters(n, arguments)

% Gets input parameters for RobustICA algorithm.
%
%
% SYNTAX: argin = get_parameters(nsrc, arguments)
%
%
% OUTPUT: 
%         argin    : a structure with the following fields:
%                 
%          .deftype  deflation type: 'orthogonalization', 'regression'
%                    * default: 'orthogonalization'
%
%          .dimred'   dimensionality reduction in regression if parameter different from zero;
%                    (not used in deflationary orthogonalization)
%                    * default: false
%
%          .kurtsign source kurtosis signs (one element per source);
%                    maximize absolute normalized kurtosis if corresponding element = 0;
%                    * default: zero vector (maximize absolute normalized kurtosis for all sources)
%
%          .maxiter  maximum number of iterations per extracted source;
%                    * default: 1000
%
%          .prewhi   prewhitening (via SVD of the observed data matrix);
%                    * default: true
%
%          .tol      threshold for statistically-significant termination test of the type
%                         ||wn - p*w||/||w|| < tol/sqrt(sample size);   (up to a phase shift p)
%                    termination is also tested by comparing the gradient norm according to: 
%                         ||g|| < tol/sqrt(sample size);
%                    termination test is not used if tol < 0, so that the algorithm runs the maximum
%                    number of iterations (except if optimal step size reaches a null value);
%                    * default: 1e-3
%
%          .verbose  verbose operation if true
%                    * default: false (quiet operation).
%
%          .Wini     extracting vectors initialization for RobustICA iterative search;
%                    if empty, identity matrix of suitable dimensions is used
%                    * default: empty
%
%          .wreal    if different from zero, keep extracting vector real valued by retaining only the
%                    real part of the gradient; useful, for instance, in separating real-valued mixtures
%                    in the frequency domain, as in the RobustICA-f algorithm
%                    * default: false.
%
%
% INPUTS:
%        nsrc      : number of sources 
%
%        arguments : the cell array used as input to the RobustICA algorithm
%                    (see 'robustica.m' for details).
%
%
% HISTORY:
%
%    <Please add modification date here>: - <please add modification details here>
%
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
% 2014/11/21: - created and added to RobustICA package by Vicente Zarzoso
%               (I3S Laboratory, University of Nice Sophia Antipolis, CNRS, France).



%%% default input parameters

argin = struct('deftype', 'orthogonalization', 'dimred',  false,  'kurtsign', zeros(1, n), ...
               'maxiter', 1e3,                 'prewhi',  true,   'tol',      1e-3, ...
               'verbose', false,               'Wini',    eye(n), 'wreal',    false);         

         
%%% parse argument cell array for input parameters

if ~exist('arguments', 'var')
    arguments = {};
end

numargin = length(arguments);
if mod(numargin, 2), error('Input-argument cell array must contain an even number of elements. Please type ''help robustica'' for details.'), end

for i = 1:(numargin/2)

arglabel = arguments{2*i-1};
argvalue = arguments{2*i};

switch arglabel
    case 'deftype'
        switch(argvalue)
            case {'orthogonalization', 'regression'}
                argin.deftype = argvalue;
            otherwise
                error('''deftype'' input parameter must be either ''orthogonalization'' or ''regression''. Please type ''help robustica'' for details.')
        end
    case 'dimred'
        argin.dimred = argvalue;
    case 'kurtsign'
        if length(argvalue) ~= n, error('''kurtsign'' input parameter must contain one element per source. Please type ''help robustica'' for details.'), end
        argin.kurtsign = sign(argvalue);
    case 'maxiter'
        argin.maxiter = argvalue;
    case 'prewhi'
        argin.prewhi = argvalue;
    case 'tol'
        argin.tol = argvalue;
    case 'verbose'
        argin.verbose = argvalue;
    case 'Wini'
        if ~isempty(argvalue)
            [nrowwini, ncolwini] = size(argvalue);
            if (nrowwini == n) & (ncolwini == n), argin.Wini = argvalue;
            else
                error('''Wini'' must be empty or have dimensions nxn, where n is the number of sources. Please type ''help robustica'' for details.')
            end
        end    
    case 'wreal'
        argin.wreal = argvalue;
end % switch
  
end % for i
