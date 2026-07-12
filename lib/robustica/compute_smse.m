function SMSE = compute_smse(S, Se)

% Computes the signal mean square error between original and estimated sources:
%
%      			SMSE = (1/K) sum_{k=1}^K E{ (s_k - se_k)^2 }
%
% after ordering, scaling and phase correction.
%
% Ordering is performed by the "greedy" algorithm described in Section IV.A of
%
% V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust independent component analysis by iterative maximization</a>
% <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the kurtosis contrast with algebraic optimal step size"</a>, 
% IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
%
%
% SYNTAX: SMSE = compute_smse(S, Se);
%
%
% OUPTPUT:
%	 
%    SMSE  : mean square error (source waveform reconstruction error).
%
% INPUTS:
%
%	 S  : actual sources (one source per row)
%
%	 Se : estimated ones (one source per row).
%
%
% HISTORY:
%
%    <Please add modification date here>: - <please add modification details here>
%
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
% 2010/02/10: - created and added to RobustICA package by Vicente Zarzoso
%               (I3S Laboratory, University of Nice Sophia Antipolis, CNRS, France).


[n, T] = size(S);


%%% Perform optimal ordering (via "greedy algorithm"), as well as scaling and phase correction

ampe = sqrt(diag(Se*Se')/T);    % estimated signal amplitudes for scaling

C = 1/T*S*Se';	    % spatial cross-correlation matrix
Cabs = abs(C);

D = zeros(n, n);    % scale matrix
Ph = D;             % phase correction matrix
P = D;              % permutation matrix       


for k = 1:n;
    [maxlin, poslin] = max(Cabs);
    [xxx, poscol] = max(maxlin);
	
    orgsrc = poslin(poscol);    % original source
    estsrc = poscol;            % estimated source
    
	D(orgsrc, orgsrc) = Cabs(orgsrc, estsrc)/ampe(estsrc)^2;    % optimal scaling in the MMSE sense
    Ph(orgsrc, orgsrc) = sign(C(orgsrc, estsrc));               % phase: related to 'sign' of correlation
	P(orgsrc, estsrc) = 1;                                      % permutation
    
    Cabs(:, estsrc) = zeros(n, 1);    %do not refer to that estimated source anymore
    Cabs(orgsrc, :) = zeros(1, n);    %do not refer to that original source either
    
end % for k

% get estimated sources ready for comparison

Se = Ph*D*P*Se;


%%% Obtain signal mean square error

SMSE = mean(mean(abs(S - Se).^2));
