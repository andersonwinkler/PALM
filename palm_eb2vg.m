function [VG,EB,nVG,nEB] = palm_eb2vg(EB,PB)
% Define the variance groups based on the exchangeability blocks.
% 
% Usage:
% [VG,EB,nVG,nEB] = eb2vg(EB,PB)
% 
% Inputs:
% EB  : Vector with the exchangeability blocks.
% PB  : Logical true/false whether whole block permutation
%       should be used (PB=true) or not (PB=false).
% 
% Outputs:
% VG  : Variance groups.
% EB  : Re-indexed exchangeability blocks.
% nVG : Number of variance groups
% nEB : Number of exchangeability blocks.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jul/2013
% http://brainder.org

% See how many blocks are defined and reindex them to 1..nEB
[u,~,EB] = unique(EB);
nEB = numel(u);

if PB,
    
    % If whole block permutation, the variance groups are orderly defined
    % as using one observation per block, in the order as they appear for
    % each block
    VG = zeros(size(EB));
    for b = 1:nEB;
        if b == 1,
            nVG = sum(EB == b);
        elseif nVG ~= sum(EB == b),
            error([
                'For whole-block permutation, all blocks\n' ...
                'need to be of the same size.%s'],'')
        end
        VG(EB == b) = 1:nVG;
    end
    
else
    
    % If within-block permutation, the variance groups are the same as the
    % exchangeability blocks.
    VG = EB;
    nVG = nEB;
end
