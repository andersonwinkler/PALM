function Br = palm_reindex(varargin)
% Reindexes blocks so that each is numbered natural numbers,
% starting from 1. There are two pure possibilities: the numbering
% goes continuously for each level, crossing block definitions,
% or it restarts at 1 for each new block. A third method is a
% mixture of both, i.e., it restarts at each block, but at the
% last block it is continuous, crossing blocks.
%
% Usage:
% Br = palm_reindex(B,method)
% 
% B     : Block definitions (multilevel).
% meth  : Method for reindexing: 'continuous', 'restart',
%         or 'mixed' as described above.
%         Default: 'mixed'.
% Br    : Reindexed block definitions.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2013
% http://brainder.org

% Default
meth = 'mixed';

% Take care of input and output vars
B  = varargin{1};
Br = zeros(size(B));
if nargin > 1,
    meth = varargin{2};
end

switch meth,
    
    case 'continuous',
        
        % Renumber using sequential natural numbers starting
        % from 1 and crossing blocks. The first level is
        % treated differently as it doesn't have a
        % "previous" block.
        U = unique(B(:,1));
        for u = 1:numel(U),
            idx = B(:,1) == u;
            Br(idx,1) = u;
        end
        for b = 2:size(B,2), % 2nd onwards
            Bb = B(:,b);
            Bp = Br(:,b-1); % previous
            Up = unique(Bp);
            cnt = 1;
            for up = 1:numel(Up),
                idxp = Bp == up;
                U = unique(Bb(idxp));
                for u = 1:numel(U),
                    idx = (Bb == u) & idxp;
                    Br(idx,b) = cnt;
                    cnt = cnt + 1;
                end
            end
        end
        
    case 'restart',
        
        % Renumber using sequential natural numbers
        % starting from 1 but never crossing blocks,
        % restarting instead for each block.
        Br = renumber(B);
        
    case 'mixed'
        
        % This mixes both above
        Ba = palm_reindex(B,'restart');
        Bb = palm_reindex(B,'continuous');
        Br = horzcat(Ba(:,1:end-1),Bb(:,end));
        
    otherwise
        error('Unknown method: %s',meth)
end

% ==============================================
function Br = renumber(B)
% Note that this runs recursively.
B1 = B(:,1);
U  = unique(B1);
nU = numel(U);
Br = zeros(size(B));
for u = 1:nU,
    idx = B1 == U(u);
    Br(idx,1) = u;
    if size(B,2) > 1,
        [Br(idx,2:end)] = renumber(B(idx,2:end));
    end
end
