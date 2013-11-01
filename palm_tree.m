function Ptree = palm_tree(B,M)
% Generates a tree that takes the dependence structure
% between observations into account. From this tree, the
% permutations can be generated later.
%
% Usage:
% Ptree = palm_tree(B,M)
%
% - B       : Multi-level block definitions
% - M       : Design matrix.
% - Ptree   : Permutation tree, from which permutations are generated
%             later.
% 
% Each node is a cell with 4 elements:
% N{1,1}    : a vector with a sequence of indices that indicates the
%             current lexicographic permutation. For within-block
%             exchangeability at the first level, shows NaN.
% N{1,2}(1) : A boolean indicating if the last lexicographic
%             permutation of the branches that begin at this node
%             has been reached, i.e., if the sequence N{1,1} is
%             the last one lexicographically.
% N{1,2}(2) : A boolean indicating if the last lexicographic
%             permutation for all the sub-branches that begin at
%             this node have been reached. For the terminal branches,
%             for which there are no deeper nodes/branches, this
%             is a NaN.
% N{1,3}    : The branches that begin here. This cell will never
%             change.
% N{1,4}    : The branches that begin here, i.e., a copy of the
%             previous. However, these are the ones
%             that are shuffled and used to make the permutations.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2013
% http://brainder.org

% Order of the observations in the original data. Note
% that B should have not been sorted, and should retain
% the order as originally entered (even if the rows are
% totally scrambled.
O = (1:size(M,1))';

% Now make the actual tree.
% For the 1st level (note that the levels begin at 0)
recdepth = 0;
Ptree = cell(1,4);
[Ptree{1},Ptree{3}] = maketree(B,M,O,recdepth+1);
Ptree{2} = [0 0];
Ptree{4} = Ptree{3};

% ==============================================================
function [S,Ptree] = maketree(B,M,O,recdepth)
% Now makes the actual tree, each branch recursively.

% Unique blocks at this level
B1 = B(:,1);
U  = unique(B1);
nU = numel(U);

% Some vars for later
Ucnt    = zeros(nU,1);
if size(B,2) > 1,
    Ptree = cell(nU,4);
else
    Ptree = cell(nU,1);
end

% For each block
for u = 1:nU,
    
    % Enter into each unique block
    idx = B1 == U(u);
    
    % If this isn't the last level, continue constructing
    % the branches recursively.
    if size(B,2) > 1,
        [Ptree{u,1},Ptree{u,3}] = ...
            maketree(B(idx,2:end),M(idx,:),O(idx),recdepth+1);
        Ptree{u,4} = Ptree{u,3};
        Ptree{u,2} = [0 0];
        if size(B,2) == 2,
            % For the terminal branches, it's not meaningful to have
            % a true/false here, and this NaN is checked elsewhere.
            Ptree{u,2}(2) = NaN;
        end
    else
        % At the terminal branches, there is no more tree, so
        % just keep track of the observation indices.
        Ptree{u,1} = O(idx);
    end
    
    % This is only to later validate the subblocks, which need
    % to be all of the same size to allow permutation.
    Ucnt(u) = sum(idx);
end

% Check if all subblocks within each EB are of the same size.
if recdepth > 1 && any(diff(Ucnt)),
    
    % Note that counting the number of times Matlab/Octave
    % reports as the source of the error the line above where
    % palm_renumber is called again helps to identify which
    % is the offending block level.
    error('Not all sub-blocks within an EB are of the same size at level %d.\n',recdepth);
end

% Make the sequence of values that are the reference for the
% lexicographic permutations to be done later. This is doesn't
% apply to the highest level.
if recdepth > 1,
    % Identify repeated sets of rows, which receive then
    % the same index; these repetitions are taken care of
    % later by the Algorithm L.
    B1M = sortrows([B1 M]);
    Ms  = B1M(:,2:end);
    S   = palm_mat2seq(reshape(Ms',[numel(Ms)/nU nU])');
    
    % Put in ascending order, and shuffle the branches
    % accordingly
    [S,idx] = sort(S);
    Ptree = Ptree(idx,:);
    
elseif nU == 1,
    % For whole block starting at the second level, the
    % permutation matrix is simply the number 1.
    S = 1;
else
    % There isn't whole-block permutation at the 1st level,
    % only within-block, so this case is marked as NaN.
    S = NaN;
end
