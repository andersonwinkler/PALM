function [Ptree,P] = palm_permtree(Ptree0,nP,flagnorep)

% if this branch is marked as finished, move to next branch and permute
% if the last perm of this branch was done, then mark as finished, move to next branch and permute
% if all branches are finished, 

Ptree = permtree(Ptree);
P = pickperm(Ptree);

function [Ptree,flagfin] = permtree(Ptree)
% This will shuffle the tree, but not yet make the pemutation vector
nU   = size(Ptree,1); % number of branches
ufin = false(nU,1);   % branch finished perms?

% For each branch in the current level
for u = 1:nU,
    
    % If there are permutations to go at this node or deeper, go there
    if     ~Ptree{u,2}(3), % if [d 0 0] or [d 1 0]
        [Ptree{u,4},Ptree{u,2}(2)] = permtree(Ptree{u,4});
        break;
        
    elseif ~Ptree{u,2}(1) &&  Ptree{u,2}(2), % if [0 1]
        
        % If all deeper nodes are known to have finished, but this node
        % itself can be shuffled, permute it, then start permuting these
        % deeper nodes again.
        tmp = palm_nextperm(Ptree{u,1});
        if tmp(1),
            Ptree{u,1} = tmp;
        else
            Ptree{u,2}(1) = true;
        end
        Ptree{u,4} = Ptree{u,3}(Ptree{u,1},:);
        Ptree{u,2}(2) = false;
        break;
        
    %elseif  Ptree{u,2}(1) && ~Ptree{u,2}(2), % if [1 0]
        
        % Do the same as if [0 0]
        
    %elseif  Ptree{u,2}(1) && ~Ptree{u,2}(2), % if [1 1]
        
        % If the permutations of the branches that start at this node 
        % are all over, and also the deeper levels are known to have 
        % finished, do nothing.
        
    end
    
    % If it reached here, this means the loop wasn't broken, which means
    % there were no more permutations to do at the deeper nodes
    ufin(u) = true;
end

% Pass along to the upper level the information that all the branches at
% this node finished (or not).
flagfin = all(ufin);

function Pcell = pickperm(Ptree,Pcell)
% This takes a tree and collects the permutation

nU   = size(Ptree,1);
for u = 1:nU,
    if 
    Pcell = pickperm(Ptree,Pcell);
    else
        Pcell{end+1} = Ptree{}
    end
end