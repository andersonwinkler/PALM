function P = palm_permtree2(Ptree)

Ptree = permtree(Ptree);
P = cell2mat(pickperm(Ptree,{})');

% ==============================================================
function Ptree = permtree(Ptree)
% Permute tree branches

nU   = size(Ptree,1);
ufin = false(nU,1);

% For each branch of the current node
for u = 1:nU,
    if     ~Ptree{u,2}(3), % if [d 0 0] or [d 1 0]
        % If there are permutations to go at this node
        % or deeper, go there
        [Ptree{u,4},Ptree{u,2}(2)] = permtree(Ptree{u,4});
        break; % break this for-loop
        
    elseif ~Ptree{u,2}(2) &&  Ptree{u,2}(3), % if [d 0 1]
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

% ==============================================================
function P = pickperm(Ptree,P)
% Take a tree in a given state and return the
% permutation. This won't permute, only return the
% indices for the already permuted tree.
nU = size(Ptree,1);
if size(Ptree,2) == 4,
    for u = 1:nU,
        P = pickperm(Ptree{u,4},P);
    end
elseif size(Ptree,2) == 2,
    for u = 1:nU,
        P{numel(P)+1} = Ptree{u,1};
    end
end
