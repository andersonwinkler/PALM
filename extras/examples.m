close all; clear; clc;
% ==============================================================
% Define the blocks
% This would be loaded as a CSV or a "VEST" file, but below there a number
% of choices to test. Uncomment the respective lines to see what happens
% (some will cause a predictable error due to mixing blocks of different
% sizes).

% % EXAMPLE 1: Simple no-block permutation (N!)
% B = [ones(5,1) (1:5)'];

% % EXAMPLE 2: This is for the whole-block permutation:
% B = [
%     1 1 1;
%     1 1 1;
%     1 2 1;
%     1 2 1;
%     1 3 1;
%     1 3 1
%     1 4 1;
%     1 4 1];

% % EXAMPLE 3: This is for within-block permutation:
% B = [
%     -1 1 1;
%     -1 1 2;
%     -1 1 3;
%     -1 2 1;
%     -1 2 2;
%     -1 2 3;
%     -1 3 1;
%     -1 3 2;
%     -1 3 3];

% % EXAMPLE 4: Now mixing whole with within-block, with multiple-levels:
% B = [
%     -1 -1 1 1;
%     -1 -1 1 2;
%     -1 -1 3 1;
%     -1 -1 3 2;
%     -1  2 4 1;
%     -1  2 4 2;
%     -1  2 6 1;
%     -1  2 6 2];
% % Try also removing the last column (uncomment accordingly):
% B = B(:,1:end-1);

% EXAMPLE 5: This is a generic multi-level. Choose how many levels there
% should be and how many rows at the last level, and it will create a
% (potentially) very big block-definition matrix:
% nlevels = 4;
% nrows   = 4;
% B=(1:nrows)';
% for j = 1:(nlevels-1),
%     B = [ ...
%         kron(kron((1:nrows)',ones(nrows,1)),ones(size(B,1)/nrows,1)) ...
%         kron(ones(nrows,1),B) ];
% end
% B = [ones(size(B,1),1) B];
% % Try also removing the last column (uncomment accordingly):
% B = B(:,1:end-1);

% EXAMPLE 6: Tom's small synthetic HCP data
B = [
    -1 1 -1 1 1
    -1 1 -1 1 2
    -1 1 -1 2 1
    -1 1 -1 2 2
    -1 1 -2 1 1
    -1 1 -2 1 2
    -1 1 -2 2 1
    -1 1 -2 2 2
    -1 2 -3 1 1
    -1 2 -3 1 2
    -1 2 -3 2 1
    -1 2 -3 2 2
    -1 3 -4 1 1
    -1 3 -4 1 2
    -1 3 -4 2 1
    -1 4  5 1 1
    -1 4  5 1 2
    -1 4  5 1 3];

% ==============================================================
% The blocks don't have to be orderly defined as above. This screws up
% the order of the rows, while preserving the definitions.
% It won't affect the results, and shows that the user can supply the
% blocks in any order.
% [~,screw] = sort(rand(size(B,1),1)); 
% B = B(screw,:);

% ==============================================================
% Design matrix
% Some options for the design matrix:

% % EXAMPLE 1: Multiple groups, no nuisance
% n = size(B,1);
% f = max(factor(n));
% M = kron(eye(f),ones(n/f,1));

% % EXAMPLE 2: Two continuous regressors
% M = randn(size(B,1),2);

% % EXAMPLE 3: A column full of ones (so, nothing to permute)
% M = ones(size(B,1),1);

% % EXAMPLE 4: A monotonically increasing variable (keeps the order of perms)
M = (1:size(B,1))';

% ==============================================================
% Now the code that matters:
% - palm_reindex will renumber the blocks to natural numbers
%   starting from 1, and restarting at each block border.
%   This version will also fix the terminal branches (leaves)
%   when some simplifications are present in the user-supplied
%   block definitions (here, the variable B).
% - palm_tree will generate the tree with the structure between
%   the observations. It is the branches of this tree that will
%   be permuted or sign-flipped.
% - palm_ptree2vg will generate variance groups.
% - palm_permtree will take the tree previously generated tree
%   and return a set of permutations that respect the structure.
% - palm_fliptree will do the same, but for sign-flips.
% - palm_shuftree will do the same, for permutation and/or
%   sign-flips together. It can also work as a wrapper for
%   palm_permtree and palm_fliptree.

% See each of these functions for details.
nP0    = 40;      % requested number of permutations
N      = size(M,1); % rows in the design matrix
CMC    = true;     % Conditional Monte Carlo
idxout = true;      % output permutation indices
Br     = palm_reindex(B,'fixleaves');
Ptree  = palm_tree(Br,M);
VG     = palm_ptree2vg(Ptree);

% ====[..........Permutations..........]======================
maxP = palm_maxshuf(Ptree,'perms');
Pset = palm_permtree(Ptree,nP0,CMC,idxout);

% Test vargroups (this should never give an error)
if idxout,
    for p = 1:size(Pset,2),
        if any(VG ~= VG(Pset(:,p))),
            error('Error with vargroups');
        end
    end
else
    for p = 1:numel(Pset),
        if any(VG ~= Pset{p}*VG),
            error('Error with vargroups');
        end
    end
end

% ====[...........Sign-flips...........]======================
maxS = palm_maxshuf(Ptree,'flips');
Sset = palm_fliptree(Ptree,nP0,CMC,idxout);

% ====[..Permutations with sign-flips..]======================
maxB = palm_maxshuf(Ptree,'both');
Bset = palm_shuftree(Ptree,N,nP0,CMC,true,true,idxout);


% Show graphically the permutations
if idxout,
    subplot(1,2,1); imagesc(Pset);
    title('Permutations');
    xlabel('Permutation index'); ylabel('Observation index');
    subplot(1,2,2); imagesc(Sset);
    title('Sign-flips');
    xlabel('Sign-flip index'); ylabel('Observation index');
end

% If exhaustive, make sure there are no shufflings missing
if nP0 == 0,
    sizP = size(unique(Pset','rows'),1);
    if sizP ~= maxP,
        error('Number of permutations produced doesn''t match the expected.');
    end
    sizS = size(unique(Sset','rows'),1);
    if sizS ~= maxS,
        error('Number of sign-flips produced doesn''t match the expected.');
    end
end
