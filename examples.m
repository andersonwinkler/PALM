%% ==============================================================
% Define the blocks
% This would be loaded as a CSV or a "VEST" file, but below there a number
% of choices to test. Uncomment the respective lines to see what happens
% (some will cause a predictable error due to mixing blocks of different
% sizes).

% % EXAMPLE 1: This is for the whole-block permutation:
% B = [
%     1 1;
%     1 1;
%     1 2;
%     1 2;
%     1 3;
%     1 3;
%     1 4;
%     1 4];

% % EXAMPLE 2: This is for within-block permutation:
% B = [
%     1 1;
%     1 2;
%     1 3;
%     2 4;
%     2 5;
%     2 6;
%     3 7;
%     3 8;
%     3 9];

% % EXAMPLE 3: Now mixing whole with within-block, with multiple-levels:
% B = [
%     1 1 1 1;
%     1 1 1 2;
%     1 1 2 3;
%     1 1 2 4;
%     1 1 3 5;
%     1 1 3 6;
%     1 2 4 7;
%     1 2 4 8;
%     1 2 5 9;
%     1 2 5 10;
%     1 2 6 11;
%     1 2 6 12];
% % Try also removing the first or last columns (uncomment accordingly):
% B = B(1:end-1);
% B = B(2:end);

% % EXAMPLE 4: This is a generic multi-level. Choose how many levels there
% % should be and how many rows at the last level, and it will create a
% % (potentially) very big block-definition matrix:
% nlevels = 4;
% nrows   = 3;
% B=(1:nrows)';
% for j = 1:(nlevels-1),
%     B = [ ...
%         kron(kron((1:nrows)',ones(nrows,1)),ones(size(B,1)/nrows,1)) ...
%         kron(ones(nrows,1),B) ];
% end
% Try also adding a column of ones (make whole-block shuffles at the
% highest level), and remove the last column (whithin-block shuffles at the
% lowest level).
% B = [ones(size(B,1),1) B];
% B = B(:,1:end-1);

% % EXAMPLE 5: This is similar to 3, except that some nodes have just
% % one branch.
% B = [
%     1 1 1 1;
%     1 1 2 2;
%     1 1 3 3;
%     1 1 4 4;
%     1 1 5 5;
%     1 1 6 6;
%     1 2 1 7;
%     1 2 2 8;
%     1 2 3 9;
%     1 2 4 10;
%     1 2 5 11;
%     1 2 6 12];
% % Try also removing the first or last columns (uncomment accordingly):
% B = B(1:end-1);
% B = B(2:end);

%% ==============================================================
% The blocks don't have to be orderly defined as above. This screws up
% the order of the rows, while preserving the definitions.
% It won't affect the results, and shows that the user can supply the
% blocks in any order.
% % [~,screw] = sort(rand(size(B,1),1));
% % B = B(screw,:);

%% ==============================================================
% Design matrix
% Some options for the design matrix:

% % EXAMPLE 1: Multiple groups, no nuisance
% n = size(B,1);
% f = max(factor(n));
% M = kron(eye(f),ones(n/f,1));

% % EXAMPLE 2: Two continuous regressors
M = randn(size(B,1),2);

% % EXAMPLE 3: A column full of ones (so, nothing to permute)
% M = ones(size(B,1),1);

% % EXAMPLE 4: A monotonically increasing variable (keeps the order of perms)
% M = (1:size(B,1))';

%% ==============================================================
% Now the actual code that matters:
% - palm_reindex will renumber the blocks to natural numbers starting
%   from 1, and restarting at each block border.
% - palm_tree will generate the tree with the structure between the
%   observations. It is the branches of this tree that will be permuted.
% - palm_permtree will take the tree previously generated and return a set
%   of permutations that respect this tree structure.
% See each of these functions for details.
nP    = 100;  % number of permutations
Br    = palm_reindex(B);
Ptree = palm_tree(Br,M);
P     = palm_permtree(Ptree,100,true);
