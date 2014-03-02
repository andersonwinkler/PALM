% If you'd like to use the 131 subjects, uncomment the line below
PAids = [];
% tmp = csvread('twins_list.csv'); ids = tmp(2:end,1); % <=== UNCOMMENT HERE IF NEEDED, SEE THE EMAIL

% Read the HCP data (restricted data file) and make the
% block definitions
PAB = hcp2blocks('../RESTRICTED_winkler_11_30_2013_16_13_44.csv',false,PAids); % <=== CHANGE HERE

% Number of subjects
PAN = size(PAB,1);

% Some design matrix (change this for actual experimental data)
M = rand(PAN,3);  % <=== CHANGE HERE

% Number of permutations
PAnP0 = 100; % <=== CHANGE HERE

% How the set of permutations should be produced?
% If false, it's a cell array, each element a sparse
% permutation matrix. If true, it's an array N by nP
% of permutation indices.
PAidxout = true; % <=== CHANGE HERE

% Use Conditional Monte Carlo? (recommended: false)
PACMC = false;

% Reindex the block definitions to a more friendly format, without
% the subject IDs.
PABr = palm_reindex(PAB,'fixleaves');

% Make the dependence tree
PAPtree = palm_tree(PABr,M);

% Define the variance groups
PAVG = palm_ptree2vg(PAPtree);

% Generate the set of permutations
PAmaxP = palm_maxshuf(PAPtree,'perms');
PAPset = palm_permtree(PAPtree,PAnP0,PACMC,PAidxout);

% If you only need permutations, comment out from this line ==========

% Test vargroups (this should never, ever give an error)
if PAidxout,
    for p = 1:size(PAPset,2),
        if any(PAVG ~= PAVG(PAPset(:,p))),
            error('Error with the variance groups.');
        end
    end
else
    for p = 1:numel(PAPset), %#ok<UNRCH>
        if any(PAVG ~= PAPset{p}*PAVG),
            error('Error with the variance groups.');
        end
    end
end

% Generate the set of sign-flips
PAmaxS = palm_maxshuf(PAPtree,'flips');
PASset = palm_fliptree(PAPtree,PAnP0,PACMC,PAidxout);

% Generate the set of permutations with sign-flips
PAmaxB = palm_maxshuf(PAPtree,'both');
PABset = palm_shuftree(PAPtree,PAN,PAnP0,PACMC,true,true,PAidxout);
