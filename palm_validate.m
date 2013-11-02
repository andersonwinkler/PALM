function palm_validate(B)
% Check the block definitions for a number of mistakes.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2013
% http://brainder.org

% Do the work here.
checkblk(B,B(1)>=0,0);

% ==============================================================
function checkblk(B,wholeblock,recdepth)
% This will check whether:
% - the leftmost column is valid
% - the blocks are of the same size for whole-block permutation
% - the indices aren't 0 or non-integer
% Note this uses recursion.

% Vars for later
B1   = B(:,1);
U    = unique(B1);
nU   = numel(U);
Ucnt = zeros(nU,1);

% Using positive/negative indices implies that the leftmost column must
% exist and be filled entirely with the same digit.
if recdepth == 0 && numel(U) > 1,
    warning('The highest level (leftmost column) must be entirely filled by the same digit.'); %#ok
end

% For each block
for u = 1:nU,
    
    if U(u) == 0,
        
        % The index 0 cannot be considered valid (neither
        % positive or negative).
        warning('The digit 0 is not a valid block indicator (level %d).\n',recdepth); %#ok
        
    elseif rem(U(u),1) ~= 0,
        
        % Let's avoid fractional indices too.
        warning('Non-integer indices should be avoided (level %d).\n',recdepth); %#ok
        
    else
        
        % Enter into each unique block to see what else
        idx = B1 == U(u);
        if size(B,2) > 1,
            
            % Here the test for whole-block permutation allows 0, just so
            % that these otherwise invalid blocks are not ignored. But the
            % warning message above should raise the attention of the user.
            checkblk(B(idx,2:end),B(find(idx,1),1)>=0,recdepth+1);
        end
        
        % Check the size of the block
        Ucnt(u) = sum(idx);
    end
end

% Check if all subblocks within each EB are of the same size.
if wholeblock && any(diff(Ucnt)),
    
    % Note that counting the number of times Matlab/Octave
    % reports as the source of the error the line above where
    % palm_renumber is called again also helps to identify
    % which is the offending block level.
    warning('Not all sub-blocks within an EB are of the same size at level %d.\n',recdepth); %#ok
end
