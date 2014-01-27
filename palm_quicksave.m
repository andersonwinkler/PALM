function P = palm_quicksave(X,flg,opts,plm,y,c,filename)
% An intermediate function to be used many times in palm.m,
% which chooses an adequate mask, then places back the data
% into the points defined by this mask, and save.
% 
% Usage:
% P = palm_quicksave(X,flg,opts,plm,y,c,filename)
% 
% Inputs:
% X        : Data to be saved.
% flg      : A flag that can be:
%            0: meaning that X is to be just saved.
%            1: meaning that X is a P-value, so it may
%               be converted to 1-P or to -log(P) if
%               this was specified by the user.
%            2: meaning that X is a statistic (G or Z)
%               and needs to be converted to a P-value,
%               and then, perhaps, to 1-p or -log(P).
% opts,plm : Structs with options and general data.
% y        : Index for the current data in plm.Yset.
% filename : File to be created.
% 
% Outputs:
% P          : True P-value (not 1-p or -log(p)).
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Modify X accodring to the flag
if flg == 1 && opts.savecdf,
    
    % Convert P to 1-P
    X = 1 - X;
    
elseif flg == 2,
    
    % Convert a statistic to a P-value
    if opts.savecdf,
        
        % CDF (1-P)
        X = palm_gcdf(X,plm.rC(c),plm.df2{y,c});
        
        % Even saving the CDF, the true p-vals may be needed
        if nargout > 0,
            P = palm_gpval(X,plm.rC(c),plm.df2{y,c});
        end
    else
        % Just P
        X = palm_gpval(X,plm.rC(c),plm.df2{y,c});
        if nargout > 0,
            P = X;
        end
    end
end

% Convert to logarithm
if opts.savelogp && any(flg == [1 2]),
    X = -log10(X);
end

% Choose an appropriate mask struct.
if opts.NPC,
    S = plm.maskinter;
else
    if plm.nmasks == 1,
        S = plm.masks{1};
    else
        S = plm.masks{y};
    end
end

% Inject the data and save.
mask         = S.data;
S.data       = double(S.data);
S.data(mask) = X;
S.filename   = filename;
palm_miscwrite(S);
