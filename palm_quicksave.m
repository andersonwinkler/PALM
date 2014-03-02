function palm_quicksave(X,opts,plm,y,filename)
% An intermediate function to be used many times in palm.m,
% which chooses an adequate mask, then places back the data
% into the points defined by this mask, and save.
% 
% X        : Data to be saved.
% opts,plm : Structs with options and general data.
% y        : Index for the current data in plm.Yset.
% filename : File to be created.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

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
