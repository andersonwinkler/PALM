function hcpinspect(csvfile)
% Takes a "restricted data" CSV file from the HCP and generates
% shows the subjects that had more than one spouse (i.e., parents
% of half-sibs).
% 
% Usage:
% hcp2inspect(csvfile)
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

% Load the data and select what is now needed
tmp = strcsvread(csvfile);
tab = tmp(2:end,1:5);
ids = cell2mat(tab(:,1:3));

% Founders
mo = ids(:,2);
mo = unique(mo);
fa = ids(:,3);
fa = unique(fa);

% Fathers with 2 wives
for f = 1:size(fa,1),
    fidx = (ids(:,3) == fa(f,1));
    wifes = unique(ids(fidx,2));
    if numel(wifes) > 1,
        fprintf('%d: ',fa(f,1));
        fprintf('%d ',wifes);
        fprintf('\n')
    end
end

% Mothers with 2 husbands
for m = 1:size(mo,1),
    midx = (ids(:,2) == mo(m,1));
    husbands = unique(ids(midx,3));
    if numel(husbands) > 1,
        fprintf('%d: ',mo(m,1));
        fprintf('%d ',husbands);
        fprintf('\n')
    end
end