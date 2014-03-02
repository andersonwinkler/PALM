function hcp2solar(csvfile,pedfile)
% Takes a "restricted data" CSV file from the HCP and generates
% a pedigree file that can be used in SOLAR.
% 
% Usage:
% hcp2blocks(csvfile,pedfile)
% 
% csvfile : CSV file downloaded from https://db.humanconnectome.org/
%           containing the "Restricted Data"
% pedfile : File name for the pedigree file to be created.
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

% Var for later
N = size(tab,1);
z = zeros(N,1);
o = ones(N,1);

% Founders
mo = ids(:,2);
mo = horzcat(mo,z,z,o,z);   % mother: sex = 1
mo = unique(mo,'rows');
fa = ids(:,3);
fa = horzcat(fa,z,z,2*o,z); % father: sex = 2
fa = unique(fa,'rows');

% Actual subjects
mzidx = strcmpi(tab(:,5),'MZ');
mfm = horzcat(ids(:,2:3),mzidx);
mz = z;
for m = find(mzidx)',
    mz(all(bsxfun(@eq,mfm,mfm(m,:)),2)) = m;
end

% MZ subjects that have its pair missing are
% treated as non-twin
U = unique(mz(mz > 0));
for u = 1:numel(U),
    uidx = mz == U(u);
    if sum(uidx) == 1,
        mz(uidx) = 0;
    end
end

% All subjects assigned to sex = 1 until the
% true sex information becomes available
su = horzcat(ids,o,mz);

% Pedigree, ready to save
ped = vertcat(mo,fa,su);

% Now save!
fid = fopen(pedfile,'w');
fprintf(fid,'id,fa,mo,sex,mztwin\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',ped');
fclose(fid);
