B = [
    1 1 1 1;
    1 1 1 2;
    1 1 2 3;
    1 1 2 4;
    1 1 3 5;
    1 1 3 6;
    1 2 4 7;
    1 2 4 8;
    1 2 5 9;
    1 2 5 10;
    1 2 6 11;
    1 2 6 12];
% B = B(:,1:3);

% B = [
%     1 1;
%     1 1;
%     1 3;
%     1 3;
%     1 5;
%     1 5;
%     1 7;
%     1 7];

B = [
    1 1 1;
    1 1 2;
    1 2 3;
    1 2 4;
    1 3 5;
    1 3 6];
% B = B(:,2:3);

Brows = 3;
nlevels = 4;
B=(1:Brows)';
for j = 1:(nlevels-1),
    B = [ ...
        kron(kron((1:Brows)',ones(Brows,1)),ones(size(B,1)/Brows,1)) ...
        kron(ones(Brows,1),B) ];
end
B = [ones(size(B,1),1) B];
B = B(:,1:end-1);

[~,screw] = sort(rand(size(B,1),1));
% B = B(screw,:);
% M = kron(eye(2),ones(size(B,1)/2,1));
M = randn(size(B,1),3);
% M = ones(size(B,1),1);

Br    = palm_reindex(B,'restart');
Ptree = palm_tree(Br,M);
P = palm_permtree2(Ptree);
% B = B(:,1:end-1);
% Br2    = palm_reindex(B,'restart');
% Ptree2 = palm_tree(Br2,M);
