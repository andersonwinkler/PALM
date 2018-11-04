function palm_plot(data,design,contrasts,res)
% Take a vector of data, a design, a constrast, then
% plots. If the code detects that what the contrasts
% is testing is an interaction, it constructs
% interaction plots.
% 
% Usage:
% 
% palm_plot(data,design,contrast,res)
% 
% - data     : CSV file containing data.
% - design   : Design matrix file.
% - contrast : Contrast file.
% - res      : Resolution of meshes (for 2-way 
%              interactions between continuous
%              variables).
%
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Nov/2018
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Load inputs
Y = palm_miscread(data);
Y = Y.data;
M = palm_miscread(design);
M = M.data;
C = palm_miscread(contrasts);
C = C.data';
if size(Y,2) > 1
    error('Input data must have just 1 column.');
end
if size(M,1) ~= size(Y,1)
    error('Input data and designs have different number of rows.');
end
if size(C,1) ~= size(M,2)
    error('Size of contrasts don''t match the size of the design.');
end

% res = 20;
% N = 100;
% %A = randn(N,1);
% %B = randn(N,1);
% A = vertcat(...
%     +ones(N/4,1),...
%     -ones(N/4,1),...
%     +ones(N/4,1),...
%     -ones(N/4,1));
% B = vertcat(...
%     +ones(N/2,1),...
%     -ones(N/2,1));
% C = randn(N,1);
% D = randn(N,1);
% I = ones(N,1);
% Y = 10*(A.*B) + C + D + 3*I + randn(100,1)/10;
% M = [A B C D A.*B I];
% C = [0 0 0 0 1 0 ]';

% Current contrast, model partitioning
[X,Z] = palm_partition(M,C,'guttman');
b = [X Z]\Y;

% Find out what are the A and B main effects for 2-way interactions
EVset = nchoosek(1:size(Z,2),2);
A = NaN; B = NaN;
for r = 1:size(EVset,1)
    if corr(X,prod(Z(:,EVset(r,:)),2)) > 1-10*eps
        A = Z(:,EVset(r,1));
        B = Z(:,EVset(r,2));
        Z(:,EVset(r,:)) = [];
        break
    end
end

% Residual forming matrix without interaction and without main effects
Rz = eye(N) - Z*pinv(Z);

if isnan(A(1)) || isnan(B(1))
    % This is not an interaction
    scatter(Rz*X,Rz*Y)
    
else
    % This is an interaction
    uA = unique(A);
    uB = unique(B);
    if     numel(uA) == 2 && numel(uB)  > 2
        % If A has 2 categories and B is continuous
        rB = Rz*B;
        rY = Rz*Y;
        for u = 1:numel(uA)
            idx = A == uA(u);
            scatter(rB(idx),rY(idx),'.');
            hold('on')
        end
        hold('off');
        
    elseif numel(uA)  > 2 && numel(uB) == 2
        % If A is continuous and B has 2 categories
        rA = Rz*A;
        rY = Rz*Y;
        for u = 1:numel(uB)
            idx = B == uB(u);
            scatter(rA(idx),rY(idx),'.');
            hold('on')
        end
        hold('off');
        
    elseif numel(uA) == 2 && numel(uB) == 2
        % If both A and B have 2 categories
        Y = Rz*Y;
        X = zeros(2,2);
        for ua = 1:numel(uA)
            for ub = 1:numel(uB)
                idx = A == uA(ua) & B == uB(ub);
                X(ua,ub) = mean(Y(idx));
            end
        end
        bar(X);
        
    else
        % if A and B are continuous
        [xg, yg] = meshgrid(linspace(min(A),max(A),res),linspace(min(B),max(B),res));
        mesh(xg,yg,xg.*yg*b(1));
        hold('on')
        scatter3(Rz*A,Rz*B,Rz*Y,'.');
        hold('off')
    end
end
