function palm_plot(Y,I,X,Z,resol,F)
% Take a vector of data, a design, a constrast, then
% plots. If the code detects that what the contrasts
% is testing is an interaction, it constructs
% interaction plots.
%
% Usage:
%
% palm_plot(data,design,contrast,res)
%
% - Y        : Data.
% - X        : Main effects (up to 3 colums, of which
%              no more than 2 can be continuous.
% - I        : The interaction term, 1 column. Leave
%              empty or NaN if not an interaction.
% - Z        : Nuisance. It should not include the
%              interaction that is to be plotted,
%              otherwise the effect of the interaction
%              is washed out.
% - resol    : Resolution of meshes (for 2-way
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

% Check sanity of inputs
if size(Y,2) > 1
    error('Input data must have just 1 column.');
end
if size(I,2) ~= 1
    error('The interaction term must have just 1 column.');
end
if size(X,2) < 1 || size(X,2) > 3
    error('Input data must have between 1 and 3 columns (inclusive).');
end
if      size(Y,1) ~= size(X,1) || ...
        size(Y,1) ~= size(Z,1) || ...
        size(Y,1) ~= size(I,1)
    error('Input variables must all have the same number of rows.');
else
    N = size(Y,1);
end
if ~isstruct(F) && ~isempty(F) && ~isnan(F),
    error('F must be a struct.')
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

% Model fitting
b = [I X Z]\Y;

% Residual forming matrix without interaction and without main effects
Rz = eye(N) - Z*pinv(Z);

switch size(X,2)
    
    case 1
        % This is not an interaction
        figure
        scatter(Rz*X,Rz*Y)
        
    case 2
        % This is an interaction of 2 variables
        rY = Rz*Y;
        A = X(:,1);
        B = X(:,2);
        uA = unique(A);
        uB = unique(B);
        if     numel(uA) == 2 && numel(uB)  > 2
            % If A has 2 categories and B is continuous
            rB = Rz*B;
            figure
            for u = 1:numel(uA)
                idx = A == uA(u);
                scatter(rB(idx),rY(idx),'.');
                hold('on')
            end
            hold('off');
            
        elseif numel(uA)  > 2 && numel(uB) == 2
            % If A is continuous and B has 2 categories
            rA = Rz*A;
            figure
            for u = 1:numel(uB)
                idx = B == uB(u);
                scatter(rA(idx),rY(idx),'.');
                hold('on')
            end
            hold('off');
            
        elseif numel(uA) == 2 && numel(uB) == 2
            % If both A and B have 2 categories
            X = zeros(2,2);
            for ua = 1:numel(uA)
                for ub = 1:numel(uB)
                    idx = A == uA(ua) & B == uB(ub);
                    X(ua,ub) = mean(rY(idx));
                end
            end
            figure
            bar(X);
            
        else
            % if A and B are continuous
            rA = Rz*A;
            rB = Rz*B;
            [xg,yg] = meshgrid(linspace(min(A),max(A),resol),linspace(min(B),max(B),resol));
            figure
            mesh(xg,yg,xg.*yg*b(1));
            hold('on')
            scatter3(rA,rB,rY,'.');
            hold('off')
            if isstruct(F)
                title(F.title);
                xlabel(F.xlabel);
                ylabel(F.ylabel);
                zlabel(F.zlabel);
            end
        end
        
    case 3
        % This is an interaction of 3 variables
        U  = cell(3,1);
        nU = zeros(3,1);
        for j = 1:size(X,2)
            U{j}  = unique(X(:,j));
            nU(j) = numel(U{j});
        end
        idxU = find(nU == 2,1,'last');
        C = X(:,idxU);
        X(:,idxU) = [];
        U = U{idxU};
        for u = 1:2
            Yu = Y(C == U(u),:);
            Iu = I(C == U(u),:);
            Xu = X(C == U(u),:);
            Xu(:,any(abs(corr(Iu,Xu)) > 1-10*eps,1)) = [];
            Zu = Z(C == U(u),:);
            Zu(:,any(abs(corr([Iu Xu],Zu)) > 1-10*eps,1)) = [];
            Zu(:,any(triu(abs(corr(Zu)))-eye(size(Zu,2)) > 1-10*eps,2)) = [];
            palm_plot(Yu,Iu,Xu,Zu,resol,F)
        end
end
