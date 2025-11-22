function palm_plot(varargin)
% Take a vector of data, regressors from a design, then
% make an interaction plot.
%
% Usage:
%
% palm_plot(Y,X,Z,F)
%
% - Y        : Data.
% - X        : Main effects (up to 3 colums, of which
%              no more than 2 can be continuous.
% - Z        : Nuisance. It should not include the
%              interaction that is to be plotted,
%              otherwise the effect of the interaction
%              is washed out.
% - F        : (Optional) A struct with fields that are Figure
%              properties 'title', 'xlabel', 'ylabel' and 'zlabel',
%              to be applied to the plot.
%              For discrete variables, names of the categories
%              can be passed as a cell array of strings
%              in fields 'xnames' and 'ynames' of F.
%              F may also contain fields named:
%              - 'res': Resolution of meshes (for 2-way interactions between
%              continuous variables). Default = 10.
%              - 'poly22' : Creates a curvy plot if 'true'
%              (it won't match the GLM, though). Default is false.
%
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Nov/2018 (first version)
% Nov/2025 (this version)
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

% Defaults
res    = 10;
poly22 = false;
F      = struct();

% Parse arguments
narginchk(3,4);
Y = varargin{1};
X = varargin{2};
Z = varargin{3};
if nargin >= 4
    F = varargin{4};
    if ~isstruct(F)
        error('Fourth argument F must be a struct or omitted.');
    end
    if ~ isfield(F,'res')
        F.res = res;
    elseif ~ isnumeric(F.res) || ~ isscalar(F.res) || F.res < 2
       error('Resolution of meshes must be higher than 2');
    end
    if ~ isfield(F,'poly22')
        F.poly22 = poly22;
    elseif ischar(F.poly22) || isstring(F.poly22)
        F.poly22 = contains(lower(F.poly22),'on','tr','po');
    end
end

% Sanity checks
if size(Y,2) ~= 1
    error('Y must be a column vector.');
end
N = size(Y,1);
if size(X,1) ~= N || (~isempty(Z) && size(Z,1) ~= N)
    error('Y, X, and Z must have the same number of rows.');
end
if size(X,2) < 1 || size(X,2) > 3
    error('X must have 1 to 3 columns.');
end

% Add intercept to nuisance
Z = [Z ones(N,1)];

N = size(Y,1);
colorlist='brgymck';

% Ensure the intercept is always included
Z = [Z ones(size(Z,1),1)];

% Residual forming matrix without interaction and without main effects
Rz = eye(N) - Z*pinv(Z);
J  = size(X,2);
switch J
    
    case 1
        % This is not an interaction
        meY = mean(Y);
        scatter(Rz*X,Rz*Y+meY);

    case 2
        % This is an interaction of 2 variables
        meY = mean(Y);
        rY  = Rz*Y;
        A   = X(:,1);
        B   = X(:,2);
        uA  = unique(A);
        uB  = unique(B);
        if     numel(uA) == 2 && numel(uB)  > 2
            % If A has 2 categories and B is continuous
            rB = Rz*B;
            xlim = [+inf -inf];
            for u = 1:numel(uA)
                idx = A == uA(u);
                scatter(rB(idx),rY(idx)+meY,colorlist(u),'.');
                xlimc = get(gca,'Xlim');
                xlim(1) = min(xlim(1),xlimc(1));
                xlim(2) = max(xlim(2),xlimc(2));
                hold('on')
            end
            for u = 1:numel(uA)
                idx = A == uA(u);
                b = rB(idx)\rY(idx);
                yfit = xlim*b;
                plot(xlim,yfit+meY,colorlist(u));
                hold('on')
            end
            ylim = get(gca,'YLim');
            axis([xlim ylim]);

        elseif numel(uA)  > 2 && numel(uB) == 2
            % If A is continuous and B has 2 categories
            rA = Rz*A;
            xlim = [+inf -inf];
            for u = 1:numel(uB)
                idx = B == uB(u);
                scatter(rA(idx),rY(idx)+meY,colorlist(u),'.');
                xlimc = get(gca,'Xlim');
                xlim(1) = min(xlim(1),xlimc(1));
                xlim(2) = max(xlim(2),xlimc(2));
                hold('on')                
            end
            for u = 1:numel(uB)
                idx = B == uB(u);
                b = rA(idx)\rY(idx);
                yfit = xlim*b;
                plot(xlim,yfit+meY,colorlist(u));
                hold('on')
            end
            hold('off');
            ylim = get(gca,'YLim');
            axis([xlim ylim]);
            
        elseif numel(uA) == 2 && numel(uB) == 2
            % If both A and B have 2 categories
            X   = zeros(2,2);
            seX = X;
            for ua = 1:numel(uA)
                for ub = 1:numel(uB)
                    idx = A == uA(ua) & B == uB(ub);
                    X(ua,ub) = mean(rY(idx)); % Mean for each category
                    seX(ua,ub) = std(rY(idx))/sqrt(sum(idx)); % Std Error for each category
                end
            end
            bar(X+meY); hold on
            ngroups = size(X,1);
            nbars = size(X,2);

            % Calculating the width for each bar group
            groupwidth = min(0.8, nbars/(nbars + 1.5));
            for b = 1:nbars
                xpos = (1:ngroups) - groupwidth/2 + (2*b-1) * groupwidth / (2*nbars);
                errorbar(xpos,X(:,b)+meY,seX(:,b),'.','Color',[0 0 0]);
            end
            ymin = min(X(:)-seX(:)+meY);
            ymax = max(X(:)+seX(:)+meY);
            ydelta = ymax - ymin;
            set(gca,'YLim',[ymin-ydelta*0.10 ymax+ydelta*0.10]);
            hold off
            
        else
            % if A and B are continuous
            rA = Rz*A;
            rB = Rz*B;
            rI = Rz*(A.*B);
            if F.poly22
                surfit = fit([rA rB],rY+meY,'poly22');
                plot(surfit,[rA,rB],rY+meY);
            else
                b = [rA rB rI]\rY;
                [xa,xb] = meshgrid(linspace(min(rA),max(rA),res),linspace(min(rB),max(rB),res));
                xi = xa.*xb;
                mesh(xa,xb,xa*b(1)+xb*b(2)+xi*b(3)+meY);
                hold('on')
                scatter3(rA,rB,rY+meY,'k.');
            end
            hold('off')
        end
        if isfield(F,'title')  title( F.title ); end %#ok<*SEPEX>
        if isfield(F,'xlabel') xlabel(F.xlabel); end
        if isfield(F,'ylabel') ylabel(F.ylabel); end
        if isfield(F,'xnames')
            xticklabels(F.xnames);
        end
        if isfield(F,'ynames')
            legend(F.ynames{:},'Location',F.legend_location);
        end
        hold('off');

    case 3
        % This is an interaction of 3 variables
        % First, identify which one is the discrete
        U  = cell(J,1);
        nU = zeros(J,1);
        for j = 1:J
            U{j}  = unique(X(:,j));
            nU(j) = numel(U{j});
        end
        idxU = find(nU == 2,1,'first');
        D = X(:,idxU); % discrete regressor
        X(:,idxU) = [];
        U = U{idxU};
        for u = 1:numel(U)
            Yu = Y(D == U(u),:);
            Xu = X(D == U(u),:);
            Zu = Z(D == U(u),:);
            subplot(1,2,u);
            palm_plot(Yu,Xu,Zu,F);
        end
end
