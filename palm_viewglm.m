function h = palm_viewglm(Y,M,C,style)
% Plot a vector view of a GLM and its fit.
%
% Usage:
% h = palm_viewglm(Y,M,C,show_submodels)
%
% Inputs:
% Y : Dependent data (N by 1)
% M : Design matrix (N by R, for R regressors)
% C : Contrast vector or matrix (R by S, S<=R)
% style : Choose a figure style:
%         1: Shows Y, Yhat, and ehat
%         2: Same as (1) but with additional vectors
%            for the two submodels given by:
%            Y = X*b_X + e_X and Y = Z*b_Z + e_Z
%         3: Shows Z, as well as the models:
%            Y = Z*b_Z and X = Z*b_X
%
% Outputs:
% h : Figure handle
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Dec/2025
% http://brainder.org

% Argument sanity
narginchk(3,4);
if size(Y,1) ~= size(M,1)
    error('Y and M must have the same number of rows.');
end
if size(M,2) ~= size(C,1)
    error('Number of columns in M must match number of rows in C.');
end
if nargin < 4
    style = false; % Do not show submodels by default
end

% Defaults
bfac = 0.05; % extra room around the axes bounding box
tfac = 0.05; % space between tip of vector and center of vector label

% Ensure we have the relevant paths
palm_checkprogs;

% Default labels and colors
lc = cell(3,1); % one per style
lc{1} = { ...
    '$\mathbf{Y}$',     '_g';...
    '$\mathbf{X}$',     '^r';...
    '$\mathbf{Z}$',     '^b';...
    '$\mathbf{H_M Y}$', '_m';...
    '$\mathbf{R_M Y}$', '_m'};
lc{2} = { ...
    '$\mathbf{Y}$',     '_g';...
    '$\mathbf{X}$',     '^r';...
    '$\mathbf{Z}$',     '^b';...
    '$\mathbf{H_M Y}$', '_m';...
    '$\mathbf{R_M Y}$', '_m';...
    '$\mathbf{H_X Y}$', '_r';...
    '$\mathbf{R_X Y}$', '_r';...
    '$\mathbf{H_Z Y}$', '_b';...
    '$\mathbf{R_Z Y}$', '_b'};
lc{3} = { ...
    '$\mathbf{Z}$',     '_g';...
    '$\mathbf{X}$',     '^r';...
    '$\mathbf{Y}$',     '^b';...
    '$\mathbf{H_Z X}$', '_y';...
    '$\mathbf{R_Z X}$', '_y';...
    '$\mathbf{H_Z Y}$', '_c';...
    '$\mathbf{R_Z Y}$', '_c'};
addtxt = false(size(lc{style},1));

% Partition the model and ensure X and Z are vectors
% representing their respective subspaces that
% capture the projections of Y in each of them.
% Doing SVD here and recasting C isn't needed, but helps
% with rank-deficient designs entered by the user.
[u,s,v] = svd(M,'econ');
idx = diag(s) ~= 0;
M1 = u(:,idx)*s(idx,idx)*v(:,idx)';
C1 = M1'*pinv(M)'*C;
M  = M1;
C  = C1;
Cn = null(C');
b  = M\Y;
X  = M*(C*C')*b;
Z  = M*(Cn*Cn')*b;
M  = [X,Z];

% Rescale to unit norm (easier on the figures)
ssqY  = Y'*Y;
normY = sqrt(sum(Y.^2,1));
normM = sqrt(sum(M.^2,1));
Y     = Y./(normY + (normY == 0));
M     = M./(normM + (normM == 0));

% Reduce the df of the design
A    = palm_reducedf(Y,M,1);
y    = A'*Y;
for d = 1:3
    if y(d) < 0
        A(:,d) = -A(:,d); % ensure Y is in the all-positive octant
    end
end
y    = A'*Y;
m    = A'*M;
x    = m(:,1);
z    = m(:,2);

% If we want to see the fittings of the submodels:
% Y = X*b_X + e_X and Y = Z*b_Z + e_Z
switch style
    case 1
        yhat  = m*(m\y);
        ehat  = y - yhat;
        uvw   = [y x z yhat ehat]';
    case 2
        yhat  = m*(m\y);
        ehat  = y - yhat;
        yhatx = x*(x\y);
        ehatx = y - yhatx;
        yhatz = z*(z\y);
        ehatz = y - yhatz;
        uvw   = [y x z yhat ehat yhatx ehatx yhatz ehatz]';
    case 3
        xhat  = z*(z\x);
        ehatx = x - xhat;
        yhat  = z*(z\y);
        ehaty = y - yhat;
        uvw   = [z x y xhat ehatx yhat ehaty]';
end

% Prepare vectors to plot, ensure X and Z are in the main horizontal plane
if rank(uvw) == 2
    [~,~,v]=svd(uvw,0);
    uvw  = uvw*v;
    is2d = true;
else
    is2d = false;
end
org = zeros(size(uvw)); % origin of the coordinate system

% Vector lengths and distances between them
d0 = sqrt(sum(uvw.^2,2));
dv = sqrt( ...
    (uvw(:,1)-uvw(:,1)').^2 + ...
    (uvw(:,2)-uvw(:,2)').^2 + ...
    (uvw(:,3)-uvw(:,3)').^2);

% Define range of values to plot and aspect ratio
mi   = min(uvw);
ma   = max(uvw);
mi   = ~(abs(mi)<eps(10)).*sign(mi);
ma   =  (abs(ma)>eps(10)).*sign(ma);
box  = [mi-bfac;ma+bfac];
pba  = range([mi;ma]);
pba(pba==0) = 2*bfac;

% Prepare figure axes
h = figure;
axis(box(:)');
pbaspect(pba);
if is2d
    view(0,90)
    zticklabels([]);
else
    view(145,20)
end
grid on
hold on

% Place the little square indicating whether the angle
% between x and z is 90 degrees
anglesymbol(uvw(2,:),uvw(3,:),.05);

% Plot the vectors and their projections
if any([1 2] == style)
    if d0(1) > eps(10) % Y
        arrow3(org(1,:),uvw(1,:),lc{style}{1,2});
        addtxt(1) = true;
    end
    if d0(2) > eps(10) % X
        arrow3(org(2,:),uvw(2,:),lc{style}{2,2});
        addtxt(2) = true;
    end
    if d0(3) > eps(10) % Z
        arrow3(org(3,:),uvw(3,:),lc{style}{3,2});
        addtxt(3) = true;
    end
    if d0(4) > eps(10) % Hm*Y
        arrow3(org(4,:),uvw(4,:),lc{style}{4,2});
        addtxt(4) = true;
        if dv(1,4) > eps(10) % projection Y->Hm*Y
            arrow3(uvw(1,:),uvw(4,:),['--',lc{style}{4,2}],0);
        end
    end
    if d0(5) > eps(10) % Rm*Y
        arrow3(org(5,:),uvw(5,:),lc{style}{5,2});
        addtxt(5) = true;
        if dv(1,5) > eps(10) % projection Y->Rm*Y
            arrow3(uvw(1,:),uvw(5,:),['--',lc{style}{5,2}],0);
        end
    end
end
if style == 2
    if d0(6) > eps(10) % Hx*Y
        arrow3(org(6,:),uvw(6,:),lc{style}{6,2});
        addtxt(6) = true;
        if dv(1,6) > eps(10) % projection Y->Hx*Y
            arrow3(uvw(1,:),uvw(6,:),[':',lc{style}{6,2}],0);
        end
    end
    if d0(7) > eps(10) % Rx*Y
        arrow3(org(7,:),uvw(7,:),lc{style}{7,2});
        addtxt(7) = true;
        if dv(1,7) > eps(10) % projection Y->Rx*Y
            arrow3(uvw(1,:),uvw(7,:),[':',lc{style}{7,2}],0);
        end
    end
    if d0(8) > eps(10) % Hz*Y
        arrow3(org(8,:),uvw(8,:),lc{style}{8,2});
        addtxt(8) = true;
        if dv(1,8) > eps(10) % projection Y->Hz*Y
            arrow3(uvw(1,:),uvw(8,:),[':',lc{style}{8,2}],0);
        end
    end
    if d0(9) > eps(10) % Rz*Y
        arrow3(org(9,:),uvw(9,:),lc{style}{9,2});
        addtxt(9) = true;
        if dv(1,9) > eps(10) % projection Y->Rz*Y
            arrow3(uvw(1,:),uvw(9,:),[':',lc{style}{9,2}],0);
        end
    end
end
if style == 3
    if d0(1) > eps(10) % Z
        arrow3(org(1,:),uvw(1,:),lc{style}{1,2});
        addtxt(1) = true;
    end
    if d0(2) > eps(10) % X
        arrow3(org(2,:),uvw(2,:),lc{style}{2,2});
        addtxt(2) = true;
    end
    if d0(3) > eps(10) % Y
        arrow3(org(3,:),uvw(3,:),lc{style}{3,2});
        addtxt(3) = true;
    end
    if d0(4) > eps(10) % Hz*X
        arrow3(org(4,:),uvw(4,:),lc{style}{4,2});
        addtxt(4) = true;
        if dv(2,4) > eps(10) % projection X->Hz*X
            arrow3(uvw(2,:),uvw(4,:),[':',lc{style}{4,2}],0);
        end
    end
    if d0(5) > eps(10) % Rz*X
        arrow3(org(5,:),uvw(5,:),lc{style}{5,2});
        addtxt(5) = true;
        if dv(2,5) > eps(10) % projection X->Rz*X
            arrow3(uvw(2,:),uvw(5,:),[':',lc{style}{5,2}],0);
        end
    end
    if d0(6) > eps(10) % Hz*Y
        arrow3(org(6,:),uvw(6,:),lc{style}{6,2});
        addtxt(6) = true;
        if dv(3,6) > eps(10) % projection Y->Hz*Y
            arrow3(uvw(3,:),uvw(6,:),[':',lc{style}{6,2}],0);
        end
    end
    if d0(7) > eps(10) % Rz*Y
        arrow3(org(7,:),uvw(7,:),lc{style}{7,2});
        addtxt(7) = true;
        if dv(3,7) > eps(10) % projection Y->Rz*Y
            arrow3(uvw(3,:),uvw(7,:),[':',lc{style}{7,2}],0);
        end
    end
end

% Plot the vector names
textopts = { ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'Interpreter','latex'};
gap  = uvw./sqrt(sum(uvw.^2,2))*tfac;
tpos = uvw + gap;
for a = 1:size(uvw,1)
    if addtxt(a)
        text(tpos(a,1),tpos(a,2),tpos(a,3),lc{style}{a,1},...
            'Color',colortable(lc{style}{a,2}),...
            textopts{:});
    end
end
hold off

% Print the sums of squares
if any([1 2] == style)
    ssqy  = (y'*y);
    ssqyh = (yhat'*yhat);
    ssqe  = (ehat'*ehat);
    fmt_head = '%-8s %15s %15s\n';
    if ssqyh > 1e-6 && ssqe > 1e-6
        fmt_row  = '%-8s %15.6f %15.6f\n';
    else
        fmt_row  = '%-8s %15.6e %15.6e\n';
    end
    fprintf(fmt_head,'ssq','scaled','original');
    fprintf('%s', repmat('-',[1 45])); fprintf('\n');
    fprintf(fmt_row,'Y',    ssqy,  ssqy*ssqY);  %#ok<CTPCT>
    fprintf(fmt_row,'Yhat', ssqyh, ssqyh*ssqY); %#ok<CTPCT>
    fprintf(fmt_row,'ehat', ssqe,  ssqe*ssqY);  %#ok<CTPCT>
end

% ==============================================================
function rgb = colortable(str)
% Return RGB triplet for a code following the convention from arrow3.m
beta = 0.4; % default factor for lightening/darkening

% Color table from arrow3.m
vc = 'kymcrgbadefhijlnpqstuvzw';
cn = [
    0 0 0;  % k
    1 1 0;  % y
    1 0 1;  % m
    0 1 1;  % c
    1 0 0;  % r
    0 1 0;  % g
    0 0 1;  % b
    0.42 0.59 0.24;  % a
    0.25 0.25 0.25;  % d
    0.00 0.50 0.00;  % e
    0.70 0.13 0.13;  % f
    1.00 0.41 0.71;  % h
    0.29 0.00 0.51;  % i
    0.00 0.66 0.42;  % j
    0.50 0.50 0.50;  % l
    0.50 0.20 0.00;  % n
    0.75 0.75 0.00;  % p
    1.00 0.50 0.00;  % q
    0.00 0.75 0.75;  % s
    0.80 0.34 0.00;  % t
    0.50 0.00 0.13;  % u
    0.75 0.00 0.75;  % v
    0.38 0.74 0.99;  % z
    1 1 1]; % w

% Parse input string
prefix = '';   % '_' or '^'
code   = str;
if numel(str) >= 2 && any(str(1) == '_^')
    prefix = str(1);
    code   = str(2:end);
end
code = lower(code(1));

% Base color
idx = find(vc == code);
if isempty(idx)
    error('Unknown color code: "%s"', str);
end
baseRGB = cn(idx,:);

% Lighen or darken
if isequal(prefix,'^')
    rgb = brighten(baseRGB, +beta);
elseif isequal(prefix,'_')
    rgb = brighten(baseRGB, -beta);
else
    rgb = baseRGB;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function c = brighten(c,beta)
% Lighten or darker the color
% Linear scaling for "pure" colors
if all(ismember(c,[0 1])) && ~all(c==0) && ~all(c==1)
    if beta < 0 % darkening
        c = (1 + beta) * c;
    else        % lightening
        c(c==0) = beta;
    end
else % Power-law for other cases
    expn = (1 - min(1-eps, abs(beta)))^sign(beta);
    c = c.^expn;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function anglesymbol(u,v,s)
% Draw symbol for 90 degrees angle
% u, v : two 3D vectors (row or column)
% s    : length of side of the little square
u = u/norm(u);
v = v/norm(v);
if abs(u*v') > eps(10)
    return
end
corners = [ ...
    0 0 0;
    s*u;
    s*u + s*v;
    s*v];
fill3( ...
    corners([1:4 1],1), ...
    corners([1:4 1],2), ...
    corners([1:4 1],3), ...
    [1 1 0],... % fill color
    'EdgeColor','k');