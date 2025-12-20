function ax = palm_viewglm(Y,M,C,style)
% Plot a vector view of a GLM and its fit.
%
% Usage:
% palm_viewglm(Y,M,C,show_submodels)
%
% Inputs:
% Y     : Dependent data (N by 1)
% M     : Design matrix (N by R, for R regressors)
% C     : Contrast vector or matrix (R by S, S<=R).
%         The contrast vector implicitly defines the
%         regressors of interest X and the nuisance
%         regressors Z.
% style : Choose a figure style:
%         1: Shows Y, X, Z, Yhat, and ehat, i.e.,
%            it shows Y, X, Z, Hm*Y, and Rm*Y, where
%            the design matrix M is [X Z].
%         2: Same as (1) but with additional vectors
%            for the two submodels given by:
%            Y = X*b_X + e_X and Y = Z*b_Z + e_Z,
%            i.e., it shows additionally:
%            Hx*Y and Rx*Y, as well as Hz*Y and Rz*Y.
%         3: Shows Z, as well as the models:
%            Y = Z*b_Z + e_Y and X = Z*b_X + e_X
%            i.e., it shows:
%            Z, X, Y, Hz*Y, Rz*Y, Hz*X, and Rz*X.
%
% Outputs:
% ax    : Axes handle. Use it to set custom axes properties.
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
bfac   = 0.05;  % extra room around the axes bounding box
tfac   = 0.075; % space between tip of vector and center of vector label
radius = 0.2;   % radius of the arcs that represent the angles

% Ensure we have the relevant paths
palm_checkprogs;

% Default vector labels and colors
vecs = cell(3,1); % one per style
%    LaTeX label        % Color code
vecs{1} = { ...
    '$\mathbf{Y}$',       '_g';...
    '$\mathbf{X}$',       '_r';...
    '$\mathbf{Z}$',       '_b';...
    '$\mathbf{H_M Y}$',    'm';...
    '$\mathbf{R_M Y}$',    'm'};
vecs{2} = { ...
    '$\mathbf{Y}$',       '_g';...
    '$\mathbf{X}$',       '_r';...
    '$\mathbf{Z}$',       '_b';...
    '$\mathbf{H_M Y}$',    'm';...
    '$\mathbf{R_M Y}$',    'm';...
    '$\mathbf{H_X Y}$',   '^r';...
    '$\mathbf{R_X Y}$',   '^r';...
    '$\mathbf{H_Z Y}$',   '^b';...
    '$\mathbf{R_Z Y}$',   '^b'};
vecs{3} = { ...
    '$\mathbf{Z}$',       '_g';...
    '$\mathbf{X}$',       '_r';...
    '$\mathbf{Y}$',       '_b';...
    '$\mathbf{H_Z X}$',   '^r';...
    '$\mathbf{R_Z X}$',   '^r';...
    '$\mathbf{H_Z Y}$',   '^b';...
    '$\mathbf{R_Z Y}$',   '^b'};
addtxt = false(size(vecs{style},1));

% Default arc labels and colors
arcs = cell(3,1); % one per style
%    LaTeX label                   Color code   Position   Value
arcs{1} = { ...
    '$\cos(\theta_{\mathbf{YX}})$',      '_y',    [],    [];
    '$\cos(\theta_{\mathbf{XZ}})$',      '_m',    [],    []};
arcs{2} = arcs{1};
arcs{3} = { ...
    '$\cos(\theta_{\mathbf{YX}})$',      '_m',    [],    [];
    '$\cos(\theta_{\mathbf{YR_{Z}X}})$',  'm',    [],    []};

% Options for arrow3
arrowopts = {...
    .75,... % arrow width (default = 1)
    1.5 };  % arrow height (in relation to width, default = 3*width)

% Keep track of the original sum of squares
ssqY  = Y'*Y;

% Recast the model into only 2 columns
[X,Z,C] = palm_partition(M,C,'visualization');

% If Z is not provided or has rank = 0, there is no point
% in using style "2". Also, it would cause issues with style "3".
if rank(Z) == 0
    if style == 2
        style = 1;
    elseif style == 3
        emptysymbol;
        warning('Cannot use style "3" if Z is empty or has rank = 0.')
        return
    end
end

% Find A that reduces the df of the model
if any([1,2] == style)

    % Rescale to unit norm (easier on the figures)
    Y = normed(Y);
    X = normed(X);
    Z = normed(Z);

    % Reduce the df proper
    A = palm_reducedf(Y,[X Z],1,true);

elseif style == 3

    % Now that we've recasted the model into X and Z,
    % let's recast again with Z as dependent; if we don't
    % do this, the figure won't represent the fact that
    % in this style, since we residualize Z, it is Z that
    % is the dependent variable in the figure
    [X,Y,~] = palm_partition([X Y],C,'visualization');

    % Rescale to unit norm
    Z = normed(Z);
    X = normed(X);
    Y = normed(Y);

    % Reduce the df proper
    A = palm_reducedf(Z,[X Y],1,false);
end

% Apply A
y = A'*Y;
x = A'*X;
z = A'*Z;
if rank(z) == 0
    m = x;
else
    m = [x z];
end

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
    [~,~,v] = svd(uvw,0);
    uvw     = uvw*v;
    is2d    = true;
else
    is2d    = false;
end
org = zeros(size(uvw)); % origin of the coordinate system

% Vector lengths and distances between them
d0 = sqrt(sum(uvw.^2,2));
dv = sqrt( ...
    (uvw(:,1)-uvw(:,1)').^2 + ...
    (uvw(:,2)-uvw(:,2)').^2 + ...
    (uvw(:,3)-uvw(:,3)').^2);

% Define range of values to plot and aspect ratio
mi   = floor(min(uvw));
ma   = ceil (max(uvw));
box  = [mi-bfac;ma+bfac];
pba  = range([mi;ma]);
pba(pba==0) = 2*bfac;

% Prepare figure axes
ax = newplot;
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

% Plot the vectors and their projections
if any([1 2] == style)

    % Vectors
    if d0(1) > eps(10) % Y
        arrow3(org(1,:),uvw(1,:),vecs{style}{1,2},arrowopts{:});
        addtxt(1) = true;
    end
    if d0(2) > eps(10) % X
        arrow3(org(2,:),uvw(2,:),vecs{style}{2,2},arrowopts{:});
        addtxt(2) = true;
    end
    if d0(3) > eps(10) % Z
        arrow3(org(3,:),uvw(3,:),vecs{style}{3,2},arrowopts{:});
        addtxt(3) = true;
    end
    if d0(4) > eps(10) % Hm*Y
        arrow3(org(4,:),uvw(4,:),vecs{style}{4,2},arrowopts{:});
        addtxt(4) = true;
        if dv(1,4) > eps(10) % projection Y->Hm*Y
            arrow3(uvw(1,:),uvw(4,:),['--',vecs{style}{4,2}],0);
        end
    end
    if d0(5) > eps(10) % Rm*Y
        arrow3(org(5,:),uvw(5,:),vecs{style}{5,2},arrowopts{:});
        addtxt(5) = true;
        if dv(1,5) > eps(10) % projection Y->Rm*Y
            arrow3(uvw(1,:),uvw(5,:),['--',vecs{style}{5,2}],0);
        end
    end

    % Angle arcs
    [arcs{style}{1,3},arcs{style}{1,4}] = ...
        anglesymbol(uvw(1,:),uvw(2,:),1.25*radius,arcs{style}{1,2});
    [arcs{style}{2,3},arcs{style}{2,4}] = ...
        anglesymbol(uvw(2,:),uvw(3,:),0.80*radius,arcs{style}{2,2});
end
if style == 2

    % Vectors only (extra in relation to style 1
    if d0(6) > eps(10) % Hx*Y
        arrow3(org(6,:),uvw(6,:),vecs{style}{6,2},arrowopts{:});
        addtxt(6) = true;
        if dv(1,6) > eps(10) % projection Y->Hx*Y
            arrow3(uvw(1,:),uvw(6,:),[':',vecs{style}{6,2}],0);
        end
    end
    if d0(7) > eps(10) % Rx*Y
        arrow3(org(7,:),uvw(7,:),vecs{style}{7,2},arrowopts{:});
        addtxt(7) = true;
        if dv(1,7) > eps(10) % projection Y->Rx*Y
            arrow3(uvw(1,:),uvw(7,:),[':',vecs{style}{7,2}],0);
        end
    end
    if d0(8) > eps(10) % Hz*Y
        arrow3(org(8,:),uvw(8,:),vecs{style}{8,2},arrowopts{:});
        addtxt(8) = true;
        if dv(1,8) > eps(10) % projection Y->Hz*Y
            arrow3(uvw(1,:),uvw(8,:),[':',vecs{style}{8,2}],0);
        end
    end
    if d0(9) > eps(10) % Rz*Y
        arrow3(org(9,:),uvw(9,:),vecs{style}{9,2},arrowopts{:});
        addtxt(9) = true;
        if dv(1,9) > eps(10) % projection Y->Rz*Y
            arrow3(uvw(1,:),uvw(9,:),[':',vecs{style}{9,2}],0);
        end
    end
end
if style == 3

    % Vectors
    if d0(1) > eps(10) % Z
        arrow3(org(1,:),uvw(1,:),vecs{style}{1,2},arrowopts{:});
        addtxt(1) = true;
    end
    if d0(2) > eps(10) % X
        arrow3(org(2,:),uvw(2,:),vecs{style}{2,2},arrowopts{:});
        addtxt(2) = true;
    end
    if d0(3) > eps(10) % Y
        arrow3(org(3,:),uvw(3,:),vecs{style}{3,2},arrowopts{:});
        addtxt(3) = true;
    end
    if d0(4) > eps(10) % Hz*X
        arrow3(org(4,:),uvw(4,:),vecs{style}{4,2},arrowopts{:});
        addtxt(4) = true;
        if dv(2,4) > eps(10) % projection X->Hz*X
            arrow3(uvw(2,:),uvw(4,:),['--',vecs{style}{4,2}],0);
        end
    end
    if d0(5) > eps(10) % Rz*X
        arrow3(org(5,:),uvw(5,:),vecs{style}{5,2},arrowopts{:});
        addtxt(5) = true;
        if dv(2,5) > eps(10) % projection X->Rz*X
            arrow3(uvw(2,:),uvw(5,:),['--',vecs{style}{5,2}],0);
        end
    end
    if d0(6) > eps(10) % Hz*Y
        arrow3(org(6,:),uvw(6,:),vecs{style}{6,2},arrowopts{:});
        addtxt(6) = true;
        if dv(3,6) > eps(10) % projection Y->Hz*Y
            arrow3(uvw(3,:),uvw(6,:),['--',vecs{style}{6,2}],0);
        end
    end
    if d0(7) > eps(10) % Rz*Y
        arrow3(org(7,:),uvw(7,:),vecs{style}{7,2},arrowopts{:});
        addtxt(7) = true;
        if dv(3,7) > eps(10) % projection Y->Rz*Y
            arrow3(uvw(3,:),uvw(7,:),['--',vecs{style}{7,2}],0);
        end
    end

    % Angle arcs
    [arcs{style}{1,3},arcs{style}{1,4}] = ...
        anglesymbol(uvw(2,:),uvw(3,:),1.25*radius,arcs{style}{1,2});
    [arcs{style}{2,3},arcs{style}{2,4}] = ...
        anglesymbol(uvw(5,:),uvw(7,:),0.80*radius,arcs{style}{2,2});
end

% Add vector names to the plot
textopts = { ...
    'VerticalAlignment',   'middle',...
    'HorizontalAlignment', 'center',...
    'Interpreter',         'latex' };
gap  = normed(uvw')'*tfac;
tpos = uvw + gap;
for v = 1:size(uvw,1)
    if addtxt(v)
        text(tpos(v,1),tpos(v,2),tpos(v,3),vecs{style}{v,1},...
            'Color',colortable(vecs{style}{v,2}),...
            textopts{:});
    end
end

% Add cosines to the plot
for a = 1:size(arcs{style},1)
    gap  = normed(arcs{style}{a,3}')'*tfac;
    tpos = arcs{style}{a,3} + gap;
    text(tpos(1),tpos(2),tpos(3),...
        sprintf('%s: %g',arcs{style}{a,1},arcs{style}{a,4}),...
        'Color',colortable(arcs{style}{a,2}),...
        textopts{:});
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
function M = normed(M)
% Normalize columns to unit norm
normM = sqrt(sum(M.^2,1));
M     = M./(normM + (normM == 0));

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [b,cos_angle] = anglesymbol(u,v,r,rgb)
% Draw symbol for 90 degrees angle
% u, v  : two 3D vectors (row or column)
% r     : radius of the arc
% rgb   : color of the arc
% b     : bisecting vector, scaled by the radius (to add text later)
% angle : angle between u and v

% Line color
if ischar(rgb)
    rgb = colortable(rgb);
end

% Radius
normu = norm(u);
normv = norm(v);
r     = min([r,normu,normv]);

% Angle between u and v
u = u/normu;
v = v/normv;
cos_angle = dot(u,v);
sin_angle = norm(cross(u,v));
angle = atan2(sin_angle,cos_angle); % signed angle, shortest path

% Normalize to [0, pi]
abs_angle = abs(angle);
if abs_angle > pi
    abs_angle = 2*pi - abs_angle;
end

if abs(abs_angle - pi) <= eps(10)

    % If u and v are at 180 deg, we want two squares for the two right angles
    b   = [u(2) v(1) 0]; % bisecting vector
    fac = 4; % scale factor for the square and text
    anglesymbol(u,b,r,rgb);
    anglesymbol(b,v,r,rgb);
    cos_angle = -1;

elseif  abs(mod(abs_angle,pi/2)     ) <= eps(10) || ...
        abs(mod(abs_angle,pi/2)-pi/2) <= eps(10)

    % If u and v are at 90 deg, we show a square filled in yellow
    b = u + v;
    s = r*sqrt(pi)/2; % side of square so that 4 squares together have area of full circle
    fac = 4; % scale factor for the square and text
    s = s/fac; % let's make it further smaller
    corners = [ ...
        0 0 0;
        s*u;
        s*u + s*v;
        s*v];
    fill3( ...
        corners([1:4 1],1), ...
        corners([1:4 1],2), ...
        corners([1:4 1],3), ...
        [1 1 0],... % fill color, yellow
        'EdgeColor',rgb);
    cos_angle = 0;

else

    % If u and v are in any other angle, we show an arc
    b     = u + v;
    fac   = 1; % keep at one
    a     = 0:(pi/180):angle; % one segment every pi/360 = .5 degree
    arc   = r * (sin(angle - a')*v + sin(a')*u)/sin_angle;
    plot3(arc(:,1),arc(:,2),arc(:,3),'Color',rgb);
end
b = b/norm(b)*r/sqrt(fac);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function emptysymbol
newplot;
text(0.5, 0.5,           ...
    '$\emptyset$',       ...
    'FontSize', 80,      ...
    'Color', [1 .5 .5],  ...
    'VerticalAlignment',   'middle', ...
    'HorizontalAlignment', 'center', ...
    'Interpreter',         'latex');
set(gca, ...
    'Visible', 'off', ...
    'XTick',      [], ...
    'YTick',      [], ...
    'XTickLabel', [], ...
    'YTickLabel', []);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
if all(ismember(c,[0 1])) && ~all(c==0) && ~all(c==1)
    % Linear scaling for pure colors
    if beta < 0 % darkening
        c = (1 + beta) * c;
    else        % lightening
        c(c==0) = beta;
    end
else
    % Power-law for other cases
    expn = (1 - min(1-eps, abs(beta)))^sign(beta);
    c = c.^expn;
end