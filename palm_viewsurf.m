function fighandle = palm_viewsurf(varargin)
% Simple visualisation of surface results still in Octave/MATLAB.
% The principles are generally similar to those of the tool "dpx2map.m",
% described here:
% https://brainder.org/2013/07/28/displaying-vertexwise-and-facewise-brain-maps/
% 
% Usage:
% h = palm_viewsurf( ...
%                   { lh_dat_file, rh_dat_file }, ...
%                   { lh_srf_file, rh_srf_file }, ...
%                   Name, Value, ...)
% 
% The first argument is a cell-array with two elements: the first is the
% filename of the data to visualise for the left hemisphere, the second is
% the filename of the data for the right hemisphere.
% 
% The second argument is also a cell-array with two elements: the first
% is the surface file for the left hemisphere, the second is the surface
% file for the right hemisphere.
%
% File formats for these inputs are those accepted by PALM.
% 
% All subsequent arguments are supplied as pairs "Name,Value", as typical
% in many MATLAB/Octave commands. Accepted Names and descriptions are:
% 
% MapName     : A MATLAB/Octave colourmap (default: 'viridis').
% DataRange   : Interval to be used to define the colourscale [min max].
%               If not specified, it uses the min and max of the DPX file.
% ShowRange   : Interval to be shown (coloured) [min max]. If not specified,
%               it uses the same as datarange.
% Dual        : True/False. If true, applies the map to the values of no
%               overlap between datarange and showrange. Useful for
%               thresholded positive+negative maps. Default is false.
% ColourGap   : Colour for values that off the colourscale, including NaN.
%               Default is light gray, 25%, i.e. [.75 .75 .75];
% COption     : True/False. The behavior varies if 'dual' is true or not.
%               - For 'dual' = true:
%                 If coption is false, don't rescale the extremities of the
%                 colourmap. Default is true, i.e., produce a higher contrast.
%               - For 'dual' = false:
%                 If coption is false, show the out-of-range values with the
%                 colour specified by 'colourgap'. Default is true, so the
%                 out-of-range values are shown with the extremities of the
%                 colourmap.
% MapSize     : Maximum number of colours in the colourmap. Default 2^16.
% Title       : Title of the figure (to appear at the top of the colourbar).
% Layout      : How to distribute the brain views in the page? The
%               following layouts are available:
%               - 'Simple'      : A simple layout meant to be interactively 
%                                 explored using the MATLAB/Octave figure tools.
%               - 'Cardinal'    : Show 4 views (lateral and medial of left and right
%                                 hemispheres), plus colorbar in the centre.
%               - 'Left'        : Left hemisphere only (lateral and medial).
%               - 'Right'       : Right hemisphere only (lateral and medial).
%               - 'Strip'       : Similar to cardinal but shows all views in a
%                                 single row. Useful for stacking subjects.
%               - 'Publication' : A well rounded layout, meant to be (nearly) ready
%                                 for publication. It always includes a colourbar.
%               - 'Worsley'     : The layout proposed by K. Worsley. Always includes
%                                 a colourbar.
% Background : A 3-element vector of the RGB color of the background.
%              Default is [1 1 1] (white).
% CamLight   : Boolean indicating whether a light should be placed in the
%              same location as the camera.
% Lightning  : Lighting mode. Valid options are:
%              - 'None'
%              - 'Flat'
%              - 'Gouraud'
% Material   : Material mode. Valid options are:
%              - 'Dull'
%              - 'Shiny'
%              - 'Metal'
% Shading    : Shading mode. Valid options are:
%              - 'Flat'
%              - 'Interp'
%              - 'Faceted'
% ColourBar  : Boolean indicating whether the colour bar should be shown.
%              Note that for the 'publication' and 'worsley' layouts, the
%              colour bar is always shown.
% 
% Example:
% freesurfer_dir = getenv('FREESURFER_DIR');
% data = { ...
%     fullfile(freesurfer_dir,'subjects','fsaverage','surf','lh.thickness'), ...
%     fullfile(freesurfer_dir,'subjects','fsaverage','surf','rh.thickness') };
% surfs = { ...
%     fullfile(freesurfer_dir,'subjects','fsaverage','surf','lh.pial'), ...
%     fullfile(freesurfer_dir,'subjects','fsaverage','surf','rh.pial') };
% palm_viewsurf(data,surfs,'layout','cardinal','background',[.5 .5 .5],'colormap','inferno');
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Mar/2021
% http://brainder.org

% Defaults
opts.mapname    = 'viridis';
opts.datarange  = []; % clim, datarange
opts.showrange  = [];
opts.dual       = false;
opts.colourgap  = [.75 .75 .75];
opts.coption    = true;
opts.mapsize    = 2^14;
opts.title      = '';
opts.layout     = 'publication';
opts.background = [1 1 1];
opts.camlight   = true;
opts.lightning  = 'gouraud'; % {'none', 'flat', 'gouraud'}
opts.material   = 'dull';    % {'shiny', 'dull', 'metal'}
opts.shading    = 'interp';  % {'flat', 'interp', 'faceted'}
opts.colourbar  = true;

% Parse inputs and check for most common errors.
dat = varargin{1};
srf = varargin{2};
for a = 3:2:nargin
    opts.(lower(varargin{a})) = varargin{a + 1};
end
if numel(opts.datarange) ~= 2 && ~ isempty(opts.datarange)
    error('The supplied DataRange has to be a two-element vector, or be empty.');
end
if numel(opts.showrange) ~= 2 && ~ isempty(opts.showrange)
    error('The supplied ShowRange has to be a two-element vector, or be empty.');
end
if numel(opts.background) ~= 3 || numel(opts.colourgap) ~= 3
    error('The supplied ColourGap or Background has to be a two-element vector, or be empty.');
end
if ~ islogical(opts.dual) || ~ islogical(opts.coption)
    error('The supplied Dual or COption must be a logical type (i.e., true/false).');
end
if ~ ischar(opts.title)
    error('The supplied Title must be a string.');
end
if ~ islogical(opts.camlight)
    error('The CamLight option must be a logical type (i.e., true/false).');
end
if ~ any(strcmpi(opts.lightning,{'none','flat','gouraud'}))
    error('Unknown Lightning style: %s',opts.lightning);
end
if ~ any(strcmpi(opts.material,{'shiny','dull','metal'}))
    error('Unknown Material style: %s',opts.material);
end
if ~ any(strcmpi(opts.shading,{'flat','interp','faceted'}))
    error('Unknown Shading style: %s',opts.shading);
end
if ~ islogical(opts.colourbar)
    error('The Colourbar option must be a logical type (i.e., true/false).');
end
if strcmpi(opts.layout,'worsley')
    opts.mapname = 'spectral';
end
if palm_isoctave && opts.mapsize > 2^14
    error('Octave cannot plot more than 2^14 = 32768 colors.');
end

% Ensure we gave the programs needed to read the data
palm_checkprogs;

% Read input files
dat{1} = palm_miscread(dat{1});
dat{2} = palm_miscread(dat{2});
srf{1} = palm_miscread(srf{1});
srf{2} = palm_miscread(srf{2});

% Merge left and right hemispheres (to be treated as a "third" hemisphere)
srf{3}.data.vtx = vertcat(srf{1}.data.vtx,srf{2}.data.vtx);
srf{3}.data.fac = vertcat(srf{1}.data.fac,srf{2}.data.fac+size(srf{1}.data.vtx,1));
dat{3}.data = vertcat(dat{1}.data,dat{2}.data);
for h = 1:3
    dat{h}.data = squeeze(dat{h}.data);
end

% Open a new figure window. The handle will be returned.
fighandle = figure;

% Ensure data and colour ranges are well behaved
if isempty(opts.datarange)
    iinf = isinf(dat{3}.data);
    inan = isnan(dat{3}.data);
    iidx = ~(iinf | inan);
    opts.datarange = [min(dat{3}.data(iidx)) max(dat{3}.data(iidx))];
end
if isempty(opts.showrange)
    opts.showrange = opts.datarange;
end
opts.showrange(1) = max(opts.showrange(1),opts.datarange(1));
opts.showrange(2) = min(opts.showrange(2),opts.datarange(2));
for h = 1:numel(dat)
    infneg = isinf(dat{h}.data) & dat{h}.data < 0;
    infpos = isinf(dat{h}.data) & dat{h}.data > 0;
    dat{h}.data(infneg) = opts.datarange(1);
    dat{h}.data(infpos) = opts.datarange(2);
end

% Define the colourmap
if opts.dual
    if opts.coption
        neg = opts.showrange(1) - opts.datarange(1);
        cen = opts.showrange(2) - opts.showrange(1);
        pos = opts.datarange(2) - opts.showrange(2);
        tot = opts.datarange(2) - opts.datarange(1);
        fracs = round([neg cen pos]./tot*opts.mapsize);
        maptmp = eval(sprintf('%s(%d)',opts.mapname,fracs(1)+fracs(3)));
        map = nan(opts.mapsize,3);
        map(1:fracs(1),:) = maptmp(1:fracs(1),:);
        map(fracs(1)+1:fracs(1)+fracs(2),:) = repmat(opts.colourgap,[fracs(2) 1]);
        map(fracs(1)+fracs(2)+1:sum(fracs),:) = maptmp(fracs(1)+1:fracs(1)+fracs(3),:);
    else
        map = eval(sprintf('%s(%d)',opts.mapname,opts.mapsize));
        idxmin = round((opts.showrange(1) - opts.datarange(1)) / (opts.datarange(2) - opts.datarange(1)) * (opts.mapsize - 1) + 1);
        idxmax = round((opts.showrange(2) - opts.datarange(1)) / (opts.datarange(2) - opts.datarange(1)) * (opts.mapsize - 1) + 1);
        map(idxmin:idxmax,:) = repmat(opts.colourgap,[idxmax-idxmin+1 1]);
    end
else
    neg = opts.showrange(1) - opts.datarange(1);
    cen = opts.showrange(2) - opts.showrange(1);
    pos = opts.datarange(2) - opts.showrange(2);
    tot = opts.datarange(2) - opts.datarange(1);
    fracs = round([neg cen pos]./tot*opts.mapsize);
    maptmp = eval(sprintf('%s(%d)',opts.mapname,fracs(2)));
    map = nan(opts.mapsize,3);
    map(fracs(1)+1:fracs(1)+fracs(2),:) = maptmp;
    if opts.coption
        map(1:fracs(1),:) = repmat(maptmp(1,:),[fracs(1) 1]);
        map(fracs(1)+fracs(2)+1:sum(fracs),:) = repmat(maptmp(fracs(2),:),[fracs(3) 1]);
    else
        map(1:fracs(1),:) = repmat(opts.colourgap,[fracs(1) 1]);
        map(fracs(1)+fracs(2)+1:sum(fracs),:) = repmat(opts.colourgap,[fracs(3) 1]);
    end
end
colormap(map)

% Render using the selected layout
switch opts.layout
    
    case 'simple'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.45;
        
        % First panel (right frontal view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(-245,20);
        useopts(opts,fighandle);
        
        % Second panel (left frontal view)
        ax(2) = axes('Position',[0.525 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(245,20);
        useopts(opts,fighandle);
        
        % Deal with the colourbar
        if opts.colourbar
            cb = colorbar('Location','SouthOutside');
            set(cb,'Position',[0.35 0.085 0.3 0.02]);
        end
         
    case 'cardinal'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.5;
        w = 0.5;
        
        % Upper-left panel (left lateral view)
        ax(1) = axes('Position',[0 0.5 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Upper-right panel (right lateral view)
        ax(2) = axes('Position',[0.5 0.5 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Lower-left panel (left medial view)
        ax(3) = axes('Position',[0 0 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Lower-right panel (right medial view)
        ax(4) = axes('Position',[0.5 0 w h]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Deal with the colourbar
        if opts.colourbar
            cb = colorbar('Location','South');
            set(cb,'Position',[0.35 0.475 0.3 0.02]);
        end
        
    case 'left'

        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.45;
        
        % First panel (left lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Second panel (left medial view)
        ax(2) = axes('Position',[0.525 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);

        % This mode has no colourbar
        opts.colourbar = false;
        
    case 'right'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.45;
        
        % First panel (right lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Second panel (right medial view)
        ax(2) = axes('Position',[0.525 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);

        % This mode has no colourbar
        opts.colourbar = false;

    case 'leftlateral'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.95;
        
        % First panel (right lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % This mode has no colourbar
        opts.colourbar = false;
        
    case 'leftmedial'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.95;
        
        % First panel (right lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % This mode has no colourbar
        opts.colourbar = false;
        
    case 'rightlateral'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.95;
        
        % First panel (right lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % This mode has no colourbar
        opts.colourbar = false;
        
    case 'rightmedial'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.95;
        
        % First panel (right lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % This mode has no colourbar
        opts.colourbar = false;

    case 'strip'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.95;
        w = 0.2;
        
        % First panel (left lateral view)
        ax(1) = axes('Position',[0.025 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Second panel (right lateral view)
        ax(2) = axes('Position',[0.275 0.025 w h]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Third panel (left medial view)
        ax(3) = axes('Position',[0.525 0.025 w h]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Fourth panel (right medial view)
        ax(4) = axes('Position',[0.775 0.025 w h]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % This mode has no colourbar
        opts.colourbar = false;

    case 'publication'
        % A more symmetric layout than Worsley's (see below)

        % Auxiliary variables for the sizes of each panel
        h = 0.39;
        w = 0.4;
        xyz_range = max(srf{3}.data.vtx,[],1) - min(srf{3}.data.vtx,[],1);
        wb = h/xyz_range(2)*xyz_range(1)*3/4;
        hb = h/xyz_range(2)*xyz_range(1);
        
        % Upper-left panel (left lateral view)
        ax(1) = axes('Position',[0.055 0.62 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Upper-center panel (top view)
        ax(2) = axes('Position',[0.3 0.58 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,90);
        useopts(opts,fighandle);
        
        % Upper-right panel (right lateral view)
        ax(3) = axes('Position',[1-0.055-h*3/4 0.62 h*3/4 w]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);

        % Middle-left panel (frontal view)
        ax(4) = axes('Position',[0.055 0.345 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(180,0);
        useopts(opts,fighandle);
        
        % Middle-right panel (occipital view)
        ax(5) = axes('Position',[1-0.055-wb 0.345 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,0);
        useopts(opts,fighandle);
        
        % Lower-left panel (left medial view)
        ax(6) = axes('Position',[0.055 -0.015 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Lower-center panel (bottom view)
        ax(7) = axes('Position',[0.3 0.03 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(180,-90);
        useopts(opts,fighandle);
        
        % Lower-right panel (right medial view)
        ax(8) = axes('Position',[1-0.055-h*3/4 -0.015 h*3/4 w]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);

        % This layout always includes a colourbar
        cb = colorbar('Location','South');
        if isempty(opts.title)
            set(cb,'Position',[0.35 0.49 0.3 0.03]);
        else
            set(cb,'Position',[0.35 0.48 0.3 0.03]);
        end
        
    case 'worsley'
        % This is the layout originally proposed by Keith Worsley
        % (McGill University, circa 2005).

        % Auxiliary variables for the sizes of each panel
        h = 0.39;
        w = 0.4;
        xyz_range = max(srf{3}.data.vtx,[],1) - min(srf{3}.data.vtx,[],1);
        wb = h/xyz_range(2)*xyz_range(1)*3/4;
        hb = h/xyz_range(2)*xyz_range(1);
        
        % Upper-left panel (left lateral view)
        ax(1) = axes('Position',[0.055 0.62 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Upper-center panel (top view)
        ax(2) = axes('Position',[0.3 0.58 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,90);
        useopts(opts,fighandle);
        
        % Upper-right panel (right lateral view)
        ax(3) = axes('Position',[1-0.055-h*3/4 0.62 h*3/4 w]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Middle-left panel (left medial view)
        ax(4) = axes('Position',[0.055 0.29 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts,fighandle);
        
        % Middle-center panel (bottom view)
        ax(5) = axes('Position',[0.3 0.18 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,-90);
        useopts(opts,fighandle);
        
        % Middle-right panel (right medial view)
        ax(6) = axes('Position',[1-0.055-h*3/4 0.29 h*3/4 w]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts,fighandle);
        
        % Lower-left panel (frontal view)
        ax(7) = axes('Position',[0.055 0.02 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(180,0);
        useopts(opts,fighandle);
        
        % Lower-right panel (occipital view)
        ax(8) = axes('Position',[1-0.055-wb 0.03 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,0);
        useopts(opts,fighandle);
        
        % This layout always includes a colourbar
        cb = colorbar('Location','South');
        set(cb,'Position',[0.35 0.085 0.3 0.03]);
        
    otherwise
        error('Unknown layout style %s.\n',opts.layout)
end

% Deal with the colour limits
for i=1:length(ax)
    set(ax(i),'CLim',opts.datarange);
end

% Deal with the colourbar
if opts.colourbar
    set(cb,'XAxisLocation','Bottom');
    set(get(cb,'Title'),'String',opts.title);
end

function useopts(opts,fighandle) % ==================================================
daspect([1 1 1]);
axis tight;
axis vis3d off;
box off
if strcmpi(opts.layout,'worsley')
    set(fighandle,'Color',[1 1 1],'InvertHardcopy','Off');
    camlight;
    lighting ('gouraud');
    material ('shiny');
    shading  ('interp');
else
    set(fighandle,'Color',opts.background,'InvertHardcopy','Off');
    if opts.camlight
        camlight;
    end
    lighting (opts.lightning);
    material (opts.material);
    shading  (opts.shading);
end
