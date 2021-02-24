function palm_viewbh()
% work in progress...

% Defaults
opts.layout     = 'basic';
opts.title      = 'Average thickness (mm)';
opts.background = [1 1 1];
opts.camlight   = false;
opts.lightning  = 'phong';
opts.material   = 'shiny';
opts.shading    = 'interp';
opts.colormap   = 'magma';
opts.clim       = [1.6 3.2];
opts.colorbar   = true;

% Ensure we gave the programs needed to read the data
palm_checkprogs;

srf{1} = palm_miscread('/opt/freesurfer/6.0.0/subjects/fsaverage3/surf/lh.pial');
srf{2} = palm_miscread('/opt/freesurfer/6.0.0/subjects/fsaverage3/surf/rh.pial');
dat{1} = palm_miscread('/opt/freesurfer/6.0.0/subjects/fsaverage3/surf/lh.thickness');
dat{2} = palm_miscread('/opt/freesurfer/6.0.0/subjects/fsaverage3/surf/rh.thickness');

% Merge left and right hemispheres
srf{3}.data.vtx = vertcat(srf{1}.data.vtx,srf{2}.data.vtx);
srf{3}.data.fac = vertcat(srf{1}.data.fac,srf{2}.data.fac+size(srf{1}.data.vtx,1));
dat{3}.data = vertcat(dat{1}.data,dat{2}.data);
for h = 1:3
    dat{h}.data = squeeze(dat{h}.data);
end

% Create new or clear the current figure window,
% define the colormap and the background color
clf;
colormap(opts.colormap)
set(gcf,'Color',opts.background,'InvertHardcopy','Off');

switch opts.layout
    
    case 'basic'
        
        % Shows a left antero-lateral view. Ideal for rotation.
        a(1) = axes('Position',[0.13 0.11 0.775 0.815]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(235,15);
        useopts(opts);
        opts.colorbarposition = 'SouthOutside';
         
    case 'cardinal'
        
        % Auxiliary variables for the sizes of each panel
        h = 0.39;
        w = 0.4;
        
        % Upper-left panel (left lateral view)
        a(1) = axes('Position',[0.055 0.62 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts);
        
        % Upper-right panel (right lateral view)
        a(2)=axes('Position',[1-0.055-h*3/4 0.62 h*3/4 w]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts);
        
        % Middle-left panel (left medial view)
        a(3)=axes('Position',[0.055 0.29 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts);
        
        % Middle-right panel (right medial view)
        a(4)=axes('Position',[1-0.055-h*3/4 0.29 h*3/4 w]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts);
        
    case 'left'
        
    case 'right'
        
    case 'strip'
        
    case 'worsley'
        % This is the layout originally proposed by Keith Worsley
        % (McGill University, circa 2005).
        
        % This layout always include a colorbar
        opts.colorbar = true;
        opts.colorbarposition = 'South';
        
        % Auxiliary variables for the sizes of each panel
        h = 0.39;
        w = 0.4;
        r  = max(srf{3}.data.vtx,[],1)-min(srf{3}.data.vtx,[],1);
        wb = h/r(2)*r(1)*3/4;
        hb = h/r(2)*r(1);
        
        % Upper-left panel (left lateral view)
        a(1) = axes('Position',[0.055 0.62 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(-90,0);
        useopts(opts);
        
        % Upper-center panel (top view)
        a(2) = axes('Position',[0.3 0.58 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,90);
        useopts(opts);
        
        % Upper-right panel (right lateral view)
        a(3) = axes('Position',[1-0.055-h*3/4 0.62 h*3/4 w]);
        trisurf(...
            srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts);
        
        % Middle-left panel (left medial view)
        a(4) = axes('Position',[0.055 0.29 h*3/4 w]);
        trisurf(...
            srf{1}.data.fac,...
            srf{1}.data.vtx(:,1),srf{1}.data.vtx(:,2),srf{1}.data.vtx(:,3),...
            dat{1}.data,...
            'EdgeColor','None');
        view(90,0);
        useopts(opts);
        
        % Middle-center panel (bottom view)
        a(5) = axes('Position',[0.3 0.18 w h]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,-90);
        useopts(opts);
        
        % Middle-right panel (right medial view)
        a(6) = axes('Position',[1-0.055-h*3/4 0.29 h*3/4 w]);
        trisurf(srf{2}.data.fac,...
            srf{2}.data.vtx(:,1),srf{2}.data.vtx(:,2),srf{2}.data.vtx(:,3),...
            dat{2}.data,'EdgeColor','None');
        view(-90,0);
        useopts(opts);
        
        % Lower-left panel (frontal view)
        a(7) = axes('Position',[0.055 0.02 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(180,0);
        useopts(opts);
        
        % Lower-right panel (occipital view)
        a(8) = axes('Position',[1-0.055-wb 0.03 wb hb]);
        trisurf(...
            srf{3}.data.fac,...
            srf{3}.data.vtx(:,1),srf{3}.data.vtx(:,2),srf{3}.data.vtx(:,3),...
            dat{3}.data,...
            'EdgeColor','None');
        view(0,0);
        useopts(opts);
end

% Deal with the color limits
for i=1:length(a)
    set(a(i),'CLim',opts.clim);
end

% Deal with the colorbar
if opts.colorbar
    cb = colorbar('Location',opts.colorbarposition);
    set(cb,'Position',[0.35 0.085 0.3 0.03]);
    set(cb,'XAxisLocation','Bottom');
    set(get(cb,'Title'),'String',opts.title);
end

function useopts(opts) % ==================================================
daspect([1 1 1]);
axis tight;
axis vis3d off;
if opts.camlight
    camlight;
end
lighting(opts.lightning);
material(opts.material);
shading(opts.shading);
