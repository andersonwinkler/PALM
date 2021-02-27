function palm_colourbars()
% This function loops over the m-files in this directory (assumed to be
% all colour maps) and generates a .png file for each map, which can be
% used for easy choice of a map or for composition in other software
% such as Inkscape or GIMP.
% Run it without arguments.
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Feb/2021
% http://brainder.org

[pth,nam,~] = fileparts(mfilename('fullpath'));
F = dir(fullfile(pth,'*.m'));

for f = numel(F):-1:1
    if F(f).name(1) == '.' || any(strcmpi(F(f).name,{[nam '.m'],[nam '.m~'],[nam '.asv']}))
        F(f) = [];
    end
end
for f = 1:numel(F)
    [~,mapname,~] = fileparts(F(f).name);
    cmap2png(mapname,800,[50 800],[mapname '.png'])
end


function cmap2png(mapname,mapsize,imagesize,imagename) % ======================
% Create a PNG image of a given colour map.

%Create the colour map
map = eval(sprintf('%s(%g)',mapname,mapsize));

% Define the indices, with the gradient running along the largest dimension
mapsize = size(map,1);
idx = round(repmat(linspace(1,mapsize,max(imagesize)),[min(imagesize) 1]));

% Transpose to vertical if needed
if imagesize(1) > imagesize(2)
    idx = idx';
end

% Create the PNG and assign the RGB triplets
png = zeros([size(idx) 3]);
for c = 1:3
    png(:,:,c) = reshape(map(idx(:),c),size(idx));
end

% Save to the disk
imwrite(png,imagename);