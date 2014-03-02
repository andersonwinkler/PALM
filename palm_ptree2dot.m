function palm_ptree2dot(Ptree,dotfile)
% Create a DOT file from a permutation tree, that can be used
% with GraphViz for visualisation.
%
% Usage:
% palm_ptree2dot(Ptree,dotfile)
%
% Ptree     : Permutation (dependence) tree.
% dotfile   : DOT file to be created.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
% http://brainder.org

graphname = 'Ptree';
ntop = 'x1';
fid = fopen(dotfile,'w');
fprintf(fid,'graph %s {\n',graphname);
if size(Ptree,2) > 1,
    if isnan(Ptree{1,1}(1)),
        nclr = 'red';
    else
        nclr = 'blue';
    end
    fprintf(fid,'%s [label="" shape=point color=%s fixedsize=true width=.5 height=.5 fontsize=20 penwidth=3];\n',ntop,nclr);
    downstream(fid,ntop,Ptree{1,3});
end
fprintf(fid,'}');
fclose(fid);

% ==============================================================
function downstream(fid,ntop,Ptree)
for u = 1:size(Ptree,1),
    
    % Current node
    ncur    = sprintf('%sx%s',ntop,num2str(u));
    
    % Print node
    inancur = isnan(Ptree{u,1}(1));
    if size(Ptree,2) == 1,
        nlab = sprintf('"%d"',Ptree{u,1});
        nshp = 'circle';
        nclr = 'black';
        nwid = '1';
        nhei = '1';
    elseif size(Ptree{u,3},1) == 1,
        nlab = '""';
        nshp = 'point';
        nclr = 'black';
        nwid = '.2';
        nhei = '.2';
    else
        nlab = '""';
        nshp = 'point';
        if inancur,
            nclr = 'red';
        else
            nclr = 'blue';
        end
        nwid = '1';
        nhei = '1';
    end
    fprintf(fid,'%s [label=%s shape=%s color=%s fixedsize=true width=%s height=%s fontsize=20 penwidth=3];\n',...
        ncur,nlab,nshp,nclr,nwid,nhei);
    
    % Print edge
    fprintf(fid,'%s -- %s [color=black penwidth=6];\n',ntop,ncur);
    
    % Keep going down more levels as needed
    if size(Ptree,2) > 1,
        downstream(fid,ncur,Ptree{u,3});
    end
end
