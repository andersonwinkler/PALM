function palm_saveall(plm,opts)
% Save most of the outputs from PALM.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
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

% For the MV and NPC methods in which the most significant stats are the
% smallest, rather than the largest, use reverse comparisons.
if opts.NPC,
    if plm.npcrev,
        npcextr = @min;
    else
        npcextr = @max;
    end
end

fprintf('Computing p-values.\n');
% Start with the uncorrected, but don't save them yet.
% They'll be used later for the FDR.
for y = 1:plm.nY,
    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
    for m = loopM,
        for c = 1:plm.nC(m),
            if opts.draft,
                plm.Gpperm{y}{m}{c} = (plm.Gpperm{y}{m}{c} + 1)./plm.Gppermp{y}{m}{c};
            else
                plm.Gpperm{y}{m}{c} = plm.Gpperm{y}{m}{c}/plm.nP{m}(c);
                if opts.tfce_uni.do,
                    plm.Gtfcepperm{y}{m}{c} = plm.Gtfcepperm{y}{m}{c}/plm.nP{m}(c);
                end
            end
        end
    end
end
if opts.NPC,
    if opts.npcmod && ~ opts.npccon,
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                plm.Tpperm{m}{c} = plm.Tpperm{m}{c}/numel(plm.Tmax{m}{c});
                if opts.tfce_npc.do,
                    plm.Ttfcepperm{m}{c} = plm.Ttfcepperm{m}{c}/plm.nP{m}(c);
                end
            end
        end
    elseif opts.npccon,
        for j = 1:numel(plm.Tmax),
            plm.Tpperm{j} = plm.Tpperm{j}/numel(plm.Tmax{j});
            if opts.tfce_npc.do,
                plm.Ttfcepperm{j} = plm.Ttfcepperm{j}/numel(plm.Tmax{j});
            end
        end
    end
end
if opts.MV || opts.CCA,
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            if opts.draft,
                plm.Qpperm{m}{c} = (plm.Qpperm{m}{c} + 1)./plm.Qppermp{m}{c};
            else
                plm.Qpperm{m}{c} = plm.Qpperm{m}{c}/plm.nP{m}(c);
            end
            if opts.tfce_mv.do,
                plm.Qtfcepperm{m}{c} = plm.Qtfcepperm{m}{c}/plm.nP{m}(c);
            end
        end
    end
end

% Save uncorrected & FWER-corrected within modality for this contrast.
if opts.saveunivariate,
    fprintf('Saving p-values (uncorrected, and corrected within modality and within contrast).\n');
    for y = 1:plm.nY,
        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
        for m = loopM,
            for c = 1:plm.nC(m),
                
                % Only permutation p-value and its FDR ajustment are saved in the draft mode.
                if opts.draft,
                    
                    % Permutation p-value, uncorrected
                    palm_quicksave(plm.Gpperm{y}{m}{c},1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    
                    % Permutation p-value, FDR adjusted
                    if opts.FDR,
                        palm_quicksave(fastfdr(plm.Gpperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                else
                    
                    % Permutation p-value
                    palm_quicksave(plm.Gpperm{y}{m}{c},1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    
                    % FWER-corrected
                    palm_quicksave( ...
                        palm_datapval(plm.G{y}{m}{c},plm.Gmax{y}{m}{c},false),1,opts,plm,y,m,c,...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    
                    % Permutation p-value, FDR adjusted
                    if opts.FDR,
                        palm_quicksave(fastfdr(plm.Gpperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                    
                    % Cluster extent results.
                    if opts.clustere_uni.do,
                        
                        % Cluster extent statistic.
                        palm_quicksave(plm.Gcle{y}{m}{c},0,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        
                        % Cluster extent FWER p-value
                        palm_quicksave( ...
                            palm_datapval(plm.Gcle{y}{m}{c},plm.Gclemax{y}{m}{c},false),1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                    
                    % Cluster mass results.
                    if opts.clusterm_uni.do,
                        
                        % Cluster mass statistic.
                        palm_quicksave(plm.Gclm{y}{m}{c},0,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        
                        % Cluster mass FWER p-value.
                        palm_quicksave( ...
                            palm_datapval(plm.Gclm{y}{m}{c},plm.Gclmmax{y}{m}{c},false),1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                    
                    % TFCE results
                    if opts.tfce_uni.do,
                        
                        % TFCE statistic
                        palm_quicksave(plm.Gtfce{y}{m}{c},0,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        
                        % TFCE p-value
                        palm_quicksave(plm.Gtfcepperm{y}{m}{c},1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        
                        % TFCE FWER-corrected within modality and contrast.
                        palm_quicksave( ...
                            palm_datapval(plm.Gtfce{y}{m}{c},plm.Gtfcemax{y}{m}{c},false),1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        
                        % TFCE p-value, FDR adjusted.
                        if opts.FDR,
                            palm_quicksave(fastfdr(plm.Gtfcepperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
                
                % Parametric p-value and its FDR adjustment
                if opts.savepara,
                    P = palm_quicksave(plm.G{y}{m}{c},2,opts,plm,y,m,c,...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncparap',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    if opts.FDR,
                        palm_quicksave(fastfdr(P),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrparap',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across modalities.
    if opts.corrmod,
        fprintf('Saving p-values (corrected across modalities).\n')
        
        % FWER correction (non-spatial stats)
        if opts.designperinput,
            for c = 1:plm.nC(1),
                distmax = zeros(plm.nP{1}(c),plm.nY);
                for y = 1:plm.nY,
                    m = y;
                    distmax(:,y) = plm.Gmax{y}{m}{c};
                end
                distmax = max(distmax,[],2);
                for y = 1:plm.nY,
                    m = y;
                    palm_quicksave( ...
                        palm_datapval(plm.G{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                end
            end
        else
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    distmax = zeros(plm.nP{m}(c),plm.nY);
                    for y = 1:plm.nY,
                        distmax(:,y) = plm.Gmax{y}{m}{c};
                    end
                    distmax = max(distmax,[],2);
                    for y = 1:plm.nY,
                        palm_quicksave( ...
                            palm_datapval(plm.G{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    pmerged = zeros(sum(plm.Ysiz),1);
                    for y = 1:plm.nY,
                        m = y;
                        pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                    end
                    pfdradj = fastfdr(pmerged);
                    for y = 1:plm.nY,
                        m = y;
                        palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        pmerged = zeros(sum(plm.Ysiz),1);
                        for y = 1:plm.nY,
                            pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                        end
                        pfdradj = fastfdr(pmerged);
                        for y = 1:plm.nY,
                            palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
            end
        end
        
        % Cluster extent
        if opts.clustere_uni.do && ...
                (all(plm.Yisvol) || all(plm.Yissrf)),
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    distmax = zeros(plm.nP{1}(c),plm.nY);
                    for y = 1:plm.nY,
                        m = y;
                        distmax(:,y) = plm.Gclemax{y}{m}{c};
                    end
                    distmax = max(distmax,[],2);
                    for y = 1:plm.nY,
                        m = y;
                        palm_quicksave( ...
                            palm_datapval(plm.Gcle{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        distmax = zeros(plm.nP{m}(c),plm.nY);
                        for y = 1:plm.nY,
                            distmax(:,y) = plm.Gclemax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            palm_quicksave( ...
                                palm_datapval(plm.Gcle{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
            end
        end
        
        % Cluster mass
        if opts.clusterm_uni.do && ...
                (all(plm.Yisvol) || all(plm.Yissrf)),
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    distmax = zeros(plm.nP{1}(c),plm.nY);
                    for y = 1:plm.nY,
                        m = y;
                        distmax(:,y) = plm.Gclmmax{y}{m}{c};
                    end
                    distmax = max(distmax,[],2);
                    for y = 1:plm.nY,
                        m = y;
                        palm_quicksave( ...
                            palm_datapval(plm.Gclm{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        distmax = zeros(plm.nP{m}(c),plm.nY);
                        for y = 1:plm.nY,
                            distmax(:,y) = plm.Gclmmax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            palm_quicksave( ...
                                palm_datapval(plm.Gclm{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
            end
        end
        
        % TFCE
        if opts.tfce_uni.do && ...
                (all(plm.Yisvol) || all(plm.Yissrf)),
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    distmax = zeros(plm.nP{1}(c),plm.nY);
                    for y = 1:plm.nY,
                        m = y;
                        distmax(:,y) = plm.Gtfcemax{y}{m}{c};
                    end
                    distmax = max(distmax,[],2);
                    for y = 1:plm.nY,
                        m = y;
                        palm_quicksave( ...
                            palm_datapval(plm.Gtfce{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        distmax = zeros(plm.nP{m}(c),plm.nY);
                        for y = 1:plm.nY,
                            distmax(:,y) = plm.Gtfcemax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            palm_quicksave( ...
                                palm_datapval(plm.Gtfce{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
            end
            if opts.FDR,
                if opts.designperinput,
                    for c = 1:plm.nC(1),
                        pmerged = zeros(sum(plm.Ysiz),1);
                        for y = 1:plm.nY,
                            m = y;
                            pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                        end
                        pfdradj = fastfdr(pmerged);
                        for y = 1:plm.nY,
                            m = y;
                            palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                else
                    for m = 1:plm.nM,
                        for c = 1:plm.nC(m),
                            pmerged = zeros(sum(plm.Ysiz),1);
                            for y = 1:plm.nY,
                                pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                            end
                            pfdradj = fastfdr(pmerged);
                            for y = 1:plm.nY,
                                palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across contrasts.
    if opts.corrcon,
        fprintf('Saving p-values (corrected across contrasts).\n');
        
        % FWER correction (non-spatial stats)
        for y = 1:plm.nY,
            if opts.designperinput,
                loopM = y;
                distmax = zeros(plm.nP{1}(1),plm.nC(1));
            else
                loopM = 1:plm.nM;
                distmax = zeros(plm.nP{1}(1),sum(plm.nC));
            end
            j = 1;
            for m = loopM,
                for c = 1:plm.nC(m),
                    distmax(:,j) = plm.Gmax{y}{m}{c};
                    j = j + 1;
                end
            end
            distmax = max(distmax,[],2);
            for m = loopM,
                for c = 1:plm.nC(m),
                    palm_quicksave( ...
                        palm_datapval(plm.G{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                end
            end
        end
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    pmerged = zeros(plm.nC(1),plm.Ysiz(y));
                else
                    loopM = 1:plm.nM;
                    pmerged = zeros(sum(plm.nC),plm.Ysiz(y));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        pmerged(j,:) = plm.Gpperm{y}{m}{c};
                        j = j + 1;
                    end
                end
                pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave(pfdradj(j,:),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_cfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        j = j + 1;
                    end
                end
            end
        end
        
        % Cluster extent
        if opts.clustere_uni.do,
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    distmax = zeros(plm.nP{1}(1),plm.nC(1));
                else
                    loopM = 1:plm.nM;
                    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gclemax{y}{m}{c};
                        j = j + 1;
                    end
                end
                distmax = max(distmax,[],2);
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gcle{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % Cluster mass
        if opts.clusterm_uni.do,
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    distmax = zeros(plm.nP{1}(1),plm.nC(1));
                else
                    loopM = 1:plm.nM;
                    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gclmmax{y}{m}{c};
                        j = j + 1;
                    end
                end
                distmax = max(distmax,[],2);
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gclm{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % TFCE
        if opts.tfce_uni.do,
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    distmax = zeros(plm.nP{1}(1),plm.nC(1));
                else
                    loopM = 1:plm.nM;
                    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gtfcemax{y}{m}{c};
                        j = j + 1;
                    end
                end
                distmax = max(distmax,[],2);
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gtfce{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
            if opts.FDR,
                for y = 1:plm.nY,
                    if opts.designperinput,
                        loopM = y;
                        pmerged = zeros(plm.nC(1),plm.Ysiz(y));
                    else
                        loopM = 1:plm.nM;
                        pmerged = zeros(sum(plm.nC),plm.Ysiz(y));
                    end
                    j = 1;
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            pmerged(j,:) = plm.Gtfcepperm{y}{m}{c};
                            j = j + 1;
                        end
                    end
                    pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                    j = 1;
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            palm_quicksave(pfdradj(j,:),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_cfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                            j = j + 1;
                        end
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across modalities and contrasts.
    if opts.corrmod && opts.corrcon,
        fprintf('Saving p-values (corrected across modalities and contrasts).\n')
        
        % FWER correction (non-spatial stats)
        if opts.designperinput,
            distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
        else
            distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
        end
        j = 1;
        for y = 1:plm.nY,
            if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
            for m = loopM,
                for c = 1:plm.nC(m),
                    distmax(:,j) = plm.Gmax{y}{m}{c};
                    j = j + 1;
                end
            end
        end
        distmax = max(distmax,[],2);
        for y = 1:plm.nY,
            if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
            for m = loopM,
                for c = 1:plm.nC(m),
                    palm_quicksave( ...
                        palm_datapval(plm.G{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                end
            end
        end
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            if opts.designperinput,
                pmerged = zeros(plm.nC(1),sum(plm.Ysiz));
            else
                pmerged = zeros(sum(plm.nC),sum(plm.Ysiz));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        pmerged(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave(pfdradj(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mcfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % Cluster extent
        if opts.clustere_uni.do,
            if opts.designperinput,
                distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
            else
                distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gclemax{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            distmax = max(distmax,[],2);
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gcle{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clustere',plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % Cluster mass
        if opts.clusterm_uni.do,
            if opts.designperinput,
                distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
            else
                distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gclmmax{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            distmax = max(distmax,[],2);
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gclm{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_clusterm',plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
        end
        
        % TFCE
        if opts.tfce_uni.do,
            if opts.designperinput,
                distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
            else
                distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gtfcemax{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            distmax = max(distmax,[],2);
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave( ...
                            palm_datapval(plm.Gtfce{y}{m}{c},distmax,false),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                    end
                end
            end
            if opts.FDR,
                if opts.designperinput,
                    pmerged = zeros(plm.nC(1),sum(plm.Ysiz));
                else
                    pmerged = zeros(sum(plm.nC),sum(plm.Ysiz));
                end
                j = 1;
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            pmerged(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                            j = j + 1;
                        end
                    end
                end
                pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            palm_quicksave(pfdradj(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,'_tfce',plm.Gname{m}{c},'_mcfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                        end
                    end
                end
            end
        end
    end
end

% Save NPC between modalities, corrected within contrasts
if opts.npcmod && ~ opts.npccon,
    fprintf('Saving p-values for NPC between modalities (uncorrected and corrected within contrasts).\n');
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),

            % NPC p-value
            palm_quicksave(plm.Tpperm{m}{c},1,opts,plm,[],m,c, ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncp',plm.mstr{m},plm.cstr{c}));
            
            % NPC FWER-corrected
            palm_quicksave( ...
                palm_datapval(plm.T{m}{c},plm.Tmax{m}{c},plm.npcrev),1,opts,plm,[],m,c, ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fwep',plm.mstr{m},plm.cstr{c}));
            
            % NPC FDR
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Tpperm{m}{c}),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fdrp',plm.mstr{m},plm.cstr{c}));
            end
            
            % Parametric combined pvalue
            if opts.savepara && ~ plm.nonpcppara,
                palm_quicksave(plm.Tppara{m}{c},1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncparap',plm.mstr{m},plm.cstr{c}));
            end
            
            % Cluster extent NPC results.
            if opts.clustere_npc.do,
                
                % Cluster extent statistic.
                palm_quicksave(plm.Tcle{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,plm.mstr{m},plm.cstr{c}));
                
                % Cluster extent FWER p-value
                palm_quicksave( ...
                    palm_datapval(plm.Tcle{m}{c},plm.Tclemax{m}{c},false),1,opts,plm,y,m,c,...
                    sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,'_fwep',plm.mstr{m},plm.cstr{c}));
            end
            
            % Cluster mass NPC results.
            if opts.clusterm_npc.do,
                
                % Cluster mass statistic.
                palm_quicksave(plm.Tclm{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,plm.mstr{m},plm.cstr{c}));
                
                % Cluster mass FWER p-value
                palm_quicksave( ...
                    palm_datapval(plm.Tclm{m}{c},plm.Tclmmax{m}{c},false),1,opts,plm,y,m,c,...
                    sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,'_fwep',plm.mstr{m},plm.cstr{c}));
            end
            
            % TFCE NPC results.
            if opts.tfce_npc.do,
                
                % TFCE statistic.
                palm_quicksave(plm.Ttfce{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,plm.mstr{m},plm.cstr{c}));
                
                % TFCE p-value
                palm_quicksave(plm.Ttfcepperm{m}{c},1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_uncp',plm.mstr{m},plm.cstr{c}));
                
                % TFCE FWER p-value
                palm_quicksave( ...
                    palm_datapval(plm.Ttfce{m}{c},plm.Ttfcemax{m}{c},false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_fwep',plm.mstr{m},plm.cstr{c}));
                
                % TFCE FDR p-value
                if opts.FDR,
                    palm_quicksave(fastfdr(plm.Ttfcepperm{m}{c}),1,opts,plm,[],m,c,...
                        sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_uncp',plm.mstr{m},plm.cstr{c}));
                end
            end
        end
    end
end

% Save NPC between modalities, corrected across contrasts
if opts.npcmod && ~ opts.npccon && opts.corrcon,
    fprintf('Saving p-values for NPC between modalities (corrected across contrasts).\n');
    
    % FWER correction (non-spatial stats)
    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
    j = 1;
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            distmax(:,j) = plm.Tmax{m}{c};
            j = j + 1;
        end
    end
    distmax = max(distmax,[],2);
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            palm_quicksave( ...
                palm_datapval(plm.T{m}{c},distmax,plm.npcrev),1,opts,plm,[],m,c,...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{c}));
        end
    end
    
    % FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                pmerged(:,j) = plm.Tpperm{m}{c};
                j = j + 1;
            end
        end
        pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave(pfdradj(:,j),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_cfdrp',plm.mstr{m},plm.cstr{c}));
                j = j + 1;
            end
        end
    end
    
    % Cluster extent NPC
    if opts.clustere_npc.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Tclemax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Tcle{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
    end
    
    % Cluster mass NPC
    if opts.clustere_npc.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Tclmmax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Tclm{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
    end
    
    % TFCE NPC
    if opts.tfce_npc.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Ttfcemax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Ttfce{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
        if opts.FDR,
            pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    pmerged(:,j) = plm.Ttfcepperm{m}{c};
                    j = j + 1;
                end
            end
            pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    palm_quicksave(pfdradj(:,j),1,opts,plm,[],m,c, ...
                        sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_cfdrp',plm.mstr{m},plm.cstr{c}));
                    j = j + 1;
                end
            end
        end
    end
end

% Save NPC between contrasts, corrected within modality
if opts.npccon,
    fprintf('Saving p-values for NPC between contrasts (uncorrected and corrected within modality).\n');
    
    for j = 1:numel(plm.Tmax),
        
        % NPC p-value
        palm_quicksave(plm.Tpperm{j},1,opts,plm,j,[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncp',plm.jstr{j}));
        
        % NPC FWER-corrected
        palm_quicksave( ...
            palm_datapval(plm.T{j},plm.Tmax{j},plm.npcrev),1,opts,plm,j,[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
        
        % NPC FDR
        if opts.FDR,
            palm_quicksave(fastfdr(plm.Tpperm{j}),1,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fdrp',plm.jstr{j}));
        end
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara,
            palm_quicksave(plm.Tppara{j},1,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncparap',plm.jstr{j}));
        end
        
        % Cluster extent NPC results.
        if opts.clustere_npc.do,
            
            % Cluster extent statistic.
            palm_quicksave(plm.Tcle{j},0,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,plm.jstr{j}));
            
            % Cluster extent FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Tcle{j},plm.Tclemax{j},false),1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
        end
        
        % Cluster mass NPC results.
        if opts.clusterm_npc.do,
            
            % Cluster mass statistic.
            palm_quicksave(plm.Tclm{j},0,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,plm.jstr{j}));
            
            % Cluster mass FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Tclm{j},plm.Tclmmax{j},false),1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
        end
        
        % TFCE NPC results.
        if opts.tfce_npc.do,
            
            % TFCE statistic.
            palm_quicksave(plm.Ttfce{j},0,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,plm.jstr{j}));
            
            % TFCE p-value
            palm_quicksave(plm.Ttfcepperm{j},1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_uncp',plm.jstr{j}));
            
            % TFCE FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Ttfce{j},plm.Ttfcemax{j},false),1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
            
            % TFCE FDR p-value
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Ttfcepperm{j}),1,opts,plm,j,[],[], ...
                    sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_fdrp',plm.jstr{j}));
            end
        end
    end
end

% Save the NPC over contrasts, corrected for modalities
if ~ opts.npcmod && opts.npccon && opts.corrmod,
    fprintf('Saving p-values for NPC over contrasts (corrected across modalities).\n')
    
    % NPC FWER-corrected across modalities.
    distmax = npcextr(cat(2,plm.Tmax{:}),2);
    for y = 1:numel(plm.nY),
        palm_quicksave( ...
            palm_datapval(plm.T{y},distmax,plm.npcrev),1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara,
            palm_quicksave(plm.Tppara{y},1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_uncparap',plm.ystr{y}));
        end
    end
    
    % NPC FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.Ysiz),1);
        for y = 1:plm.nY,
            pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Tpperm{y};
        end
        pfdradj = fastfdr(pmerged);
        for y = 1:plm.nY,
            palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.Tname,'_mfdrp',plm.ystr{y}));
        end
    end
    
    % NPC FDR-corrected across modalities.
    distmax = npcextr(cat(2,plm.Tmax{:}),2);
    for y = 1:numel(plm.nY),
        palm_quicksave( ...
            palm_datapval(plm.T{y},distmax,plm.npcrev),1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara,
            palm_quicksave(plm.Tppara{y},1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_uncparap',plm.ystr{y}));
        end
    end
    
    % Cluster extent NPC results.
    if opts.clustere_npc.do,
        distmax = npcextr(cat(2,plm.Tclemax{:}),2);
        for y = 1:numel(plm.nY),
            
            % Cluster extent statistic.
            palm_quicksave(plm.Tcle{y},0,opts,plm,y,[],[], ...
                sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,plm.ystr{y}));
            
            % Cluster extent FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Tcle{y},distmax,false),1,opts,plm,y,[],[],...
                sprintf('%s',opts.o,'_clustere',plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        end
    end
    
    % Cluster mass NPC results.
    if opts.clusterm_npc.do,
        distmax = npcextr(cat(2,plm.Tclmmax{:}),2);
        for y = 1:numel(plm.nY),
            
            % Cluster extent statistic.
            palm_quicksave(plm.Tclm{y},0,opts,plm,y,[],[], ...
                sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,plm.ystr{y}));
            
            % Cluster extent FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Tclm{y},distmax,false),1,opts,plm,y,[],[],...
                sprintf('%s',opts.o,'_clusterm',plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        end
    end
    
    % TFCE NPC results.
    if opts.tfce_npc.do,
        distmax = npcextr(cat(2,plm.Ttfce{:}),2);
        for y = 1:numel(plm.nY),
            
            % Cluster extent statistic.
            palm_quicksave(plm.Ttfce{y},0,opts,plm,y,[],[], ...
                sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,plm.ystr{y}));
            
            % Cluster extent FWER p-value
            palm_quicksave( ...
                palm_datapval(plm.Ttfce{y},distmax,false),1,opts,plm,y,[],[],...
                sprintf('%s',opts.o,'_tfce',plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        end
        
        % NPC FDR correction TFCE
        if opts.FDR,
            pmerged = zeros(sum(plm.Ysiz),1);
            for y = 1:plm.nY,
                pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Ttfcepperm{y};
            end
            pfdradj = fastfdr(pmerged);
            for y = 1:plm.nY,
                palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,'_tfce',plm.Tname,'_mfdrp',plm.ystr{y}));
            end
        end
    end
end

% Save the MV results for each contrast
if opts.MV || opts.CCA,
    fprintf('Saving p-values for classical multivariate (uncorrected and corrected within contrast).\n')
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            % MV p-value
            palm_quicksave(plm.Qpperm{m}{c},1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_uncp',plm.mstr{m},plm.cstr{c}));
            
            % MV FWER-corrected within modality and contrast.
            palm_quicksave( ...
                palm_datapval(plm.Q{m}{c},plm.Qmax{m}{c},plm.mvrev{m}{c}),1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{c}));
            
            % MV FDR
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Qpperm{m}{c}),1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_fdrp',plm.mstr{m},plm.cstr{c}));
            end
            
            % Parametric MV pvalue
            if opts.savepara && ~ plm.nomvppara,
                palm_quicksave(plm.Qppara{m}{c},1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_uncparap',plm.mstr{m},plm.cstr{c}));
            end
            
            % Cluster extent MV results.
            if opts.clustere_mv.do,
                
                % Cluster extent statistic.
                palm_quicksave(plm.Qcle{m}{c},0,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,'_clustere',plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}));
                
                % Cluster extent FWER p-value
                palm_quicksave( ...
                    palm_datapval(plm.Qcle{m}{c},plm.Qclemax{m}{c},false),1,opts,plm,[],[],[],...
                    sprintf('%s',opts.o,'_clustere',plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{c}));
            end
            
            % Cluster mass MV results.
            if opts.clusterm_mv.do,
                
                % Cluster mass statistic.
                palm_quicksave(plm.Qclm{m}{c},0,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,'_clusterm',plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}));
                
                % Cluster mass FWER p-value
                palm_quicksave( ...
                    palm_datapval(plm.Qclm{m}{c},plm.Qclmmax{m}{c},false),1,opts,plm,[],[],[],...
                    sprintf('%s',opts.o,'_clusterm',plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{c}));
            end
            
            % TFCE MV results.
            if opts.tfce_mv.do,
                
                % TFCE statistic.
                palm_quicksave(plm.Qtfce{m}{c},0,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}));
                
                % TFCE p-value
                palm_quicksave(plm.Qtfcepperm{m}{c},1,opts,plm,[],[],[],...
                    sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname,'_uncp',plm.mstr{m},plm.cstr{c}));
                
                % TFCE FWER p-value
                palm_quicksave(palm_datapval( ...
                    plm.Qtfce{m}{c},plm.Qtfcemax{m}{c},false),1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{c}));
                
                % TFCE MV FDR
                if opts.FDR,
                    palm_quicksave(fastfdr(plm.Qtfcepperm{m}{c}),1,opts,plm,[],[],[], ...
                        sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname{m}{c},'_fdrp',plm.mstr{m},plm.cstr{c}));
                end
            end
        end
    end
end

% Save FWER corrected across contrasts for MV.
if ( opts.MV  || opts.CCA ) && opts.corrcon,
    fprintf('Saving p-values for MANOVA/MANCOVA (corrected across contrasts).\n')
    
    % FWER correction (non-spatial stats)
    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
    j = 1;
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            distmax(:,j) = plm.Qmax{m}{c};
            j = j + 1;
        end
    end
    if plm.mvrev{m}{c}, mvextr = @min; else mvextr = @max; end % FIXME FIXME FIXME
    distmax = mvextr(distmax,[],2);
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            palm_quicksave( ...
                palm_datapval(plm.Q{m}{c},distmax,plm.mvrev{m}{c}),1,opts,plm,[],m,c,...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{c}));
        end
    end
    
    % FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                pmerged(:,j) = plm.Qpperm{m}{c};
                j = j + 1;
            end
        end
        pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave(pfdradj(:,j),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Tname,'_cfdrp',plm.mstr{m},plm.cstr{c}));
                j = j + 1;
            end
        end
    end
    
    % Cluster extent MV
    if opts.clustere_mv.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Qclemax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Qcle{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_clustere',plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
    end
    
    % Cluster mass MV
    if opts.clusterm_mv.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Qclmmax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Qclm{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_clusterm',plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
    end
    
    % TFCE MV
    if opts.tfce_mv.do,
        
        % FWER correction
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Qtfcemax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave( ...
                    palm_datapval(plm.Qtfce{m}{c},distmax,false),1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{c}));
            end
        end
        
        % FDR correction TFCE
        if opts.FDR,
            pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    pmerged(:,j) = plm.Qtfcepperm{m}{c};
                    j = j + 1;
                end
            end
            pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    palm_quicksave(pfdradj(:,j),1,opts,plm,[],m,c, ...
                        sprintf('%s',opts.o,'_tfce',plm.mvstr,plm.Qname{m}{c},'_cfdrp',plm.mstr{m},plm.cstr{c}));
                    j = j + 1;
                end
            end
        end
    end
end

% ==============================================================
function padj = fastfdr(pval)
% Compute FDR-adjusted p-values

V = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR]   = sort(oidx);
padj = zeros(size(pval));
prev = 1;
for i = V:-1:1,
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);
