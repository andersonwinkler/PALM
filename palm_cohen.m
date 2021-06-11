function palm_cohen(M,psi,res,plm,opts,y,m,c,o,p)
% Computes and saves effect sizes and a few other GLM estimates.
% This is intended to be invoked if "-saveglm" is used.
% If "-saveperms" is also used, then it will be executed for every permutation.
% 
% This function isn't intended to be executed by the user, but be called
% from palm_core.m.
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health (NIH)
% Jun/2021
% http://brainder.org
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2021 Anderson M. Winkler
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

if opts.saveperms
    pstr = sprintf('_perm%06d',p);
else
    pstr = '';
end

% Fork for t or F contrasts
if plm.rC0{m}(c) == 1
    cope  = plm.eC{y}{m}{c}{o}'*psi;
    Xr    = range(vertcat(M*plm.eC{y}{m}{c}{o},0),1);
    sigsq = sum(res.^2,1)./(plm.N-plm.rM{y}{m}{c}{o});
    cohen = cope*Xr./sigsq.^.5;
    cfvar = 1./cohen;
    if opts.evperdat
        Theta = zeros(1,size(psi,2));
        for t = 1:size(psi,2)
            Theta(t) = plm.mrdiv(plm.eC{y}{m}{c}{o}',(M(:,:,t)'*M(:,:,t)))*plm.eC{y}{m}{c}{o};
        end
        varcope = Theta ./ sigsq;
    else
        varcope = plm.mrdiv(plm.eC{y}{m}{c}{o}',(M'*M))*plm.eC{y}{m}{c}{o} * sigsq;
    end
    palm_quicksave(cope,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cope',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(varcope,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_varcope',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(cohen,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cohen',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(cfvar,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cfvar',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
else
    error('The option -saveglm hasn''t yet been implemented for F-tests.')
    cope  = plm.eC{y}{m}{c}{o}'*psi;
    if opts.evperdat
        Theta = zeros(1,size(psi,2));
        for t = 1:size(psi,2)
            Theta(t) = pinv(plm.mrdiv(plm.eC{y}{m}{c}{o}',(M(:,:,t)'*M(:,:,t)))*plm.eC{y}{m}{c}{o});
        end
        varcope = Theta .* sigsq;
    else
        % varcope = plm.mrdiv(plm.eC{y}{m}{c}{o}',(M'*M))*plm.eC{y}{m}{c}{o} * sigsq;
        varcope = zeros(size(cope));
    end
    
    
    
    Xr    = range(vertcat(M*plm.eC{y}{m}{c}{o},zeros(1,size(plm.eC{y}{m}{c}{o},2))),1);
    cope  = sqrt(mean(cope.^2,1));
    sigsq = sum(res.^2,1)./(plm.N-plm.rM{y}{m}{c}{o});
    cohen = cope./sigsq;
    cfvar = 1./cohen.^.5;

    palm_quicksave(cope,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cope',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(varcope,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_varcope',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(cohen,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cohenf2',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
    palm_quicksave(cfvar,0,opts,plm,y,m,c, ...
        sprintf('%s',opts.o,plm.Ykindstr{y},'_cfvar',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c},pstr));
end
