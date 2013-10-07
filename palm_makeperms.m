function [Pset,nP] = palm_makeperms(opts,plm)
% This is just a wrapper to 'plexi.m' and 'signflip.m', and it
% generates a sigle set of permutations based on the options
% and parameters given in 'opts' and 'plm' struct.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

if opts.EE,
    Pset = palm_plexi(plm.tmp.seq,opts.nP0,plm.tmp.EB,opts.PB,~opts.lx);
    if numel(Pset) == 1,
        opts.EE = false;
        opts.ISE = true;
    end
end
if opts.ISE,
    if opts.PB,
        Sset = palm_signflip(plm.tmp.EB,opts.nP0,~opts.lx);
    else
        Sset = palm_signflip(plm.tmp.N,opts.nP0,~opts.lx);
    end
end
if opts.EE && opts.ISE,
    Bset = palm_drawproducts(Pset,Sset,opts.nP0);
    Pset = Bset;
elseif ~opts.EE && opts.ISE,
    Pset = Sset;
end
nP = numel(Pset);
