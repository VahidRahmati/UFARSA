function [fluor_deflectFree, opt] = deflectionRemoval_UAFASAR(opt, fluor)

% [fluor_deflectFree, opt] = deflectionRemoval_UAFASAR(opt, fluor)
%    
% Removing the large positive impulses and large short-lasting deflections 
% from fluorescence trace. Currently this function should only be applied to 
% drift-free traces or traces de-trended by e.g. BEADS method in UFASAR.
%
% INPUT:
%   fluors: the drift-free or drift-corrected fluorescence trace 
%
%   opt: an structure variable, where we need
%        opt.remove_posDeflections --> 1: apply large-impulse (deflection) removal step, 0:skip this step
%        opt.posPrctileLevel:  positive percentile level, used in the large-impulse removal step
%        opt.scale_posDefl: a constant
%
%        opt.remove_negDeflections --> 1: remove large short-lasting negative deflections, 0:skip this step
%        opt.negPrctileLevel:  negative percentile level, used in the negative-deflections removal step
%        opt.scaleRep_negDefl: a constant
%
% OUTPUT:
%
%   fluor_deflectFree: the drift-free or detrended trace, after removal of the positive and/or negative deflections
%   opt.defl_flag --> 1: at least one deflection was detected, and then corrected, 0: no deflection was found
%
% Author: Vahid Rahmati (December, 2017)

fluor_deflectFree = fluor;

%% remove negative deflection artefacts
if opt.remove_negDeflections
    prctile_neg = prctile(fluor_deflectFree(fluor_deflectFree<0),opt.negPrctileLevel); % this is equivalent to -prctile(fluor(fluor<0),opt.negPrctileLevel)
    cond_neg = fluor_deflectFree<prctile_neg;
    fluor_deflectFree(cond_neg) = opt.scaleRep_negDefl*prctile_neg;
    if sum(cond_neg) > 1; opt.defl_flag = 1; end
end


%% remove positive deflection artefacts (large impulses)
if opt.remove_posDeflections
    prctile_pos = prctile(fluor_deflectFree(fluor_deflectFree>0),opt.posPrctileLevel);
    cond_pos = fluor_deflectFree>opt.scale_posDefl*prctile_pos;
    
    if sum(cond_pos)>0 % true when there is at least one positive impulse
        
        % find the location of impulses
        impulses_idx = find([NaN diff(diff(cond_pos))] == -2);
        
        for iImpulse=1:numel(impulses_idx)
            idx_current = impulses_idx(iImpulse);
            fluor_deflectFree(idx_current) = 0.5*(fluor_deflectFree(idx_current-1) + fluor_deflectFree(idx_current+1));
        end
        if cond_pos(1) == 1 && impulses_idx(2) == 0; fluor_deflectFree(1) = fluor_deflectFree(2); end
        if cond_pos(end) == 1 && impulses_idx(end-1) == 0; fluor_deflectFree(end) = fluor_deflectFree(end-1); end
        
    end
    
    if sum(cond_pos) > 1; opt.defl_flag = 1; end
end

end