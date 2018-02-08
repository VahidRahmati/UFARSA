function [detrend_sginal] = detrending_beads_UFARSA(fluor, opt)

% [detrend_sginal] = detrending_beads_UFARSA(fluor, opt)
% 
% this function prepares and pass the fluorescence trace through BEADS
% method, so that its slowly varying drifts are removed.
%
% INPUT
%   fluor: the given fluorescence trace
%   opt: an structure variable, where
%        opt.filterOrder_param: filter order parameter (d = 1 or 2)
%        opt.asymmetery_param: asymmetry ratio parameter
%        opt.fc_beads: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%        opt.beads_nAppendFrames: number of frams to be appended to the start and end of fluorescence trace,
%                                 used for compensating for the potential boarder effects after applying BEADS method
%
% OUTPUT
%   detrend_sginal: the dirft-corrected (i.e. detrended) fluorescence trace
%
% Author: Vahid Rahmati (December, 2017)

%% compensate for the potential boarder effects
m     = opt.beads_nAppendFrames;
fluor = [mean(fluor(1:10))*ones(m,1); fluor; mean(fluor(end-10:end))*ones(m,1)]; 

%% set the parameters for BEADS method
d = opt.filterOrder_param;          % d : filter order parameter (d = 1 or 2)
r = opt.asymmetery_param;           % r : asymmetry parameter

% Regularization parameters
lam0 = 1; lam1 = 1; lam2 = 1; 

%% run BEADS
switch opt.beads_fast
    case 1 % use our fast implementation
        [x1, f1] = beadsFast_UFARSA(fluor, opt.fc_beads, lam0);
    case 0 % use the original code
        [x1, f1, cost] = beads(fluor, d, opt.fc_beads, r, lam0, lam1, lam2);
end

detrend_sginal = fluor  - f1; % drift_corrected trace with appended frames
detrend_sginal = detrend_sginal(m+1:end-m); % drift_corrected trace without appended frames


