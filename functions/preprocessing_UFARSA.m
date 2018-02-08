function [fluors, opt] = preprocessing_UFARSA(opt,fluor)

% [fluors, opt] = preprocessing_UFARSA(opt,fluor)
%
% Pre-processing of the given fluorescence trace. The main steps considered
% in this function are: Normalization, Slow drifts removal, Deflections
% removal, Smoothing, and Noise level estimation.
%
% INPUT:
%   fluor: the given raw fluorescence trace
%   opt: an structure variable, where
%        
%   parameters of "Slow drifts removal"
%        opt.remove_drifts --> 1: remove slowly varying drifts, 0: skip the drift removal step
%        opt.detrending_method: select the drift removal method, 1: BEADS, 2:Polynomial
%        opt.fc_beads: in [cycles/sample], filter's normalized cut-off frequency in BEADS method (0<fc_beads<0.5)
%        opt.asymmetery_param: asymmetry ratio parameter in BEADS method
%        opt.beads_nAppendFrames: number of appended frames to the start and end of fluorescence trace, used for compensating for the potential boarder effects after applying BEADS method
%        opt.beads_fast --> 1: use our fast implementation of BEADS method, 0: use the original code of BEADS method
%        opt.filterOrder_param: filter's order in BEADS method (set it to 1 or 2)
%        opt.poly_degree: degree of polynomial function used by polynomial detrending method
%
%
%   parameters of "Deflections removal"
%        opt.remove_posDeflections --> 1: apply large-impulse (deflection) removal step, 0:skip this step
%        opt.posPrctileLevel:  positive percentile level, used in the large-impulse removal step
%        opt.scale_posDefl: a constant
%
%        opt.remove_negDeflections --> 1: remove large short-lasting negative deflections, 0:skip this step
%        opt.negPrctileLevel:  negative percentile level, used in the negative-deflections removal step
%        opt.scaleRep_negDefl: a constant
%
%
%   parameters of "smoothing & noise-level estimation"
%        opt.test_smoothing --> 1: check smoothing quality (recommended), 0: skip the check
%        opt.test_smallS --> 1: check smoothing quality when smoothing parameter is smaller than 1 (recommended for real data), 0: skip the check
%        opt.smoothing_param: smoothing parameter; []: automatic, [NUMBER]: the given NUMBER will be used for smoothing (non-automatic)
%        opt.Svalue_when_fail: the value of smoothing parameter to be used when smoothing was not proper
%        opt.samples_randomly --> 1: select the frames for computing the STD randomly, 0: select them in a non-random manner
%        opt.nSamplesForSTD --> []: consider the whole trace for estimating noise STD
%                            % [NUMBER]: consider the first, or randomly selected, NUMBER of frames to estimate noise STD
%                            % [NUMBER1 NUMBER2]: consider frames from frame NUMBER1 up to frame NUMBER2, to estimate noise STD
%        opt.warn_pms --> 1: show the potential warning messages in command window, 0: don't show them
%        opt.noiseLevel_method: noise-STD estimation method: 1 ---> Power-spectrum density (PSD), 0 ---> residual
%
%
% OUTPUT
%   fluors: an structure variable, where
%           fluors.original: raw (non-normalized) fluorescence trace
%           fluors.normalized: raw normalized fluorescence trace
%           fluors.beforeSmoothing: normalized pre-processed before-smoothing fluorescence trace or raw normalized fluorescence trace (to be used for smoothing)
%           fluors.afterSmoothing: fluorescence trace after smoothing (to be used for reconstruction)
%
% Author: Vahid Rahmati (December, 2017)


%% make the trace a row vector
opt.nFrames_original = numel(fluor);
fluor = fluor(:)';
fluors.original = fluor; % the raw trace before applying any pre-processing

opt.normalization       = 1; % 1: normalize the fluorescence trace (mandatory)
opt.Min_nFrames         = 1200; % something internal (skip this)

if sum(abs(fluor)) ~= 0
    
    %% Normalization: normalize the fluorescence trace
    if opt.normalization == 1
        baseline = mean(fluor);
        fluor = (fluor-baseline);
        fluor = fluor./max( baseline, 1); % normalized fluorescence trace
    end
    fluors.normalized = fluor;
    
    %% Slow drifts removal (temporal): removing low-frequency drifts (i.e. the baseline signal)
    if opt.remove_drifts
        switch opt.detrending_method
            case 1 % BEADS method
                detrend_fluor = detrending_beads_UFARSA(fluor', opt)';
                fluor = detrend_fluor; % detrended trace
            case 2 % Polynomial detrending method
                poly_degree = opt.poly_degree;
                [p,~,mu] = polyfit([1:numel(fluor)],fluor,poly_degree);
                f_y = polyval(p,[1:numel(fluor)],[],mu);
                detrend_fluor = fluor - f_y;
                fluor = detrend_fluor; % detrended trace
        end
    end
    
    
    %% Deflections removal: remove large positive impulses and/or large negative short-lasting defection artefacts
    if opt.remove_negDeflections || opt.remove_posDeflections
        [fluor_deflectFree, opt] = deflectionRemoval_UFARSA(opt, fluor);
        fluor = fluor_deflectFree;
    end
    
    
    %% Trace size (internal): check whether the number of Frames in the given trace is above the given minimum (just for robustness considerations)
    opt.flag_nFrames = [];
    Min_nFrames  = opt.Min_nFrames;
    if numel(fluor)<Min_nFrames
        fluor = repmat(fluor, 1, floor(Min_nFrames/numel(fluor)) + 1 );
        opt.flag_nFrames = 'Low';
    end
    nFrames = numel(fluor);
    opt.nFrames = nFrames;
    fluor_preSmth = fluor;
    fluors.beforeSmoothing = fluor_preSmth; % the detrended, normalized trace to be used for smoothing
    
    
    %% Smoothing
    [fluor_smth, Noise_trace, smoothing_param] = smoothing_UFARSA(opt,fluor_preSmth);
    opt.smoothing_param = smoothing_param;
    fluors.afterSmoothing = fluor_smth; % the smoothed trace
    
    
    %% Noise level estimation
    switch opt.noiseLevel_method
        case 1 % based on the PSD of the trace used for smoothing
            std_noise = GetSn(fluor_preSmth);
        case 0 % based on the smoothing residual
            std_noise = median(abs(Noise_trace-median(Noise_trace)))*1.4826;
    end
    
    opt.std_noise = std_noise;   
    opt.ZeroFluor_falg = 0;
    
else
    opt.ZeroFluor_falg = 1;
end

end