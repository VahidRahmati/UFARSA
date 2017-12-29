function [fluor_smoothed, Noise_trace, smoothing_param] = smoothing_UFASAR(opt,fluor_preSmoothing)

% [fluor_smoothed, Noise_trace, smoothing_param] = smoothing_UFASAR(opt,fluor_preSmoothing)
%
% this function performs a high-frequency noise reduction of the fluroescence trace, based on
% the SMOOTHN method. Following this, the function will also compute the noise trace (i.e. the residual). 
%
% INPUT
%   fluor_preSmoothing: the fluorescence trace
%   opt: an structure variable, where
%        opt.test_smoothing --> 1: check smoothing quality (recommended), 0: skip the check
%        opt.test_smallS --> 1: check smoothing quality when smoothing parameter is smaller than 1 (recommended for real data), 0: skip the check
%        opt.smoothing_param: smoothing parameter; []: automatic, [NUMBER]: the given NUMBER will be used for smoothing (non-automatic)
%        opt.Svalue_when_fail: the value of smoothing parameter to be used when smoothing was not proper
%        opt.samples_randomly --> 1: select the frames for computing the STD randomly, 0: select them in a non-random manner
%        opt.nSamplesForSTD --> []: consider the whole trace for estimating noise STD
%                            % [NUMBER]: consider the first, or randomly selected, NUMBER of frames to estimate noise STD
%                            % [NUMBER1 NUMBER2]: consider frames from frame NUMBER1 up to frame NUMBER2, to estimate noise STD
%        opt.warn_pms --> 1: show the potential warning messages in command window, 0: don't show them
%        opt.nFrames: number of frames of the fluorescence trace
%
% OUTPUT
%   fluor_smoothed: the smoothed fluorescence trace (this trace will be passed through the reconstruction process)
%   Noise_trace: the computed residual trace, which will be used to estimate the STD of noise in the fluor_preSmoothing
%   smoothing_param: the smoothing parameter (S) used by the SMOOTHN method to filter out the high frequency noise 
%   
% Author: Vahid Rahmati (December, 2017)

    
%% smoothing
if opt.warn_pms == 0; warning off; end
lastwarn('')
[fluor_smoothed, smoothing_param] = smoothn(fluor_preSmoothing,opt.smoothing_param);


%% compute the noise trace (i.e. the residual)
nSamplesForSTD = opt.nSamplesForSTD; 
empty_nSamples = isempty(nSamplesForSTD);
switch empty_nSamples
    case 0
        try
            if nSamplesForSTD(2)> opt.nFrames; nSamplesForSTD(2) = opt.nFrames; end
        catch
            if nSamplesForSTD> opt.nFrames; nSamplesForSTD = opt.nFrames; end
        end
        
    otherwise
        nSamplesForSTD = opt.nFrames;
end

% select samples for estimating the noise STD
if opt.samples_randomly == 1 && size(nSamplesForSTD,2)==1      % select "nSamplesForSTD" number of frames randomly
    Samples_indices = randperm(opt.nFrames, nSamplesForSTD); 
elseif opt.samples_randomly == 0 && size(nSamplesForSTD,2)==1  % select the first "nSamplesForSTD" number of frames
    Samples_indices = 1:nSamplesForSTD;
elseif size(nSamplesForSTD,2)==2  % select the frames within the interval determined by "nSamplesForSTD"
    Samples_indices = nSamplesForSTD(1):nSamplesForSTD(2);
end
Noise_trace = fluor_preSmoothing(Samples_indices) - fluor_smoothed(Samples_indices); % Noise trace (residual)


%% check whether smoothing was performed reliably
if opt.test_smoothing ==1 && isempty(opt.smoothing_param) 
    
    [~, msgid] = lastwarn;
    
    cond1 = strcmp('MATLAB:smoothn:SLowerBound',msgid); % sub-smoothing condition
    cond2 = strcmp('MATLAB:smoothn:SUpperBound',msgid); % over-smoothing condition
    cond3 = strcmp('MATLAB:smoothn:MaxIter',msgid); % non-converged condition
    
    if smoothing_param<=1 && opt.test_smallS==1 && opt.sim_flag==0, cond4 = 1; else cond4 = 0; end % when computed smoothing parameter is less than 1
    
    if cond1 || cond2 || cond3 || cond4
        
        smoothing_param = opt.Svalue_when_fail;
        warning('failed_smoothing condition was used')

        % re-run the smoothing
        OPTIONS.TolZ = 0.1;
        OPTIONS.MaxIter = 100;
        [fluor_smoothed, smoothing_param] = smoothn(fluor_preSmoothing,smoothing_param,OPTIONS);
        
        % re-compute the residual trace (i.e. the noise)
        Noise_trace = fluor_preSmoothing(Samples_indices) - fluor_smoothed(Samples_indices);
        
    end
    
end

if opt.warn_pms == 0; warning on; end

end