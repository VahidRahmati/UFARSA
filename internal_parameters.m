% Internal_parameters
% 
% This script enables the user to access the internal parameters and options
% of UFARSA. After changing the value of any parameter or option, just SAVE 
% the script (CTRL+S), and run again the main code (e.g. demo_UFARSA_realData.m).
% The default parameter values used in this script are according to the paper in which
% we have presented the UFARSA. Please read the "user_guide.pdf" file for more details.
%
% Author: Vahid Rahmati (December, 2017)


%% options for plotting and saving figures: applicable when you've set "opt.plt = 1"
opt.plt                 = 1;       % 1 --> plot the results, 0 ---> skip plotting results
opt.save_plt            = 0;       % 1 --> save the plotted figure, 0 --> don't save it
opt.save_plt_format     = '.png';  % the format of saved figures; e.g. '.png' or '.jpg' or '.tif' or etc. 
opt.plt_eTrainOnFt      = 1;       % 1 --> plot the user's selected type of reconstructed event train on the fluorescence trace(s), 0: don't show it

opt.plt_eTrain          = 1;       % 1 --> plot the reconstructed event train (E(t)), 0 --> don't show it
opt.plt_cTrain          = 1;       % 1 --> plot the reconstructed spike-count train (C(t)), 0 --> don't show it
opt.plt_cFR             = 1;       % 1 --> plot the estimated firing rate vector based on C(t), 0 --> don't show it

opt.plt_eTrain_dem      = 1;       % 1 --> plot the reconstructed demerged event train (E_dem(t)), 0 --> don't show it
opt.plt_cTrain_dem      = 1;       % 1 --> plot the reconstructed demerged spike-count train (C_dem(t)), 0 --> don't show it
opt.plt_cFR_dem         = 1;       % 1 --> plot the estimated firing rate vector based on C_dem(t), 0 --> don't show it

opt.plt_eTrain_sim      = 1;       % 1 --> plot the simulated event train, 0 --> don't show it
opt.plt_cTrain_sim      = 1;       % 1 --> plot the simulated spike-count train, 0 --> don't show it

opt.plt_Ft_orig_norm    = 1;       % 1 --> plot the normalized version of the raw fluorescence trace, 0 --> don't show it
opt.plt_Ft_beforeSmth   = 1;       % 1 --> plot the fluorescence trace before smoothing, 0 --> don't show it
opt.plt_Ft_afterSmth    = 1;       % 1 --> plot the pre-processed fluorescence trace after smoothing, 0 --> don't show it

opt.display_NoiseSTD    = 1;       % 1 --> 1: display the estimated noise_std for each processed trace in the command window, 0 --> skip it


%% options for saving the reconstruction results and pre-processed fluorescence traces: applicable when you've set "opt.save_outputs = 1"
opt.save_outputs        = 0;       % 1 --> enable saving the reconstruction results and different versions of the fluorescence trace, 0 --> skip it
opt.save_output_format  = '.txt';  % '.txt' or '.mat'; the files will be saved in this given format

opt.save_eTrains        = 1;       % 1 --> save the reconstructed event train (E(t)), 0 --> don't save it
opt.save_cTrains        = 1;       % 1 --> save the reconstructed spike-count train (C(t)), 0 --> don't save it
opt.save_cFR            = 1;       % 1 --> save the estimated firing rate vector based on C(t), 0 --> don't save it

opt.save_eTrains_dem    = 1;       % 1 --> save the reconstructed demerged event train (E_dem(t)), 0 --> don't save it
opt.save_cTrains_dem    = 1;       % 1 --> save the reconstructed demerged spike-count train (C_dem(t)), 0 --> don't save it
opt.save_cFR_dem        = 1;       % 1 --> save the estimated firing rate vector based on C_dem(t), 0 --> don't save it

opt.save_Ft_orig_norm   = 1;       % 1 --> save the normalized version of the raw fluorescence trace, 0 --> don't save it
opt.save_Ft_beforeSmth  = 1;       % 1 --> save the fluorescence trace before smoothing, 0 --> don't save it
opt.save_Ft_afterSmth   = 1;       % 1 --> save the pre-processed fluorescence trace after smoothing, 0 --> don't save it

opt.save_ROIs_all       = 0;       % 1 --> save the cell variable containing the UFARSA's outputs for all fluorescence traces (ROIs), 0: don't save it


%% parameters of pre-processing (detrending): applicable when you've set "opt.remove_drifts = 1"
opt.detrending_method   = 1;       % select the drift removal method, 1: BEADS, 2:Polynomial 
opt.fc_beads            = 0.002;   % in [cycles/sample], filter's normalized cut-off frequency in BEADS method (0<fc_beads<0.5)
opt.asymmetery_param    = 1;       % asymmetry ratio parameter in BEADS mehtod
opt.beads_nAppendFrames = 10000;   % number of appended frames to the start and end of fluorescence trace, used for compensating for the potential boarder effects after applying BEADS method
opt.beads_fast          = 1;       % 1: use our fast implementation of BEADS method, 0: use the original code of BEADS method
opt.filterOrder_param   = 1;       % filter's order in BEADS method (set it to 1 or 2)
opt.poly_degree         = 3;       % degree of polynomial function used by polynomial detrending method


%% parameters of pre-processing (deflection removal)
opt.posPrctileLevel     = 98;      % positive percentile level, used in the larg-impulse removal step
opt.scale_posDefl       = 3;

opt.negPrctileLevel     = 2;       % negative percentile level, used in the negative-deflections removal step
opt.scaleRep_negDefl    = 0.5; 


%% parameters of pre-processing (smoothing + noise-level estimation)                      
opt.smoothing_param     = [];      % smoothing parameter; []: automatic, [NUMBER]: the given NUMBER will be used for smoothing (non-automatic) 
opt.test_smoothing      = 1;       % 1: check smoothing quality (recommended), 0: skip the check  
opt.test_smallS         = 1;       % 1: check smoothing quality when smoothing parameter is smaller than 1 (recommended for real data), 0: skip the check
opt.Svalue_when_fail    = 350;     % the value of smoothing parameter to be used when smoothing was not proper
opt.samples_randomly    = 0;       % 1: select the frames for computing the STD randomly, 0: select them in a non-random manner
opt.nSamplesForSTD      = 1000;    % [] --> consider the whole trace for estimating noise STD
%                                  % [NUMBER] --> consider the first, or randomly selected, NUMBER of frames to estimate noise STD     
%                                  % [NUMBER1 NUMBER2] --> consider frames from frame NUMBER1 up to frame NUMBER2, to estimate noise STD  
opt.noiseLevel_method   = 0;       % noise-STD estimation method: 1 --> Power-spectrum density (PSD), 0 ---> residual
opt.warn_pms            = 0;       % 1: show the potential warning messages in command window, 0: don't show them                            


%% parameters of reconstruction
opt.min_leading_amp     = [];      % [] --> leading threshold will be estimated, and be lower-bounded with the minimum leading threshold (default)
%                                  % 0  --> minimum leading threshold will not be used to lower-bound the leading threshold (recommended when user estimated leading threshold from data)
%                                  % [scalar] --> an estimate of the minimum amplitude of non-within-burst transients (e.g. isolated or leading); see user_guide.pdf file, for more details.     

opt.denominator_prct    = 6;       % Min-leading-threshold scaling constant. 
%                                  % To determine the lower bound for leading threshold, the 98th percentile of the smoothed fluorescence trace is divided by this number
                             
opt.scale_burstAmp      = 0.2;     % Min-interior-amplitude scaling constant
opt.plateau_interior    = 0.25;    % Interior-plateau scaling constant
opt.scale_min_amp       = 0.75;    % Min-leading-amplitude scaling constant (see the user_guide.pdf file), this is used when " opt.min_leading_amp = [scalar] "  
opt.scale_min_interior_thr = 0.75; % Min-interior-threshold scaling constant (a number within the range [0 1]): decreasing this (e.g. to 0.5) can lead to a better reconstruction of within-burst spikes 
opt.max_nBurstSpikes    = 5;       % Maximum spike-count constant; maximum spike-count which can be reconstructed from a detected leading transient
                                   % [] --> no restriction (i.e. maximum bound) will be applied to the raw estimated spike-count (not recommended)
opt.round_border        = 0.75;    % Border used for the imbalanced rounding of the estimated spike-counts: a scaler between 0.5-1 (default: 0.75).
opt.onset_shift         = 1;       % Onset-shifting parameter (in [frame]); the corrected reconstructed event times will be shifted for "onset_shift" number of frames
opt.sigma_gauss         = 4;       % in [frame], the STD of Gaussian kernel used to be convolved with the spiking activity train, in order to generate estimated firing rate vector
