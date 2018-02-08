% The script "demo_UFARSA_simData.m" simulates a Poisson spike train and the corresponding fluorescence trace (i.e.  
% time-course). Simply run this script and see UFARSA's spiking activity reconstruction results for this trace. You can  
% change the simulation parameters below in block 1. To apply UFARSA to your recorded data please see "demo_UFARSA_realData.m"
% file. Below, in blocks 2 and 3, you can set the central parameter of UFARSA and decide on the main steps of UFARSA. 
% For the rest of parameters as well as the saving and plotting options please open the "internal_parameters.m"
% file (in Matlab), and adjust it as you desired. Please also read the "user_guide.pdf" file for more details. 
% Hint: If this is the first time to execute this script, please be sure that you had already run the "steup_UFARSA.m" file.
%
% Author: Vahid Rahmati (December, 2017)


clear all
clc;
close all

%% 1: Simulation: simulate the fluorescence trace (i.e. fluorescence time-course)
fs            = 30;   % in [Hz], sampling frequency
firing_rate   = 0.4;  % in [Hz], expected mean firing rate of the simulated neuron
nFrames       = 3000; % number of simulated frames (for a single ROI)
tau           = 0.8;  % in [sec], decay time-constant of calcium transients
std_noise_sim = 0.25; % standard deviation of the additive white noise in fluorescence trace (note, the expected transient amplitude is 1) 
seed          = [15]; % random seed; using a different seed will lead to a different spike train and additive white noise 

% simulate the spike train and the corresponding fluorescence trace
[fluor, true_spikes] = sim_fluor_UFARSA(firing_rate, fs, tau, std_noise_sim, nFrames, seed);


%% 2: The central reconstruction parameter 
opt.scale_NoiseSTD   = [2.25];   % Leading-threshold scaling constant (by default 2.25). 
        

%% 3: Decision on UFARSA's steps
opt.remove_drifts         = 0; % 1: remove slowly varying drifts, 0: skip the drift removal step
opt.remove_posDeflections = 0; % 1: apply large-impulse (deflection) removal step, 0:skip this step
opt.remove_negDeflections = 0; % 1: remove large short-lasting negative deflections, 0:skip this step
opt.demerging             = 1; % 1 --> apply the demerging step, 0 ---> skip it
opt.gen_FR_count          = 1; % 1 --> generate the estimated firing rate vector based on the reconstructed spike-count train, 0 --> skip it.
opt.gen_FR_count_dem      = 1; % 1 --> generate the estimated firing rate vector based on the reconstructed demerged spike-count train, 0 --> skip it.


%% 4: run UFARSA
[output_UFARSA,opt_out] = run_UFARSA(opt,fluor,true_spikes);


