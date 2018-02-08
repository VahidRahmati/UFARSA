function [fluor, spikes] = sim_fluor_UFARSA(firing_rate,frame_rate,tau, std_noise_sim, nFrames, seed)

% [fluor, spikes] = sim_fluor_UFARSA(firing_rate,frame_rate,tau, std_noise_sim, nFrames, seed)
%
% This function simulates a Poisson spike train and its corresponding fluorescence trace.
%
% INPUT: 
%   firing_rate: in [Hz], expected mean firing rate of the simulated neuron
%   frame_rate: in [Hz], sampling frequency
%   tau: in [sec], decay time-constant of calcium transients
%   std_noise_sim: standard deviation of the additive white noise in fluorescence trace (note, the expected transient amplitude is 1) 
%   nFrames: number of simulated frames
%   seed: random seed; using a different seed will lead to a different spike train and additive white noise 
%
% Ouput:
%   fluor: simulated fluorescence trace, acquired at the user-determined sampling frequency
%   spikes: the underlying simulated Poisson spike-count at the given sampling frequency
%
% Author: Vahid Rahmati (December, 2017)


if ~isempty(seed); rng(seed); end

frame_size = 1/frame_rate;  % time step size
gam_tau   = 1-frame_size/tau;

spikes = poissrnd(firing_rate*frame_size*ones(nFrames,1)); % spike train simulation
calcium = filter(1,[1 -gam_tau],spikes); % [Ca2+] trace
fluor = calcium + std_noise_sim*randn(nFrames,1);   % fluorescence trace
fluor = fluor';
spikes = spikes';


