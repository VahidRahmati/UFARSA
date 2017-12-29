% The script "demo_UFASAR_realData.m" applies UFASAR to the "sample_data.mat" file containing a number of fluorescence
% traces (i.e. time-courses), which can be assumed to be extracted from different region of interests (ROIs). Simply run 
% this script and see UFASAR's spiking activity reconstruction results for these traces.
% This script can be adapted to apply UFASAR to new recorded fluorescence traces. The main blocks to be adapted by 
% user are those corresponding to setting the data directory, and the central reconstruction parameter; i.e. blocks 1
% and 3. The other lines in this script mainly provide the options for deciding on the steps of UFASAR. 
% For the rest of parameters as well as the saving and plotting options please open the "internal_parameters.m"
% file (in Matlab), and adjust it as you desired. Please also read the "user_guide.pdf" file for more details. 
% Hint: If this is the first time to execute this script, please be sure that you had already run the "steup_UFASAR.m" file. 
%
% Author: Vahid Rahmati (December, 2017)


clear all
clc; 
close all

%% 1: set the data directory    
opt.FluorFile_name = 'sample_data.mat'; % name of the file containing the fluorescence data; e.g. 'mydata.txt'
opt.FluorFile_dir  = [fileparts(mfilename('fullpath')),filesep]; % directory (folder) where your fluorescence data exist; e.g. 'D:\data\'


%% 2: select the ROIs 
opt.which_ROIs     = [];       % select the fluorescence traces of which ROIs in the file should be processed e.g. set [1 3 6] 
                               % [] --> all RIOs will be processed

                            
%% 3: set the central reconstruction parameter 
opt.scale_NoiseSTD = [2.25];   % Leading-threshold scaling constant (by default 2.25). 
%                              % We strongly recommended to estimate this parameter easily from a couple of you fluorescence traces (see user_guide.pdf
%                              % file). Following this estimation, we recommend to set the " opt.min_leading_amp = 0 " in the "internal_parameters.m" 
%                              % file, in order to remove the internally determined lower-bound used for this threshold.
        
                                                        
%% 4: decide on UFASAR's steps (for the rest of parameters see "internal_parameters.m" file)
opt.remove_drifts         = 0; % 1: remove slowly varying drifts, 0: skip the drift removal step
opt.remove_posDeflections = 0; % 1: apply large-impulse (deflection) removal step, 0:skip this step
opt.remove_negDeflections = 0; % 1: remove large short-lasting negative deflections, 0:skip this step
opt.demerging             = 1; % 1 --> apply the demerging step, 0 ---> skip it
opt.gen_FR_count          = 1; % 1 --> generate the estimated firing rate vector based on the reconstructed spike-count train, 0 --> skip it.
opt.gen_FR_count_dem      = 1; % 1 --> generate the estimated firing rate vector based on the reconstructed demerged spike-count train, 0 --> skip it.


%% 5: run UFASAR
[output_UFASAR,opt_out] = run_UFASAR(opt);


