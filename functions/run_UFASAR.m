function [output_UFASAR, opt_out, ROIs_all] = run_UFASAR(opt, fluor, true_spikes)

% [output_UFASAR, opt_out, ROIs_all] = run_UFASAR(opt, fluor, true_spikes)
%
% Applying UFASAR on the given fluorescence trace(s). The trace(s) will be first pre-processed, and then spiking activity trains will be reconstructed,
% plotted and saved.
%
%
% INPUT:
%   fluor: the raw "simulated" fluorescence trace
%   true_spikes: the "simulated" train of spiking event times or spike-counts (the time-unit is [frame])
%
%   #Hint#: for real data don't assign any of "fluor" and "true_spikes" variables as input to this function; i.e. use instead: ... = run_UFASAR(opt)
%
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
%   parameters of "reconstruction"  
%        opt.demerging --> 1: apply the demerging step, 0: skip it
%        opt.gen_FR_count --> 1: generate the estimated firing rate vector based on the reconstructed spike-count train, 0: skip it.
%        opt.gen_FR_count_dem --> 1: generate the estimated firing rate vector based on the reconstructed demerged spike-count train, 0: skip it
%        opt.min_leading_amp --> []: leading threshold will be estimated, and be lower-bounded with the minimum leading threshold (default)
%                              % 0: minimum leading threshold will not be used to lower-bound the leading threshold (recommended when user estimated leading threshold from data)
%                              % [scalar]: an estimate of the minimum amplitude of non-within-burst transients (e.g. isolated or leading); see user_guide.pdf file, for more details.     
%
%        opt.scale_NoiseSTD: Leading-threshold scaling constant (by default 2.25). 
%                           % We strongly recommended to estimate this parameter easily from a couple of you fluorescence traces (see user_guide.pdf
%                           % file). Following this estimation, we recommended to set the " opt.min_leading_amp = 0 " in the "internal_parameters.m" 
%                           % file, in order to remove the internally determined lower-bound used for this threshold.
%
%        opt.denominator_prct: Min-leading-threshold scaling constant. 
%                             % To determine the lower bound for leading threshold, the 98th percentile of the smoothed fluorescence trace is divided by this number
%                             
%        opt.scale_burstAmp: Min-interior-amplitude scaling constant
%        opt.plateau_interior: Interior-plateau scaling constant
%        opt.scale_min_amp: Min-leading-amplitude scaling constant (see the user_guide.pdf file), this is used when " opt.min_leading_amp = [scalar] "  
%        opt.scale_min_interior_thr: Min-interior-threshold scaling constant (a number within the range [0 1]): decreasing this (e.g. to 0.5) can lead to a better reconstruction of within-burst spikes 
%        opt.max_nBurstSpikes --> [Number]: Maximum spike-count constant; maximum spike-count which can be reconstructed from a detected leading transient
%                               % [] --> no restriction (i.e. maximum bound) will be applied to the raw estimated spike-count (not recommended)
%
%        opt.onset_shift: Onset-shifting parameter (in [frame]); the corrected reconstructed event times will be shifted by "onset_shift" number of frames
%        opt.sigma_gauss: in [frame], the STD of Gaussian kernel used to be convolved with the spiking activity train, in order to generate estimated firing rate vector
%
%
%   parameters of "figures": plotting and saving
%        opt.plt  --> 1: plot the results, 0: skip plotting results
%        opt.save_plt --> 1: save the plotted figure, 0 --> don't save it
%        opt.save_plt_format: the format of saved figures; e.g. '.png' or '.jpg' or '.tif' or etc. 
%        opt.plt_eTrainOnFt --> 1: plot the user's selected type of the reconstructed event train on fluorescence trace(s), 0: don't show it
% 
%        opt.plt_eTrain --> 1: plot the reconstructed event train (E(t)), 0: don't show it
%        opt.plt_cTrain --> 1: plot the reconstructed spike-count train (C(t)), 0: don't show it
%        opt.plt_cFR --> 1: plot the estimated firing rate vector based on C(t), 0: don't show it
% 
%        opt.plt_eTrain_dem --> 1: plot the reconstructed demerged event train (E_dem(t)), 0: don't show it
%        opt.plt_cTrain_dem --> 1: plot the reconstructed demerged spike-count train (C_dem(t)), 0: don't show it
%        opt.plt_cFR_dem --> 1: plot the estimated firing rate vector based on C_dem(t), 0: don't show it
% 
%        opt.plt_eTrain_sim --> 1: plot the simulated event train, 0: don't show it
%        opt.plt_cTrain_sim --> 1: plot the simulated spike-count train, 0: don't show it
% 
%        opt.plt_Ft_orig_norm -->1: plot the normalized version of the raw fluorescence trace, 0: don't show it
%        opt.plt_Ft_beforeSmth -->1: plot the fluorescence trace before smoothing, 0: don't show it
%        opt.plt_Ft_afterSmth --> 1: plot the pre-processed fluorescence trace after smoothing, 0: don't show it
%
%        opt.display_NoiseSTD --> 1: display the estimated noise_std for each processed trace in the command window, 0: skip it
%
%   #Hint#: To save the figures, first a new folder called "UFASAR_results", and then new sub-folder called "Figures" will be created in the folder of data. 
%         If these folders already existed, then they won't be created and the already existed saved figures will be replaced by the new ones.
%
%
%   parameters of "saving UFASAR's results": reconstructions and fluorescence traces
%       opt.save_outputs --> 1: enable saving the reconstruction results and different versions of the fluorescence trace, 0: skip it
%       opt.save_output_format --> '.txt' or '.mat'; the files will be saved in this given format
% 
%       opt.save_eTrains --> 1: save the reconstructed event train (E(t)), 0: don't save it
%       opt.save_cTrains --> 1: save the reconstructed spike-count train (C(t)), 0: don't save it
%       opt.save_cFR --> 1: save the smoothed firing rate trace based on C(t), 0: don't save it
% 
%       opt.save_eTrains_dem --> 1: save the reconstructed demerged event train (E_dem(t)), 0: don't save it
%       opt.save_cTrains_dem --> 1: save the reconstructed demerged spike-count train (C_dem(t)), 0: don't save it
%       opt.save_cFR_dem --> 1: save the estimated firing rate vector based on C_dem(t), 0: don't save it
% 
%       opt.save_Ft_orig_norm --> 1: save the normalized version of the raw fluorescence trace, 0: don't save it
%       opt.save_Ft_beforeSmth --> 1: save the fluorescence trace before smoothing, 0: don't save it
%       opt.save_Ft_afterSmth --> 1:save the pre-processed fluorescence trace after smoothing, 0: don't save it
% 
%       opt.save_ROIs_all --> 1: save the cell variable containing the UFASAR's outputs for all fluorescence traces (ROIs), 0: don't save it
%
%   #Hint#: To save these variables, first a new folder called "UFASAR_results", and then a new sub-folder called "Reconstructions" will be created in the folder of data. 
%         If these folders already exist, then they won't be created and the already existed saved files will be replaced by the new ones.
%         
%   parameters for "setting the data directory"
%       opt.FluorFile_name --> name of the file containing your recorded fluorescence data; e.g. 'mydata.txt'
%       opt.FluorFile_dir --> directory where your recorded fluorescence data were saved; e.g. 'D:\data\'
%
%   #Hint#: when applying UFASAR to your recorded data, give only "opt" as an input to the run_UFASAR function; i.e. as in the demo_UFASAR_realData.m code use: ...= run_UFASAR(opt)
%
%
% OUTPUT:
%
%   output_UFASAR: an structure variable, where
%        output_UFASAR.eTrain: reconstructed "Event train" (in paper, denoted by E(t))
%        output_UFASAR.eTrain_dem: reconstructed "Demerged event train" (in paper, denoted by E_dem(t))
%        output_UFASAR.cTrain: reconstructed "Spike-count train" (in the paper, denoted by C(t))
%        output_UFASAR.cTrain_dem: reconstructed "Demerged spike-count train" (in the paper, denoted by C_dem(t))
%        output_UFASAR.cFR: estimated firing rate vector based on the reconstructed spike-count train, and the user-determined STD of Gaussian kernel
%        output_UFASAR.cFR_dem: estimated firing rate vector based on the reconstructed demerged spike-count train, and the user-determined STD of Gaussian kernel
%        output_UFASAR.leading_thr: Leading threshold
%        output_UFASAR.std_noise: estimated std of noise in the (normalized) fluorescence trace (it is computed in pre-processing step)
%        output_UFASAR.min_leading_thr_est: Minimum leading threshold
%        output_UFASAR.interior_thr: Interior threshold
%           
%        output_UFASAR.fluors.original: raw (non-normalized) fluorescence trace
%        output_UFASAR.fluors.normalized: raw normalized fluorescence trace
%        output_UFASAR.fluors.beforeSmoothing: normalized pre-processed before-smoothing fluorescence trace or raw normalized fluorescence trace (to be used for smoothing)
%        output_UFASAR.fluors.afterSmoothing: fluorescence trace after smoothing (to be used for reconstruction)       
%
%   opt_out: an structure variable, including all parameters assigned to the "opt" structure variable by user, preprocessing_UFASAR function, and within this current function.
%
%   ROIs_all: a cell variable including the UFASAR's outputs (i.e. both output_UFASAR and opt_out) for each of the fluorescence traces; 
%             e.g. to access the outputs for ROI number 3, write ROIs_all{3} in the Matlab Command Window prompt.
%             To remove this output (recommended) use "[output_UFASAR,opt_out] = run_UFASAR(...", instead of "[output_UFASAR,opt_out, ROIs_all] = run_UFASAR(..." 
%   
%   #Hint#: Note that allowing "ROIs_all" to be an output of UFASAR may significantly reduce the speed of the computer, especially in the case of long-term data. This is because 
%           all different versions of all fluorescence traces (i.e. of all ROIs) as well as all reconstruction variables will be stored in the working memory. Instead, we suggest 
%           using the abovementioned saving options to save the outputs of interest. 
%   #Hint#: If there are more than one fluorescence traces need to be processed, then the output variables "output_UFASAR" and "opt_out" will correspond to the last trace.
%
% Author: Vahid Rahmati (December, 2017)

%% import/load data
sim_flag = exist('true_spikes') && exist('fluor');
opt.sim_flag = sim_flag;

ROIs_all_flag = 0;
if nargout>2; ROIs_all_flag = 1; end

if ~sim_flag
    FluorFile_format = lower(opt.FluorFile_name(end-2:end)); % only .txt and .mat files are accepted (see README file)
    
    switch FluorFile_format
        case 'mat'
            fluor_data = importdata([opt.FluorFile_dir,opt.FluorFile_name]);
            nROIs = size(fluor_data,2);
        case 'txt'
            All_data = importdata([opt.FluorFile_dir,opt.FluorFile_name]);
            fluor_data = All_data.data(:,2:end);
            nROIs = size(fluor_data,2) - 1;
    end
    
    % select fluorescence trace belonging to which RIOs should be processed
    if ~isempty(opt.which_ROIs)
        ROIs_vec = opt.which_ROIs;
    else
        ROIs_vec = 1:nROIs;
    end
    
else
    ROIs_vec = 1;
    fluor_data = fluor(:);
end


%% read the rest of UFASAR's parameter values
internal_parameters % If desired, user can open the "internal_parameters.m" file, and adjust the values of the rest of UFASAR's parameters
s_orig = opt.smoothing_param;


%% run UFSAR
opt.warn_flag = 1;
for iROI = ROIs_vec
    opt.iROI = iROI;
    fluor = fluor_data(:,iROI);
    
    %% Pre-processing
    [fluors, opt] = preprocessing_UAFASAR(opt,fluor);
    
    %% Reconstruction
    [output_UFASAR,opt_out] = reoncstruction_UFASAR(fluors,opt);
    
    %% store all UFASAR's outputs for all ROIs 
    if ROIs_all_flag || opt.save_ROIs_all
        ROIs_all{iROI}.output_UFASAR = output_UFASAR;
        ROIs_all{iROI}.opt_out = opt_out;
    end
    
    %% Plotting && saving
    opt.ZeroEvent_flag = opt_out.ZeroEvent_flag;
    if opt.plt
        
        switch sim_flag
            case 1
                plot_UFASAR(output_UFASAR,opt,true_spikes);
            case 0
                h_fig = plot_UFASAR(output_UFASAR,opt);
                
                % save figure
                if opt.save_plt; save_figure(h_fig, opt, iROI); end
                
                % save reconstructions and fluorescence traces
                if opt.save_outputs; save_recVectors(output_UFASAR,opt,iROI), end
                
        end  
        opt.warn_flag = 0;
    end
    
    opt.smoothing_param = s_orig; % useful when more than one ROIs need to be processed
    
    % show the estimated noise level
    if opt.display_NoiseSTD; display(sprintf('estimated normalized std_noise for ROI(%d) =  %d', iROI, output_UFASAR.std_noise)), end
    
end

% save the UFASAR's outputs for each fluorescence trace
if opt.save_ROIs_all && sim_flag==0

    if exist([opt.FluorFile_dir,'UFASAR_results'])==0
        mkdir([opt.FluorFile_dir,'UFASAR_results'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\ROIs_all','.mat']) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\ROIs_all','.mat'])
    end
    save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\ROIs_all','.mat'],'ROIs_all')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% -----------------------> Local Functions <-----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function save_figure(h_fig, opt, iROI)
    % save the plotted figures
    
set(h_fig, 'PaperPositionMode', 'auto')

if exist([opt.FluorFile_dir,'UFASAR_results'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results'])
end

if exist([opt.FluorFile_dir,'UFASAR_results\Figures'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Figures'])
end

saveas(h_fig, [opt.FluorFile_dir,'UFASAR_results\Figures\','ROI(',num2str(iROI),')',opt.save_plt_format])

end

function save_recVectors(output_UFASAR, opt,iROI)
    % save the reconstruction results and pre-processed fluorescence traces

if exist([opt.FluorFile_dir,'UFASAR_results'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results'])
end

if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions'])
end


%% save the reconstructed event train
if opt.save_eTrains
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\'])
    end
    
    eTrain_times = [[1:opt.nFrames_original]' output_UFASAR.eTrain'];
    event_times_fire = find(output_UFASAR.eTrain>0)';
    event_times_fire = [event_times_fire];
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')',opt.save_output_format ])
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ])
    end
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')',opt.save_output_format ],eTrain_times,'delimiter','\t','newline','pc')
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ],event_times_fire,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')',opt.save_output_format ],'eTrain_times')
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains\eTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ],'event_times_fire')        
    end
    
end

%% save the reconstructed, demerged event train
if opt.save_eTrains_dem && opt.demerging
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\'])
    end
    
    eTrain_times_dem = [[1:opt.nFrames_original]' output_UFASAR.eTrain_dem'];
    event_times_dem_fire = find(output_UFASAR.eTrain_dem>0)';
    event_times_dem_fire = [event_times_dem_fire];
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ])
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ])
    end
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ],eTrain_times_dem,'delimiter','\t','newline','pc')
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ],event_times_dem_fire,'delimiter','\t','newline','pc')
        case '.mat'            
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ],'eTrain_times_dem')
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\eTrains_dem\eTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ],'event_times_dem_fire')          
    end
end


%% save the reconstructed spike-count train
if opt.save_cTrains

    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\'])
    end
    
    cTrain_times = [[1:opt.nFrames_original]' output_UFASAR.cTrain'];
    cTrain_times_fire = find(output_UFASAR.cTrain>0);
    cTrain_times_fire = [cTrain_times_fire' output_UFASAR.cTrain(cTrain_times_fire)'];
     
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')',opt.save_output_format ])
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ])
    end
    
    switch opt.save_output_format
        case '.txt'            
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')',opt.save_output_format ],cTrain_times,'delimiter','\t','newline','pc')
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ],cTrain_times_fire,'delimiter','\t','newline','pc')            
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')',opt.save_output_format ],'cTrain_times')
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains\cTrain_ROI(',num2str(iROI),')_fire',opt.save_output_format ],'cTrain_times_fire')
    end
end


%% save the reconstructed, demerged spike-count train
if opt.save_cTrains_dem && opt.demerging

    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\'])
    end
    
    cTrain_times_dem = [[1:opt.nFrames_original]' output_UFASAR.cTrain_dem'];
    cTrain_times_dem_fire = find(output_UFASAR.cTrain_dem>0);
    cTrain_times_dem_fire = [cTrain_times_dem_fire' output_UFASAR.cTrain_dem(cTrain_times_dem_fire)'];
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ])
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ])
    end
    
    switch opt.save_output_format
        case '.txt'            
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ],cTrain_times_dem,'delimiter','\t','newline','pc')
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ],cTrain_times_dem_fire,'delimiter','\t','newline','pc')            
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')',opt.save_output_format ],'cTrain_times_dem')
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cTrains_dem\cTrain_dem_ROI(',num2str(iROI),')_fire',opt.save_output_format ],'cTrain_times_dem_fire')
    end
end


%% save the estimated firing rate vectors based on the reconstructed spike-count train
if opt.save_cFR && opt.gen_FR_count
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\cFR_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\cFR_ROI(',num2str(iROI),')',opt.save_output_format ])
    end
    
    cFR_vector = [[1:opt.nFrames_original]' output_UFASAR.cFR'];
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\cFR_ROI(',num2str(iROI),')',opt.save_output_format],cFR_vector,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR\cFR_ROI(',num2str(iROI),')',opt.save_output_format],'cFR_vector')           
    end
end


%% save the estimated firing rate vector based on the reconstructed, demerged spike-count train
if opt.save_cFR_dem && opt.gen_FR_count_dem && opt.demerging
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\cFR_dem_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\cFR_dem_ROI(',num2str(iROI),')',opt.save_output_format ])
    end

    cFR_vector_dem = [[1:opt.nFrames_original]' output_UFASAR.cFR_dem'];
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\cFR_dem_ROI(',num2str(iROI),')',opt.save_output_format],cFR_vector_dem,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\cFR_dem\cFR_dem_ROI(',num2str(iROI),')',opt.save_output_format],'cFR_vector_dem')
    end
end


%% save the normalized version of the given raw fluorescence trace
if opt.save_Ft_orig_norm
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\Ft_normalized_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\Ft_normalized_ROI(',num2str(iROI),')',opt.save_output_format ])
    end
    
    Ft_orig_norm = [[1:opt.nFrames_original]' output_UFASAR.fluors.normalized'];
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\Ft_normalized_ROI(',num2str(iROI),')',opt.save_output_format],Ft_orig_norm,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_normalized\Ft_normalized_ROI(',num2str(iROI),')',opt.save_output_format],'Ft_orig_norm')
    end
end

%% save the (pre-processed) trace used in smoothing
if opt.save_Ft_beforeSmth
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\Ft_beforeSmoothing_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\Ft_beforeSmoothing_ROI(',num2str(iROI),')',opt.save_output_format ])
    end
    
    Ft_pre_smth = [[1:opt.nFrames_original]' output_UFASAR.fluors.beforeSmoothing'];
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\Ft_beforeSmoothing_ROI(',num2str(iROI),')',opt.save_output_format],Ft_pre_smth,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_beforeSmoothing\Ft_beforeSmoothing_ROI(',num2str(iROI),')',opt.save_output_format],'Ft_pre_smth')
    end
end



%% save the smoothed trace used in reconstruction
if opt.save_Ft_afterSmth 
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\'])==0
    mkdir([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\'])
    end
    
    if exist([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\Ft_afterSmoothing_ROI(',num2str(iROI),')',opt.save_output_format ]) ~= 0
        delete([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\Ft_afterSmoothing_ROI(',num2str(iROI),')',opt.save_output_format ])
    end
    
    Ft_after_smth= [[1:opt.nFrames_original]' output_UFASAR.fluors.afterSmoothing'];
    
    switch opt.save_output_format
        case '.txt'
            dlmwrite([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\Ft_afterSmoothing_ROI(',num2str(iROI),')',opt.save_output_format],Ft_after_smth,'delimiter','\t','newline','pc')
        case '.mat'
            save([opt.FluorFile_dir,'UFASAR_results\Reconstructions\Ft_afterSmoothing\Ft_afterSmoothing_ROI(',num2str(iROI),')',opt.save_output_format],'Ft_after_smth')          
    end
end


end
