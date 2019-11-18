function [output_UFARSA, opt_out] = reconstruction_UFARSA(fluors,opt)

% [output_UFARSA opt_out] = reconstruction_UFARSA(fluors,opt)
%
% Reconstruction of spiking activities from the smoothed fluorescence trace.
% The main steps included in this function are used for: Setting the leading threshold, Reconstruction,
% timing correction, Spike-count estimation, Demerging, and Building spiking activity trains.
%
% INPUT:
%   fluors: an structure variable, where we need
%           fluors.beforeSmoothing: drift-free or given normalized fluorescenec trace to which smoothing step was applied
%           fluors.afterSmoothing: fluorescence trace after smoothing (to be used for reconstruction)
%           
%   opt: an structure variable, where we need
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
%        opt.scale_burstAmp: Min-interior-amplitude scaling constant;
%        opt.plateau_interior: Interior-plateau scaling constant
%        opt.scale_min_amp: Min-leading-amplitude scaling constant (see the user_guide.pdf file), this is used when " opt.min_leading_amp = [scalar] "  
%        opt.scale_min_interior_thr: Min-interior-threshold scaling constant (a number within the range [0 1]): decreasing this can lead to a better reconstruction of within-burst spikes 
%        opt.max_nBurstSpikes --> [Number]: Maximum spike-count constant; maximum spike-count which can be reconstructed from a detected leading transient
%                               % [] --> no restriction (i.e. maximum bound) will be applied to the raw estimated spike-count (not recommended)
%
%        opt.round_border: Border used for the imbalanced rounding of the estimated spike-counts: a scaler between 0.5-1 (default: 0.75)
%        opt.onset_shift: Onset-shifting parameter (in [frame]); the corrected reconstructed event times will be shifted by "onset_shift" number of frames
%        opt.sigma_gauss: in [frame], the STD of Gaussian kernel used to be convolved with the spiking activity train, in order to generate estimated firing rate vector
%        opt.ZeroFluor_falg --> 1: the fluorescence trace is a zero vector, 0: the tracce is non-zero   
%        opt.nFrames_original: number of frames in the original given fluorescence trace
%
%
% OUTPUT: 
%   output_UFARSA: an structure variable, where
%        output_UFARSA.eTrain: reconstructed "Event train" (in paper, denoted by E(t))
%        output_UFARSA.eTrain_dem: reconstructed "Demerged event train" (in paper, denoted by E_dem(t))
%        output_UFARSA.cTrain: reconstructed "Spike-count train" (in the paper, denoted by C(t))
%        output_UFARSA.cTrain_dem: reconstructed "Demerged spike-count train" (in the paper, denoted by C_dem(t))
%        output_UFARSA.cFR: estimated firing rate vector based on the reconstructed spike-count train, and the user-determined STD of Gaussian kernel
%        output_UFARSA.cFR_dem: estimated firing rate vector based on the reconstructed demerged spike-count train, and the user-determined STD of Gaussian kernel
%        output_UFARSA.leading_thr: Leading threshold
%        output_UFARSA.std_noise: estimated std of noise in the (normalized) fluorescence trace (it is computed in pre-processing step)
%        output_UFARSA.min_leading_thr_est: Minimum leading threshold
%        output_UFARSA.interior_thr: Interior threshold
%           
%        output_UFARSA.fluors.original: raw (non-normalized) fluorescence trace
%        output_UFARSA.fluors.normalized: raw normalized fluorescence trace
%        output_UFARSA.fluors.beforeSmoothing: normalized pre-processed before-smoothing fluorescence trace or raw normalized fluorescence trace (to be used for smoothing)
%        output_UFARSA.fluors.afterSmoothing: fluorescence trace after smoothing (to be used for reconstruction)       
%
%   opt_out: an structure variable, including all parameters assigned to the "opt" structure variable by user, preprocessing_UFARSA function, and within this current function.
%
% Athor: Vahid Rahmati (December, 2017)

if opt.ZeroFluor_falg == 0
    
    %% Setting the leading threshold: determine the threshold for detecting leading events (thr_leading)
    [thr_leading,opt] = set_leading_thr(opt,fluors);
    
    
    %% reconstruct the onsets of all leading and interior events (includes the "Setting the interior threshold" step)
    [idx_onset_actu_All,idx_rise_mat_All,idx_rise_vec_All,Amps_actu_All,L_indices_aux,Amps_actu_L,ZeroEvent_flag, thr_wburst] = rec_events(fluors,thr_leading, opt);
    
else
    ZeroEvent_flag = 1; % 1: no event was detected for the given trace
end


if ZeroEvent_flag ~=1 % when there are non-zero amount of detected events
    
    %% Timing correction: refine the reconstructed times of both leading and interior events (except the demerged ones)
    [refined_OnsetTime,Frise_smth_orig_L,A_aux_L,refined_OnsetTime_L, refined_OnsetTime_B] = timing_correction(idx_onset_actu_All,idx_rise_mat_All,idx_rise_vec_All, fluors, Amps_actu_All, L_indices_aux );
    EventTimes_noSR = refined_OnsetTime; % the refined onset times of all events (i.e. leading or interior), except the merged ones
   
    
    %% Spike-count estimation: compute the relative spike-counts per detected leading fluorescence transient
    [burst_leaders_discrete_indices, vector_nSpikes_discrete] = spikecount_estimation(opt, Amps_actu_L, thr_leading);
    
    
    %% Demerging: compute the super-resolution event and spike-count trains, by dissecting the merged transients
    de_merging % an internal function
    
    
    %% Reconstruction of spiking activity trains: build the final output-vectors of the reconstruction process
    output_UFARSA = build_trains(EventTimes_SR,vector_nSpikes_discrete_SR,idx_onset_Actu_B_SR,ZeroEvent_flag);
    
    
else
    output_UFARSA = build_trains([],[],[],ZeroEvent_flag);
end


% write the main parameters of UFARSA into its output variables
output_UFARSA.leading_thr = opt.thr_leading;
output_UFARSA.std_noise = opt.std_noise;
output_UFARSA.min_leading_thr_est = opt.min_thr_leading_est;
output_UFARSA.interior_thr = thr_wburst;
opt_out = opt;
opt_out.ZeroEvent_flag = ZeroEvent_flag;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% -----------------------> Local Functions <-----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%
    function de_merging
        % Demerging: compute the super-resolution event- and count-trains, by dissecting the merged transients
        
        vector_nSpikes_discrete_SR = vector_nSpikes_discrete;
        if opt.demerging == 1 % when the user has selected to apply the Demerging step

            if ~isempty(burst_leaders_discrete_indices)
                
                refined_OnsetTime_aux = refined_OnsetTime_L(burst_leaders_discrete_indices,:);
                Frise_smth_orig_aux = Frise_smth_orig_L(burst_leaders_discrete_indices,:);
                A_aux2 = A_aux_L(burst_leaders_discrete_indices,:);
                
                sub_indices3 = coloncatrld(A_aux2(:,1),refined_OnsetTime_aux-1);
                nan_loc4 = ismember(A_aux2,sub_indices3);
                Frise_smth_orig_aux(nan_loc4) = NaN;
                A_aux2(nan_loc4) = NaN;
                Rising_diffdiff = diff(Frise_smth_orig_aux,2,2); % second derivative of the refined rising phase
                
                condneg_mat = zeros(size(Rising_diffdiff));
                condneg_mat(Rising_diffdiff<=0) = 1;
                condpos_mat = zeros(size(Rising_diffdiff));
                condpos_mat(Rising_diffdiff>0) = 1;
                
                if sum(sum(condneg_mat))~=0 && sum(sum(condpos_mat))~=0
                    
                    [i_neg j_neg] = find(Rising_diffdiff<=0);
                    [i_neg idx_sort_neg] = sort(i_neg);
                    j_neg = j_neg(idx_sort_neg);
                    
                    [~, i_neg_first] = unique(i_neg(end:-1:1));
                    i_neg_first = numel(i_neg) - i_neg_first + 1;
                    i_neg = i_neg(i_neg_first);
                    j_neg_first = j_neg(i_neg_first);
                    A_aux3 = A_aux2(:,1:size(A_aux2,2)-2); % becasue of diffdiff
                    sub_indices4 = A_aux3(sub2ind(size(A_aux3),i_neg,j_neg_first));
                    sub_indices5 = coloncatrld(sub_indices4-j_neg_first+1,sub_indices4-1);
                    invalid_loc = ismember(A_aux3,sub_indices5);
                    condpos_mat(invalid_loc) = 0;
                    
                    mergedEvents_mat = [zeros(size(Rising_diffdiff,1),1),condpos_mat(:,1:end-1)] & condneg_mat;                    
                    mergedEvents_indices = sort(A_aux3(mergedEvents_mat));
                    nNewEvents_vec = sum(mergedEvents_mat,2);
                    
                    
                    %% update spike-counts
                    CountRefineCands_indices_actu = burst_leaders_discrete_indices(nNewEvents_vec~=0);
                    CoundRefineCands_nNewEvents = nNewEvents_vec(nNewEvents_vec~=0);
                    vector_nSpikes_discrete_SR(CountRefineCands_indices_actu) = max(vector_nSpikes_discrete_SR(CountRefineCands_indices_actu) - CoundRefineCands_nNewEvents,1);
                    
                else
                    mergedEvents_indices = [];
                end
                
            else
                mergedEvents_indices = [];
            end
            
            EventTimes_SR = unique([EventTimes_noSR;mergedEvents_indices(:)]);
            idx_onset_Actu_B_SR = [refined_OnsetTime_B; mergedEvents_indices(:)];
            
        else % when opt.demerging == 0
            
            EventTimes_SR = [];
            idx_onset_Actu_B_SR = [];
            
        end % End of "if opt.demerging"
        
    end


%% 
    function  output_UFARSA = build_trains(EventTimes_SR,vector_nSpikes_discrete_SR,idx_onset_Actu_B_SR,ZeroEvent_flag)
        % Reconstruction of spiking activity trains: build the final output-vectors of the reconstruction process
        
        nFrames = opt.nFrames_original;
        
        EventTrain_rec = zeros(1,nFrames);
        EventTrain_rec_SuperResol = zeros(1,nFrames);
        CountTrain_discrete_rec = zeros(1,nFrames);
        CountTrain_discrete_rec_SuperResol = zeros(1,nFrames);
        
        if ZeroEvent_flag==0 % when non-zero number of events could be detected
            shift_time = opt.onset_shift;
            EventTrain_rec(EventTimes_noSR + shift_time ) = 1; % a binary vector of the reconstructed times of all leading and interior events
            EventTrain_rec(nFrames) = 0;
            
            EventTrain_rec_SuperResol(EventTimes_SR + shift_time ) = 1;
            EventTrain_rec_SuperResol(nFrames) = 0;
            
            CountTrain_discrete_rec(refined_OnsetTime_L + shift_time ) = vector_nSpikes_discrete;
            CountTrain_discrete_rec(refined_OnsetTime_B + shift_time ) = 1;
            CountTrain_discrete_rec(nFrames) = 0;
            
            CountTrain_discrete_rec_SuperResol(refined_OnsetTime_L + shift_time ) = vector_nSpikes_discrete_SR;
            CountTrain_discrete_rec_SuperResol(idx_onset_Actu_B_SR + shift_time ) = 1;
            CountTrain_discrete_rec_SuperResol(nFrames) = 0;
        end
        
        output_UFARSA.eTrain = EventTrain_rec(1:nFrames);
        output_UFARSA.eTrain_dem = EventTrain_rec_SuperResol(1:nFrames);
        output_UFARSA.cTrain = CountTrain_discrete_rec(1:nFrames);
        output_UFARSA.cTrain_dem = CountTrain_discrete_rec_SuperResol(1:nFrames);        
        
        if opt.demerging == 0 % when no Demerging was applied
            output_UFARSA.eTrain_dem = [];
            output_UFARSA.cTrain_dem = [];
        end
            
        % fluorescence trace in different versions
        output_UFARSA.fluors.beforeSmoothing = fluors.beforeSmoothing(1:nFrames);
        output_UFARSA.fluors.original = fluors.original(1:nFrames);
        output_UFARSA.fluors.afterSmoothing = fluors.afterSmoothing(1:nFrames);
        output_UFARSA.fluors.normalized = fluors.normalized(1:nFrames);
        
        if opt.gen_FR_count == 1
            output_UFARSA.cFR = gen_smoothed_FR(opt.sigma_gauss, output_UFARSA.cTrain); 
        else
            output_UFARSA.cFR = []; 
        end
        
        if opt.gen_FR_count_dem == 1 && opt.demerging == 1
            output_UFARSA.cFR_dem = gen_smoothed_FR(opt.sigma_gauss, output_UFARSA.cTrain_dem);
        else
            output_UFARSA.cFR_dem = [];
        end
        
    end


end


%%
function [thr_leading,opt] = set_leading_thr(opt,fluors)
    % Setting the leading threshold: determine the threshold for detecting leading events (thr_leading)

if opt.min_leading_amp == 0 % minimum leading threshold will not used to lower-bound the leading threshold (recommended when user estiamted leading threshold from data)
    
    min_thr_leading = 0; % the minimum-leading threshold (i.e. the lower-bound for leading threshold) is removed 
    thr_leading = opt.scale_NoiseSTD*opt.std_noise;  % leading threshold
    
elseif isempty(opt.min_leading_amp) % leading threshold will be estimated, and be lower-bounded with the minimum leading threshold (default)
 
    fluor_main = fluors.beforeSmoothing(1:opt.nFrames_original);
    min_thr_leading = prctile(fluor_main(fluor_main>0),98)/opt.denominator_prct; % minimum leading threshold
    thr_leading = max([min_thr_leading, opt.scale_NoiseSTD*opt.std_noise]); % leading threshold

else

    min_thr_leading = 0; % the minimum-leading threshold (i.e. the lower-bound for leading threshold) is removed 
    thr_leading = opt.min_leading_amp*opt.scale_min_amp; % leading threshold

end

opt.thr_leading = thr_leading;
opt.min_thr_leading_est = min_thr_leading;

end


%%
function [idx_onset_actu_All,idx_rise_mat_All,idx_rise_vec_All,Amps_actu_All,L_indices_aux,Amps_actu_L,ZeroEvent_flag, thr_wburst] = rec_events(fluors,thr_leading,opt)
    % reconstruct the onsets of all leading and interior events (includes "Setting the interior threshold" step)

cond1 = find(diff(fluors.afterSmoothing)>0);
cond2 = find(diff(fluors.afterSmoothing)<=0);
idx_onset = intersect(cond1, cond2+1)';

if (fluors.afterSmoothing(2)> fluors.afterSmoothing(1))
    idx_onset = [1; idx_onset];
end

% find the end frames of all potential rising phases
idx_end = intersect(cond1+1, cond2)';
if numel(idx_onset)>numel(idx_end); idx_onset(end)=[]; end

% compute the amplitudes of the potential rising phases
F_onset = fluors.afterSmoothing(idx_onset)';
F_max   = fluors.afterSmoothing(idx_end)';
Amps     = abs(F_max - F_onset);

nTargets = numel(idx_onset); % number of all potential events

idx_onset_L = idx_onset;
idx_onset_L(Amps < thr_leading) = 0;

if sum(idx_onset_L) ~= 0
    ZeroEvent_flag = 0;
    
    idx_end_L = idx_end;
    idx_end_L(Amps < thr_leading) = 0;
    F_onset_L = F_onset;
    F_onset_L(idx_onset_L==0) = 0;
    Amps_L = Amps;
    Amps_L(idx_onset_L==0) = 0;
    
    idx_onset_actu_L = idx_onset_L(idx_onset_L~=0);
    idx_end_actu_L = idx_end_L(idx_end_L~=0);
    Amps_actu_L = Amps_L(Amps_L~=0);
    
    rec_flag = ones(nTargets,1); % flag of reconstructed events
    rec_flag(idx_onset_L==0) = 0;
    
    idx_onset_B = idx_onset;
    idx_onset_B(idx_onset_L~=0) = 0;
    
    if sum(idx_onset_B)~=0
        
        % Setting the interior threshold step
        thr_wburst = opt.scale_burstAmp * compute_minAmpL( Amps_actu_L, thr_leading, opt );
        thr_wburst =  max( min(opt.scale_min_interior_thr,opt.scale_NoiseSTD) * opt.std_noise, thr_wburst ); % the interior threshold
        
        % check the amplitudes of the candidates for interior events
        idx_onset_B(Amps <thr_wburst) = 0;
        
        Amps_L_aux = Amps_L;
        F_onset_L_aux = F_onset_L;
        indices_new = [];
        indices_old = NaN;
        
        while ~isequal(indices_new, indices_old)
            indices_old = indices_new;
            indices_new = find(([ 0; rec_flag(1:end-1)] & idx_onset_B)~=0);
            rec_flag(indices_new) = 1;
            
            Amps_L_aux(indices_new) = Amps_L_aux(indices_new-1);
            F_onset_L_aux(indices_new) = F_onset_L_aux(indices_new-1);
        end
        
        idx_onset_B(rec_flag==0) = 0;
        Amps_L_aux(idx_onset_B==0) = 0;
        F_onset_L_aux(idx_onset_B==0) = 0;
        
        F_onset_B = F_onset;
        F_onset_B(idx_onset_B==0) = 0;
        
        %
        extra_cond = F_onset_B > opt.plateau_interior*Amps_L_aux + F_onset_L_aux;
        idx_onset_B(extra_cond==0) = 0;
        idx_onset_Actu_B = idx_onset_B(idx_onset_B~=0);
        idx_end_actu_B = idx_end(idx_onset_B~=0);
        Amps_actu_B = Amps(idx_onset_B~=0);
        
    else
        
        idx_end_actu_B = [];
        idx_onset_Actu_B = [];
        Amps_actu_B = [];
        thr_wburst = [];
    end
    
    idx_end_actu_All = [idx_end_actu_L;idx_end_actu_B];
    Amps_actu_All = [Amps_actu_L; Amps_actu_B];
    
    [idx_onset_actu_All, idx_sort_all] = unique([idx_onset_actu_L;idx_onset_Actu_B]);
    idx_end_actu_All = idx_end_actu_All(idx_sort_all);
    Amps_actu_All = Amps_actu_All(idx_sort_all);
    
    max_len = max(idx_end_actu_All-idx_onset_actu_All);
    idx_rise_vec_All = coloncatrld(idx_onset_actu_All, idx_end_actu_All);
    idx_rise_mat_All = reshape(coloncatrld(idx_onset_actu_All, idx_onset_actu_All+max_len),max_len+1,size(idx_onset_actu_All,1))';
    
    L_indices_aux = find(ismember(idx_onset_actu_All,idx_onset_actu_L));
    
else
    
    ZeroEvent_flag = 1;
    thr_wburst = [];
    idx_onset_actu_All = [];
    idx_rise_mat_All = [];
    idx_rise_vec_All = [];
    Amps_actu_All = [];
    L_indices_aux = [];
    Amps_actu_L = [];
    
end

end


%%
function  [refined_OnsetTime,Frise_smth_orig_L,A_aux_L,refined_OnsetTime_L, refined_OnsetTime_B] = timing_correction(idx_onset_actu_All,idx_rise_mat_All,idx_rise_vec_All, fluors, Amps_actu_All, L_indices_aux )
    % Timing correction: refine the reconstructed times of both leading and interior events (except the demerged ones)
    
nCols = size(idx_rise_mat_All,2);
nRows = size(idx_rise_mat_All,1);
rows = [1:nRows]';

A = idx_rise_mat_All;
A_vec = reshape(A', 1, nCols*nRows);
[~, z] = unique(A_vec, 'last');
extra_indices = 1:nCols*nRows;
extra_indices(z) = [];
A_vec(extra_indices) = NaN;
A = reshape(A_vec, nCols, nRows)';
A(A> numel(fluors.beforeSmoothing)) = NaN;

idx_rise_mat_All(idx_rise_mat_All> numel(fluors.beforeSmoothing)) = numel(fluors.beforeSmoothing);

nan_loc = ~ismember(A,idx_rise_vec_All);
A_aux = A;
A_aux(nan_loc) = NaN;
A_aux_L = A_aux(L_indices_aux,:); % we need this for demerging step


% find the onset based on noisy (but pre-processed) trace:
Frise_noisy = fluors.beforeSmoothing(idx_rise_mat_All);
Frise_noisy(nan_loc) = NaN;
Frise_smth = fluors.afterSmoothing(idx_rise_mat_All);
Frise_smth(nan_loc) = NaN;
Frise_smth_orig_L = Frise_smth(L_indices_aux,:);
Actualmin_NoisyRising = repmat(nanmin(Frise_noisy,[],2), 1, nCols);
min_zero_NoisyRising = Frise_noisy - Actualmin_NoisyRising;
min_zero_SmoothedRising = Frise_smth - Actualmin_NoisyRising;
[~, max_diff_idx] = nanmax(diff(min_zero_SmoothedRising,[],2),[],2);


sub_indices = coloncatrld(idx_rise_mat_All(:,1),idx_rise_mat_All(sub2ind(size(idx_rise_mat_All),rows,max_diff_idx)));
nan_loc2 = ismember(idx_rise_mat_All,sub_indices);
Frise_smth(nan_loc2) = NaN;
[~, maxNoisy_idx] = nanmax(Frise_smth,[],2);
idx_max_noisy_transient = maxNoisy_idx;


Thr_LowerRange = min_zero_SmoothedRising(:,1) + 0.25*Amps_actu_All;
min_zero_SmoothedRising(min_zero_SmoothedRising>=repmat(Thr_LowerRange,1,nCols))= NaN;
MedianLow_intensity = nanmedian(min_zero_SmoothedRising,2);
Thr_LowerRange_refined = MedianLow_intensity + 0.25*Amps_actu_All;
RisingPhase_meanFree = min_zero_NoisyRising - repmat(Thr_LowerRange_refined,1,nCols);

RisingPhase_meanFree_orig = RisingPhase_meanFree;
RisingPhase_meanFree(RisingPhase_meanFree>=0) = NaN;

[i_end j_end] = find(isnan(A_aux));
[i_end idx_sort] = sort(i_end);
j_end = j_end(idx_sort);

if ~isempty(i_end)
    
    i1_first_All = [1 find(diff(i_end')>0)+1];
    i1_last_All = j_end(i1_first_All) - 1;
    rows_nan = unique(i_end);
    idx_rise_end = idx_rise_mat_All(:,end);
    idx_rise_end(rows_nan) = idx_rise_mat_All(sub2ind(size(idx_rise_mat_All),rows_nan,i1_last_All));
    
else    
    idx_rise_end = idx_rise_mat_All(:,end);
end

sub_indices2 = coloncatrld( idx_rise_mat_All(sub2ind(size(idx_rise_mat_All),rows,idx_max_noisy_transient)), idx_rise_end);
nan_loc3 = ismember(A,sub_indices2);
RisingPhase_meanFree(nan_loc3) = NaN;

idx_rise_mat_aux = idx_rise_mat_All;
idx_rise_mat_aux(isnan(RisingPhase_meanFree)) = NaN;
[~, j_yes] = max(idx_rise_mat_aux,[],2);
cols_aux = idx_rise_mat_All(sub2ind(size(idx_rise_mat_All),rows,j_yes+1));
empty_rows = sum(~isnan(idx_rise_mat_aux),2)==0;
cols_aux(empty_rows) = 0;
sub_indices3 = coloncatrld( idx_rise_mat_All(:,1), cols_aux );
rem_loc = ismember(A,sub_indices3);

RisingPhase_meanFree_orig(~rem_loc) = NaN;
diff_rise = diff(RisingPhase_meanFree_orig,[],2);
diff_rise(diff_rise>0)=NaN;
diff_rise(diff_rise<=0)=1; 
diff_rise = diff_rise(:,end:-1:1);

empty_rows2 = sum(~isnan(diff_rise),2)==0;
[~, i_max] = max(diff_rise,[],2);
i_max(empty_rows2) = NaN;
i_onset = nCols - i_max + 1;
refined_OnsetTime = idx_onset_actu_All;

refined_OnsetTime(~isnan(i_onset)) = idx_rise_mat_All( sub2ind(size(idx_rise_mat_All), rows(~isnan(i_onset)),i_onset(~isnan(i_onset)) ) );
refined_OnsetTime_L = refined_OnsetTime(L_indices_aux);

refined_OnsetTime_B = refined_OnsetTime(ismember(refined_OnsetTime,refined_OnsetTime_L)==0);

if numel(refined_OnsetTime)~=numel(idx_onset_actu_All)
    error('there is a bug, please contact the author (Vahid Rahmati)')
end

end


%% 
function  [burst_leaders_discrete_indices, vector_nSpikes_discrete] = spikecount_estimation(opt, Amps_actu_L,thr_leading)
    % Spike-count estimation: compute the relative spike-counts per detected leading fluorescence transient

min_amp = compute_minAmpL( Amps_actu_L, thr_leading, opt ); % estimate the minimum interior amplitude 
frac_amp = Amps_actu_L./min_amp;
frac_amp_floor = floor(frac_amp);
vector_nSpikes_discrete = frac_amp_floor + [(frac_amp-frac_amp_floor)>=opt.round_border];

burst_leaders_discrete_indices = find(vector_nSpikes_discrete>1);
vector_nSpikes_discrete(vector_nSpikes_discrete==0) = 1;

end


%%
function min_amp = compute_minAmpL( Amps_actu_L, thr_leading, opt )
    % estimate the minimum amplitude of interior transients (i.e. the minimum interior amplitude) 
    
min_amp = mean(Amps_actu_L( Amps_actu_L < median(Amps_actu_L) ));

if min_amp ==0 || isnan(min_amp)
    min_amp = mean(Amps_actu_L);
end

max_nBurstSpikes = opt.max_nBurstSpikes; % the maximum spike-count constant 

if ~isempty(max_nBurstSpikes) % when a maximum spike-count was set by user
    if (max(Amps_actu_L)/min_amp) > max_nBurstSpikes % true, when there are some spike-counts exceeding the maximum spike-count
        
        % re-compute min_amp so that the mximum spike-count does not exceed the max_nBurstSpikes 
        min_amp = max(Amps_actu_L)/max_nBurstSpikes;
        
    end
end

min_amp = max(thr_leading, min_amp); % the estimated minimum interior amplitude

end


%%
function  smoothed_FR_trace = gen_smoothed_FR(sigma_gauss, spiking_activity_train)
    % generate the estimated firing rate vector based on the given spiking activity train
    
% Time points 
t = 0:1:length(spiking_activity_train)-1; 

if sigma_gauss<1 % a sigma_gauss may not be meaningful, if it is lower than temporal resolution of fluorescence data 
    warning('sigma_gauss was <1 frame, thus it was set internally to 1')
    sigma_gauss = 1; % in [frame]
end

% define the Gaussian kernel used for generating the estimated firing rate vector 
h_Gauss = gaussmf( -0.5*t(end):1:0.5*t(end),[sigma_gauss 0]); 

% convolve the spike-count train with the Gaussian kernel to produce the estimated firing rate vector 
smoothed_FR_trace = conv(h_Gauss,spiking_activity_train);

% Align the firing rate trace to the spiking activity train
smoothed_FR_trace = smoothed_FR_trace(ceil(length(h_Gauss)/2): (end-floor(length(h_Gauss)/2)) );

end