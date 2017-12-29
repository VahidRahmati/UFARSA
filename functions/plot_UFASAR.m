function h_fig = plot_UFASAR(output_UFASAR,opt, true_spikes)

% h_fig = plot_UFASAR(output_UFASAR,opt, true_spikes)
%
% This function plots the UFASAR's reconstruction results as well as different versions  
% of the processed fluorescence trace. For a complete list of plotting options, and saving
% them please see the "internal_parameters.m" file or the "user_guide.pdf" file. 
%
% INPUT: 
%   output_UFASAR: an structure variable containing UFASAR's reuslts (see the description of run_UFASAR.m function)
%
%   opt: an structure variable containing the plotting options (for a complete list see the description of run_UFASAR.m function)
%        opt.plt  --> 1: plot the results, 0: skip plotting results
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
%   true_spikes: the "simulated" train of spiking event times or spike-counts (the time-unit is [frame])
%
% OUTPUT:
%   h_fig: the current figure handle.
%
% Author: Vahid Rahmati (December, 2017)


%% initialize the figure and check the options
try
    count_train_sim = true_spikes;
    sim_flag = 1;
catch
    sim_flag = 0;
end

if opt.plt_cTrain && opt.plt_cTrain_dem && opt.demerging && opt.ZeroEvent_flag==0
   text1 = '#Warning A#: Currently, the code cannot plot simultaneously both rec. spike-count trains with and without demerging';
   text2 = '\nin the same figure. Therefore, only the reconstructed demerged spike-counts are shown.';
   if opt.warn_flag; display(sprintf([text1,text2])), end
   opt.plt_cTrain = 0;
end

if opt.plt_eTrain && opt.plt_eTrain_dem && opt.demerging && opt.ZeroEvent_flag==0
   text1 = '#Warning B#: Currently, the code cannot plot simultaneously both rec. event trains with and without demerging';
   text2 = '\nin the same figure. Therefore, only the reconstructed demerged event times are shown.';
   if opt.warn_flag; display(sprintf([text1,text2])), end
   opt.plt_eTrain = 0;
end

if opt.plt_cFR && opt.plt_cFR_dem && opt.demerging && opt.ZeroEvent_flag==0
   text1 = '#Warning C#: Currently, the code cannot plot simultaneously both estimated FR vectors with and without demerging';
   text2 = '\nin the same figure. Therefore, only the estimated FR vector under demerging case has been shown.';
   if opt.warn_flag; display(sprintf([text1,text2])), end
   opt.plt_cFR = 0;
end

col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0], [86 180 233],...
     [204 121 167], [64 224 208], [200 208 0], [240 228 66], [204 204 204], [229 191 0],  [255 20 200]*0.92, [1 0.4 0.6]*255};

height_p = 0.16;
width_p = 0.98/2.2; 
leftDiff_p = 0.038;

set(0, 'Units' , 'normalized' )
h_fig = figure;
scale_width = 0.9;
scale_heigth = 0.6;
set(gcf,'Units','normalized','position',[0,0,scale_width,scale_heigth])

set(gcf, 'color', 'w', ...
    'defaultAxesFontSize', 15, ...
    'defaultlinelinewidth',2, ...
    'Name', ['ROI (',num2str(opt.iROI),')'],'NumberTitle','off');


%% plot different versions of the given fluorescence trace
if (opt.plt_cTrain || opt.plt_cTrain_dem || opt.plt_cFR || opt.plt_cFR_dem  || opt.plt_cTrain_sim) && opt.ZeroEvent_flag==0
    axes('Units', 'normalized','position', [leftDiff_p, .53, width_p*2.12, height_p*2.25]);
    set(gca,'xtick',[])
else
    axes('Units', 'normalized','position', [leftDiff_p, .13, width_p*2.12, height_p*5]);
end
alpha(.7);

hold on
legends = {};
plt_idx = 0;
plt_notAssigned = 0;
if (opt.plt_Ft_orig_norm + opt.plt_Ft_beforeSmth + opt.plt_Ft_afterSmth)==0
    plt_notAssigned = 1;
end

orig_f = 0;
if opt.plt_Ft_orig_norm || plt_notAssigned % plot the original, normalized fluorescence trace
  plt_idx = plt_idx + 1;
  h(plt_idx) = plot(1:opt.nFrames_original, output_UFASAR.fluors.normalized,'color', col{8}/255);
  max_fluor_orig_norm = max(output_UFASAR.fluors.normalized);
  min_fluor_orig_norm = min(output_UFASAR.fluors.normalized);
  orig_f = 1;
  legends = [legends 'Norm. fluor'];
end

pre_f = 0;
if opt.plt_Ft_beforeSmth % plot the pre-processed fluorescence trace used for smoothing step
  plt_idx = plt_idx + 1;
  h(plt_idx) =plot(1:opt.nFrames_original, output_UFASAR.fluors.beforeSmoothing,'color', col{10}/255);
  max_fluor_beforeSmoothing = max(output_UFASAR.fluors.beforeSmoothing);
  min_fluor_beforeSmoothing = min(output_UFASAR.fluors.beforeSmoothing);
  pre_f = 1;
  legends = [legends 'Pre-. fluor'];
end

if opt.plt_Ft_afterSmth % plot the smoothed fluorescence trace used for reconstruction step
   plt_idx = plt_idx + 1; 
  h(plt_idx) =plot(1:opt.nFrames_original, output_UFASAR.fluors.afterSmoothing,'color', col{3}/255);
  legends = [legends 'Smth. fluor'];
end


if orig_f && pre_f
    max_fluor_plt = max([max_fluor_orig_norm,max_fluor_beforeSmoothing]);
    min_fluor_plt = min([min_fluor_orig_norm,min_fluor_beforeSmoothing]);
elseif orig_f 
    max_fluor_plt = max([max_fluor_orig_norm]);
    min_fluor_plt = min([min_fluor_orig_norm]);
elseif pre_f || smth_f
    max_fluor_plt = max([max_fluor_beforeSmoothing]);
    min_fluor_plt = min([min_fluor_beforeSmoothing]);
end

%% plot the reconstructed and/or simulated spiking activity trains
if opt.plt_eTrainOnFt % plot the user's selected type of reconstructed event train on the fluorescence trace(s)
     if opt.plt_eTrain 
        evet_times_rec = find(output_UFASAR.eTrain>0);
        plot(evet_times_rec,max_fluor_plt*ones(1,numel(evet_times_rec)), 'x', 'LineWidth',2, 'MarkerEdgeColor',col{12}/255,'MarkerSize',10); 
     elseif opt.plt_eTrain_dem
        evet_times_rec = find(output_UFASAR.eTrain_dem>0);
        plot(evet_times_rec, max_fluor_plt*ones(1,numel(evet_times_rec)), 'x', 'LineWidth',2, 'MarkerEdgeColor',col{12}/255,'MarkerSize',10);   
     end
    clear evet_times_rec
end

set(gca,'TickDir','out');
xlim([1 opt.nFrames_original])
ylim([min_fluor_plt 1.2*max_fluor_plt])

size_tick = get(gca,'TickLength');
set(gca,'TickLength',size_tick*0.5)

hold off

if (opt.plt_cTrain || opt.plt_cTrain_dem || opt.plt_cFR || opt.plt_cFR_dem  || opt.plt_cTrain_sim) && opt.ZeroEvent_flag==0
    
    axes('position', [leftDiff_p+0.0, .13, width_p*2.12, height_p*2.25]);
    hold on
    
    if opt.plt_cTrain % plot reconstructed spike-count train without demerging
        count_times_rec = find(output_UFASAR.cTrain>0);
        for ii = 1:numel(count_times_rec)
            if ii==1
                plt_idx = plt_idx + 1;
                h(plt_idx) = line([count_times_rec(ii) count_times_rec(ii)],...
                    [0 output_UFASAR.cTrain(count_times_rec(ii))], 'color', col{5}/255, 'linewidth', 1.5);
            else
                h(plt_idx) = line([count_times_rec(ii) count_times_rec(ii)],...
                    [0 output_UFASAR.cTrain(count_times_rec(ii))], 'color', col{5}/255, 'linewidth', 1.5);
            end
        end
        legends = [legends 'Rec. count'];
    end
    alpha(.7);
        
    if opt.plt_cTrain_dem && opt.demerging % plot the reconstructed spike-count train with demerging
        count_times_rec_SR = find(output_UFASAR.cTrain_dem>0);
        for ii = 1:numel(count_times_rec_SR)
            if ii==1
                plt_idx = plt_idx + 1;
                h(plt_idx) = line([count_times_rec_SR(ii) count_times_rec_SR(ii)],...
                    [0 output_UFASAR.cTrain_dem(count_times_rec_SR(ii))], 'color', col{5}/255, 'linewidth', 1.5);
            else
                h(plt_idx) = line([count_times_rec_SR(ii) count_times_rec_SR(ii)],...
                    [0 output_UFASAR.cTrain_dem(count_times_rec_SR(ii))], 'color', col{5}/255, 'linewidth', 1.5);
            end
        end
        legends = [legends 'Rec. count (dem)'];
        
    end
    
    
    if opt.plt_eTrain_dem && opt.demerging % plot the reconstructed event train with demerging
        evet_times_rec = find(output_UFASAR.eTrain_dem>0);
        plt_idx = plt_idx + 1;
        h(plt_idx) = plot(evet_times_rec, ones(1,numel(evet_times_rec)), 'x', 'LineWidth',2, 'MarkerEdgeColor',col{12}/255,'MarkerSize',10);
        legends = [legends 'Rec. event (dem)'];
    end
    
    if opt.plt_eTrain % plot the reconstructed event train without demerging
        clear evet_times_rec
        evet_times_rec = find(output_UFASAR.eTrain>0);
        plt_idx = plt_idx + 1;
        h(plt_idx) = plot(evet_times_rec,ones(1,numel(evet_times_rec)), 'x', 'LineWidth',2, 'MarkerEdgeColor',col{12}/255,'MarkerSize',10);
        legends = [legends 'Rec. event'];
    end
    
    if opt.plt_cFR && opt.gen_FR_count % plot the estimated firing rate vector, based on reconstructed spike-count without demerging
        plt_idx = plt_idx + 1;
        h(plt_idx)= plot(1:opt.nFrames_original,output_UFASAR.cFR, 'LineWidth',2, 'Color',col{2}/255);
        legends = [legends 'Rec. FR'];
    end
    
    if opt.plt_cFR_dem && opt.demerging && opt.gen_FR_count_dem % plot the estimated firing rate vector, based on reconstructed spike-counts with demerging
        plt_idx = plt_idx + 1;
        h(plt_idx)= plot(1:opt.nFrames_original,output_UFASAR.cFR_dem, 'LineWidth',2, 'Color',col{2}/255);
        legends = [legends 'Rec. FR (dem)'];
    end
    
    
    if opt.plt_cTrain_sim && sim_flag % plot the simulated spike-count train
        event_times_sim = find(true_spikes>0);
        plt_idx = plt_idx + 1;
        h(plt_idx) = plot(event_times_sim,true_spikes(event_times_sim), 'o', 'LineWidth',2, 'MarkerEdgeColor',col{4}/255,'MarkerSize',10);
        legends = [legends 'True. count'];
    end
    
    if opt.plt_eTrain_sim && sim_flag % plot the simulated spike-count train
        event_times_sim = find(true_spikes>0);
        plt_idx = plt_idx + 1;
        h(plt_idx) = plot(event_times_sim,ones(1,numel(event_times_sim)), 's', 'LineWidth',2, 'MarkerEdgeColor','k','MarkerSize',15);
        legends = [legends 'True. event'];
    end
    
    hold off
    
    xlim([1 opt.nFrames_original])
    axis tight
    set(gca,'TickDir','out');
    
    ylim_2 = get(gca,'ylim');
    ylim(ylim_2*1.08)
    
elseif (opt.plt_cTrain || opt.plt_cTrain_dem || opt.plt_cFR || opt.plt_cFR_dem  || opt.plt_cTrain_sim) && opt.ZeroEvent_flag==1
   
   if opt.warn_flag; display(sprintf('#Warning D#: No reconstruction was shown, as no event was reconstructed')), end
    
end % End of plotting spiking activities

leg_g = legend(h,legends,'FontSize',13,'orientation','horizental');
set(leg_g,'Units','normalized');
leg_pos = get(leg_g,'Position');
set(leg_g,'Position',[ 0.01 0.95  leg_pos(3) leg_pos(4)])

set(leg_g, 'box', 'off')

set(gca,'TickLength',size_tick*0.5)
xl = xlabel('Time [frame]','FontSize',15);
set(xl,'Units','normalized');
xl_pos = get(xl, 'position');
set(xl,'Units', 'normalized','Position',[xl_pos(1), 0.95*xl_pos(2), xl_pos(3)])

movegui(h_fig,'center')

hold off
