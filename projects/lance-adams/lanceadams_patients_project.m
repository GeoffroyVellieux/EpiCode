% analysis for Lance-Adams patients' EEG

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end

ft_defaults

config = lanceadams_patients_setparams;

%% read marker, verify markers, and add merged markers
for ipatient = 1:size(config, 2)
    
    %read muse markers
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, true);
    
    %check that there are no error in markers
    wrong_markers{ipatient} = verifymarkers(config{ipatient}, MuseStruct{ipatient}, 1);
    
    % create a marker "ied_all" with all ied markers
    cfgtemp = [];
    cfgtemp.newmarker = 'ied_all';
    cfgtemp.markers_to_merge = {'Polyspike', 'Polyspike_wave', 'Spike', 'Spike_wave'};
    MuseStruct{ipatient} = merge_MuseMarkers(cfgtemp, MuseStruct{ipatient});
    
    % create a marker "emg_all" with all EMG markers
    cfgtemp = [];
    cfgtemp.newmarker = 'emg_all_start';
    cfgtemp.markers_to_merge = {'Myoclonus_start', 'EMGnoEEG_start'};
    MuseStruct{ipatient} = merge_MuseMarkers(cfgtemp, MuseStruct{ipatient});
    
    cfgtemp = [];
    cfgtemp.newmarker = 'emg_all_end';
    cfgtemp.markers_to_merge = {'Myoclonus_end', 'EMGnoEEG_end'};
    MuseStruct{ipatient} = merge_MuseMarkers(cfgtemp, MuseStruct{ipatient});
        
end

%% data length without artefacts
for ipatient = 1:size(config, 2)
    %go through each eeg of this patient
    for idir = 1:size(MuseStruct{ipatient}{1}, 2)
        %data
        data_start  = MuseStruct{ipatient}{1}{idir}.starttime;
        data_end    = MuseStruct{ipatient}{1}{idir}.endtime;
        data_length{ipatient}(idir) = seconds(data_end - data_start);
        %artefacts
        bad_start{ipatient}{idir}  = ft_getopt(MuseStruct{ipatient}{1}{idir}.markers.BAD__START__, 'synctime', []);
        bad_end{ipatient}{idir}    = ft_getopt(MuseStruct{ipatient}{1}{idir}.markers.BAD__END__, 'synctime', []);
        bad_length{ipatient}(idir) = sum(bad_end{ipatient}{idir} - bad_start{ipatient}{idir});
    end
end

%% frequency of events

data = table.empty;
irow = 0;

for ipatient = 1:size(config, 2)
    
    %go through each eeg of this patient
    for idir = 1:size(MuseStruct{ipatient}{1}, 2)
        %compute and store values for each marker
        for markername = ["Polyspike", "Polyspike_wave", "Spike", "Spike_wave", "ied_all", "Myoclonus", "EMGnoEEG", "emg_all"]
            
            irow = irow + 1;
            
            data.pat_ID{irow} = config{ipatient}.prefix(1:end-1);
            data.EEG_ID{irow} = config{ipatient}.directorylist{1}{idir};
            data.markername(irow) = markername;
            
            %freq
            n_events = numel(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime);
            t        = data_length{ipatient}(idir) - bad_length{ipatient}(idir);
            data.n_events(irow) = n_events;
            data.nb_per_min(irow) = n_events/t * 60;
            
            %cv of inter event intervals
            events_intervals = diff(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime);
            toremove = [];
            %remove interval which starts just before the artefact
            for ibad = 1:size(bad_start{ipatient}{idir}, 2)
                toremove(end+1) = find(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime < bad_start{ipatient}{idir}(ibad), 1, 'last');
                toremove(end+1) = find(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime < bad_end{ipatient}{idir}(ibad), 1, 'last');
            end
            toremove = unique(toremove);
            events_cv = events_intervals;
            events_cv(toremove) = [];
            data.cv(irow) = std(events_cv) / mean(events_cv);
            
            %cv2 of inter events intervals
            cv2_data = cv2(events_intervals);
            % remove intervals with artefacts
            toremove = [];
            for ibad = 1:size(bad_start{ipatient}{idir}, 2)
                toremove(end+1) = find(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime < bad_start{ipatient}{idir}(ibad), 1, 'last');
                toremove(end+1) = find(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime < bad_end{ipatient}{idir}(ibad), 1, 'last');
            end
            toremove = unique(toremove);
            toremove = [toremove, toremove+1];
            toremove = sort(toremove);
            for i=1:2
                if toremove(end) > length(cv2_data)
                    toremove = toremove(1:end-1);
                end
            end
            cv2_data(toremove) = [];
            data.cv2(irow) = mean(cv2_data);
        end
    end
end

tablename = fullfile(config{ipatient}.tablesavedir, 'allpatients-events_frequency.xlsx');
writetable(data, tablename);

%% plot event occurences of 15 randomly-chosen minutes
ipatient = 1;
idir = 1;
markername = "ied_all";
fig = figure;
subplot(5,1,1); hold on;
%plot events
for ievent= 1:size(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).clock,2)
    t = MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).clock(ievent);
    plot([t,t],[-1,1],'-k','LineWidth',2);
end
%plot artefacts
x = xlim;
for ibad = 1:size(MuseStruct{ipatient}{1}{idir}.markers.BAD__START__.clock,2)
    y = ylim;
    s = MuseStruct{ipatient}{1}{idir}.markers.BAD__START__.clock(ibad);
    e = MuseStruct{ipatient}{1}{idir}.markers.BAD__END__.clock(ibad);
    xpatch = [s e e s];
    ypatch = [y(1) y(1) y(2) y(2)];
    patch('XData',xpatch,'YData',ypatch,'facecolor',[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
end

%select a random 15 minutes
plot_length = minutes(diff(xlim));
i = rand * (round(plot_length) - 15);
xlim([x(1)+minutes(i) x(1)+minutes(i+15)]);

%plot settings
yticks([]);
ylim([-2 2]);
set(gca, 'FontSize', 15);
xlabel('Time');
fname = fullfile(config{ipatient}.imagesavedir, 'examples_of_ied_occurences', sprintf('%s_dir%d_markername_examples_of_occurence', config{ipatient}.prefix, idir, markername));
dtx_savefigure(fig,fname,'pdf','png','close');

%% EMG timings

data = table.empty;
irow = 0;

for ipatient = 1:size(config, 2)
    
    %go through each eeg of this patient
    for idir = 1:size(MuseStruct{ipatient}{1}, 2)
        
        for markername = ["Myoclonus", "EMGnoEEG", "emg_all"]
            irow = irow+1;
            data.pat_ID{irow} = config{ipatient}.prefix(1:end-1);
            data.EEG_ID{irow} = config{ipatient}.directorylist{1}{idir};
            data.markername{irow} = markername;
            
            data.n_events(irow) = size(MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime, 2);
            t        = data_length{ipatient}(idir) - bad_length{ipatient}(idir);
            data.nb_per_min(irow) = data.n_events(irow)/t * 60;
            
            
            emgstart = MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.startmarker.(markername)).synctime;
            emgend   = MuseStruct{ipatient}{1}{idir}.markers.(config{ipatient}.muse.endmarker.(markername)).synctime;
            
            emg_duration = emgend - emgstart;
            data.emg_duration_mean(irow) = mean(emg_duration);
            data.emg_duration_std(irow)  = std(emg_duration);
            
            if markername == "Myoclonus"
                t_emg = MuseStruct{ipatient}{1}{idir}.markers.Myoclonus_start.synctime;
                for ievent = 1:size(t_emg, 2)
                    %fixme : nearest ?
                    idx = find(MuseStruct{ipatient}{1}{idir}.markers.ied_all.synctime < t_emg(ievent), 1, 'last');
                    t_eeg(ievent) = MuseStruct{ipatient}{1}{idir}.markers.ied_all.synctime(idx);
                end
                
                eeg_emg_delay = t_emg - t_eeg;
                
                data.eeg_emg_delay_mean(irow) = mean(eeg_emg_delay);
                data.eeg_emg_delay_std(irow) = std(eeg_emg_delay);
            else
                data.eeg_emg_delay_mean(irow) = nan;
                data.eeg_emg_delay_std(irow) = nan;
            end
        end
    end
end

tablename = fullfile(config{ipatient}.tablesavedir, 'allpatients-emg_durations.xlsx');
writetable(data, tablename);

%% plot raw data of each marker, before alignment
for ipatient = 1:size(config, 2)
    %readLFP
    %fixeme removeme : fix an error in marker annotation
    config{ipatient}.muse.endmarker.Myoclonus = 'Myoclonus_start';
    config{ipatient}.muse.endmarker.emg_all = 'emg_all_start';
    %
    LFP = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
    
    %remove artefacted trials
    for markername = string(fieldnames(LFP{1}))'
        cfgtemp        = [];
        cfgtemp.trials = ~LFP{1}.(markername).trialinfo.BAD_cnt;
        LFP{1}.(markername) = ft_selectdata(cfgtemp, LFP{1}.(markername));
    end
    
    %plot chaque pattern
    for markername = string(fieldnames(LFP{1}))'
        
        cfgtemp         = [];
        cfgtemp.channel = config{1}.align.channel{:};
        data            = ft_selectdata(cfgtemp, LFP{1}.(markername));
        
        fig = figure; hold on;
        for itrial = 1:size(data.trial, 2)
            p = plot(data.time{itrial}, data.trial{itrial}, 'k');
            p.Color(4) = 0.2;
        end
        
        dataavg = ft_timelockanalysis([], data);
        plot(dataavg.time, dataavg.avg, 'r', 'linewidth', 2);
        
        xlim([-1 1])
        title(sprintf('%s \n%s : %d trials', config{ipatient}.prefix(1:end-1), markername, size(data.trial, 2)), 'fontsize', 20, 'interpreter', 'none');
        ylabel(sprintf('%s (µV)', config{1}.align.channel{:}));
        xlabel('Time (sec. from Muse marker)');
        set(gca, 'tickdir', 'out', 'fontsize', 15)
        
        figname = fullfile(config{ipatient}.imagesavedir, 'events_raw', sprintf('%s%s_raw', config{ipatient}.prefix, markername));
        %     dtx_savefigure(fig, figname, 'png', 'close');
        
    end
end

%% align data 
clear LFP
for ipatient = 1:size(config, 2)
    
    %perform the alignment
    MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
    
    %fixeme removeme : fix an error in marker annotation
    config{ipatient}.muse.endmarker.Myoclonus = 'Myoclonus_start';
    config{ipatient}.muse.endmarker.emg_all = 'emg_all_start';
    %
    %read LFP and cut data according to the aligned markers
    LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
end

%% plot aligned data
for ipatient = 1:size(config, 2)
       
    %plot each pattern
    for markername = ["Polyspike", "Polyspike_wave", "Spike", "Spike_wave", "ied_all"]
        
        %check that artefacts were well removed
        assert(sum(LFP{ipatient}{1}.(markername).trialinfo.BAD_cnt) == 0);
        
        cfgtemp         = [];
        cfgtemp.channel = config{1}.align.channel{:};
        data            = ft_selectdata(cfgtemp, LFP{ipatient}{1}.(markername));
        
        fig = figure; hold on;
        for itrial = 1:size(data.trial, 2)
            p = plot(data.time{itrial}, data.trial{itrial}, 'k');
            p.Color(4) = 0.2;
        end
        
        dataavg = ft_timelockanalysis([], data);
        plot(dataavg.time, dataavg.avg, 'r', 'linewidth', 2);
        
        title(sprintf('%s \n%s : %d trials', config{ipatient}.prefix(1:end-1), markername, size(data.trial, 2)), 'fontsize', 20, 'interpreter', 'none');
        ylabel(sprintf('%s (µV)', config{1}.align.channel{:}));
        xlabel('Time (sec. from Muse marker)');
        set(gca, 'tickdir', 'out', 'fontsize', 15)
        
        figname = fullfile(config{ipatient}.imagesavedir, 'events_aligned', sprintf('%s%s_alignedXcorr', config{ipatient}.prefix, markername));
        %     dtx_savefigure(fig, figname, 'png', 'close');
        
    end
    
    % plot average of each pattern
    fig = figure; 
    sgtitle(config{ipatient}.prefix(1:end-1), 'interpreter', 'none', 'fontsize', 28, 'fontweight', 'bold');
    iplot = 0;
    for markername = ["Polyspike", "Polyspike_wave", "Spike", "Spike_wave"]
        
        iplot = iplot+1;
        subplot(2,2,iplot); hold on;
        title(markername, 'fontsize', 20, 'interpreter', 'none');
        
        cfgtemp         = [];
        cfgtemp.channel = config{1}.align.channel;
        data            = ft_selectdata(cfgtemp, LFP{ipatient}{1}.(markername));
        
        dataavg = ft_timelockanalysis([], data);
        patch_std(dataavg.time, dataavg.avg, sqrt(dataavg.var), 'k');
        plot(dataavg.time, dataavg.avg, 'k', 'linewidth', 2);
        
        ylabel(sprintf('%s (µV)', config{1}.align.channel{:}));
        xlabel('Time (sec. from Muse marker)');
        set(gca, 'tickdir', 'out', 'fontsize', 15)
        axis tight;
        xlim(config{ipatient}.align.latency.(markername));
        
    end
    
    figname = fullfile(config{ipatient}.imagesavedir, 'events_aligned', sprintf('%sallmarkers_aligned', config{ipatient}.prefix));
    dtx_savefigure(fig, figname, 'png', 'close');
end
    
%% plot avg of all channels
for ipatient = 1:size(config, 2)
    
    for markername = string(config{ipatient}.LFP.name)
        
        dataavg = ft_timelockanalysis([], LFP{ipatient}{1}.(markername));
        
        fig = figure; hold on
        h = 50;
        %plot EEG
        ichan = 0;
        for channame = string(config{1}.LFP.channel )
            
            ichan = ichan+1;
            
            %select channel
            cfgtemp = [];
            cfgtemp.channel = char(channame);
            data_1chan = ft_selectdata(cfgtemp,dataavg);
            plot(data_1chan.time,data_1chan.avg+(numel(dataavg.label)+1)*h-h*ichan,'k'); %first on top
        end
       
        axis tight
        ylim([-h (numel(dataavg.label)+2)*h]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('Channel name', 'Fontsize',15);
        tick = h;
        yticks(0 : tick : (numel(config{ipatient}.LFP.channel))*h);
        %chan_name_all = [string(config{ipatient}.LFP.channel), string(EMG.label)];
        chan_name_all = string(config{ipatient}.LFP.channel);
        set(gca, 'YTickLabel',flip(chan_name_all));
        set(gca, 'FontWeight','bold', 'Fontsize',15, 'FontName', 'Arial');
        set(gca,'TickDir','out');
        
        title(sprintf('%s : %s',config{ipatient}.prefix(1:end-1), markername),'Interpreter','none','Fontsize',18);
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir, 'events_aligned',sprintf('%s%s_allchannels',config{ipatient}.prefix,markername));
        %     dtx_savefigure(fig, fname, 'png', 'pdf', 'close');
    end
end

%% ied topography

for markername = string(config{1}.name)
    
    fig = figure;
    sgtitle(markername, 'interpreter', 'none');
    iplot = 0;
    for ipatient = 1:size(config,2)
        fprintf('topoplot %s, patient %d\n', markername, ipatient);
        
        cfgtemp = [];
        cfgtemp.channel = config{ipatient}.topoplot.channel;
        to_plot = ft_timelockanalysis([], LFP{ipatient}{1}.(markername));
        toi = config{ipatient}.topoplot.toi;
        
        iplot = iplot+1;
        subplot(ceil(size(config,2)/3), min([size(config,2), 3]), iplot); hold on;
        cfgtemp = [];
        cfgtemp.layout        = 'EEG1020';
        cfgtemp.colorbar      = 'yes';
        cfgtemp.zlim          = 'zeromax';
        cfgtemp.xlim          = toi;
        cfgtemp.comment       = 'xlim';
        cfgtemp.fontsize      = 15;
        cfgtemp.renderer      = 'painters';
        cfgtemp.colormap      = jet(5000); %flip so negative peak is in yellow
        cfgtemp.comment = 'no';
        ft_topoplotER(cfgtemp,to_plot);
        
        c = caxis;
        if c(2) < 20
            caxis([c(1) 20]);
        end
        
        title(sprintf('pat %d (n=%d)',ipatient, size(LFP{ipatient}{1}.(markername).trial,2)), 'Interpreter','none');
    end
    
    %print to file
    fname = fullfile(config{ipatient}.imagesavedir, 'events_topography',sprintf('allpatients_topoplot_%s',markername));
%     dtx_savefigure(fig,fname,'pdf','png','close');

end

%% EMG : compute envelopes
clear env*
for markername = ["Myoclonus", "EMGnoEEG", "emg_all", "ied_all"]
    
    for ipatient = 1:size(config,2)
        
        cfgtemp         = [];
        cfgtemp.channel = config{ipatient}.EMG.(markername);
        cfgtemp.latency = [-2 2];
        EMG             = ft_selectdata(cfgtemp,LFP{ipatient}{1}.(markername));
        
        if isempty(EMG.label)
            continue
        end
        
        t = EMG.time{1};
        
        %compute all envelopes, create a fieldtrip structure to average
        i=0;
        for itrial = 1 : size(EMG.trial,2)
            i = i+1;
            rect_emg                    = abs(EMG.trial{itrial}(1,:));
            [env{ipatient}.(markername).trial{i}, ~] = envelope(rect_emg, config{ipatient}.EMG.envparam, config{ipatient}.EMG.envmethod);
            env{ipatient}.(markername).time{i}       = t;
        end
        env{ipatient}.(markername).label     = {'dummy'};
        env{ipatient}.(markername).trialinfo = EMG.trialinfo;
        
        cfgtemp                   = [];
        cfgtemp.demean            = 'yes';
        cfgtemp.baselinewindow    = [-2 -1];
        env{ipatient}.(markername)= ft_preprocessing(cfgtemp, env{ipatient}.(markername));
        
        env_avg{ipatient}.(markername) = ft_timelockanalysis(cfgtemp, env{ipatient}.(markername));
        
    end
end

%% EMG : plot overdraw and avg of envelopes for each patient 
for markername = ["Myoclonus", "EMGnoEEG", "emg_all", "ied_all"]
    for ipatient = 1:size(config,2)
        
        fig = figure; hold on
        sgtitle(sprintf('%s : %d trials \n%s', config{ipatient}.prefix(1:end-1), size(env{ipatient}.(markername).trial,2), markername),'Interpreter','none','FontSize',18,'FontWeight','bold');
        
        %plot each trial with baseline substraction
        for itrial = 1:size(env{ipatient}.(markername).trial,2)
            t = env{ipatient}.(markername).time{itrial};
            trial = env{ipatient}.(markername).trial{itrial};
            p = plot(t,trial,'color', [0.5 0.25 0]);
        end
        
        %plot avg
        p = plot(t,env_avg{ipatient}.(markername).avg, 'color', [0.2 0.6 1], 'LineWidth', 2);
        
        %set figure display
        ax = axis;
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca, 'FontWeight','bold', 'Fontsize',15);
        set(gca,'TickDir','out');
        
        fname = fullfile(config{ipatient}.imagesavedir,'emg_morpho', sprintf('%semg_morphology_figure_%s',config{ipatient}.prefix, markername));
%         dtx_savefigure(fig,fname, 'png','close'); %do not save in pdf as it takes 6 hours per figure
    end
end

%% EMG : plot one example of raw EMG and its envelope