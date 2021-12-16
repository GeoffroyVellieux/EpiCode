
%% Setting parameters DTX project Paul Baudin

function [config] = lanceadams_patients_setparams(config)

disp('setting parameters for Lance-Adams patients');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis       = '/network/lustre/iss01/charpier/analyses/lgi1/lance-adams_temporaire';
    rootpath_data           = '/network/lustre/iss01/charpier/echanges/geoffroy.vellieux';
elseif ispc
    rootpath_analysis       = '\\lexport\iss01.charpier\analyses\lgi1\lance-adams_temporaire';
    rootpath_data           = '\\lexport\iss01.charpier\echanges\geoffroy.vellieux';
else
    error('Platform not supported')
end


%% Parameters common for all patients

%where to save files
config{1}.datasavedir  = fullfile(rootpath_analysis, 'data'); %.mat files
config{1}.imagesavedir = fullfile(rootpath_analysis, 'image'); %images
config{1}.tablesavedir = fullfile(rootpath_analysis, 'table'); %tables

%general infos
config{1}.name                      = {'Spike_wave','Spike','Myoclonus','Polyspike_wave', 'Polyspike', 'EMGnoEEG', 'ied_all', 'emg_all'};
config{1}.LFP.channel              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'}; %the order is important. Plot first on top.

config{1}.verifymarkers.one_off     = {'Spike_wave','Spike','Myoclonus_start','Myoclonus_end','Polyspike_wave', 'Polyspike', 'EMGnoEEG_start', 'EMGnoEEG_end'};
config{1}.verifymarkers.start       = {'Myoclonus_start','EMGnoEEG_start', 'BAD__START__'};
config{1}.verifymarkers.end         = {'Myoclonus_end'  ,'EMGnoEEG_end'  , 'BAD__END__'};

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.Spike_wave    = 'Spike_wave';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.Spike_wave      = 'Spike_wave';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.Spike_wave           = [-1, 1];
config{1}.epoch.pad.Spike_wave           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.Spike    = 'Spike';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.Spike      = 'Spike';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.Spike           = [-1, 1];
config{1}.epoch.pad.Spike           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.Myoclonus    = 'Myoclonus_start';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.Myoclonus      = 'Myoclonus_end';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.Myoclonus           = [-2, 2];
config{1}.epoch.pad.Myoclonus           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.Polyspike_wave    = 'Polyspike_wave';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.Polyspike_wave      = 'Polyspike_wave';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.Polyspike_wave           = [-1, 1];
config{1}.epoch.pad.Polyspike_wave           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.Polyspike    = 'Polyspike';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.Polyspike      = 'Polyspike';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.Polyspike           = [-1, 1];
config{1}.epoch.pad.Polyspike           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.EMGnoEEG    = 'EMGnoEEG_start';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.EMGnoEEG      = 'EMGnoEEG_end';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.EMGnoEEG           = [-2, 2];
config{1}.epoch.pad.EMGnoEEG           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.ied_all    = 'ied_all';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.ied_all      = 'ied_all';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.ied_all           = [-1, 1];
config{1}.epoch.pad.ied_all           = 0;

%cut into trials according to Muse markers with readLFP.m
config{1}.muse.startmarker.emg_all    = 'emg_all_start';   % start and end Muse marker. For defining trials
config{1}.muse.endmarker.emg_all      = 'emg_all_end';   % start and end Muse marker. For defining trials
config{1}.epoch.toi.emg_all           = [-2, 2];
config{1}.epoch.pad.emg_all           = 0;

%preprocessing of the LFP and EMG in readLFP.m, and som einformations required for
%plots of the LFP and EMG
config{1}.LFP.name                  = config{1}.name;
config{1}.LFP.write                 = true;

config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 0;
config{1}.LFP.resamplefs            = 256;
config{1}.LFP.baseline              = 'no';%No baseline when reading LFP. Baseline applied for each trial after trial definition
config{1}.LFP.reref                 = 'yes';
config{1}.LFP.rerefmethod           = 'avg';
config{1}.LFP.refchannel            = config{1}.LFP.channel;
config{1}.LFP.bsfilter              = 'yes';
config{1}.LFP.bsfreq                = [49 51];
config{1}.LFP.lpfilter              = 'no';
config{1}.LFP.lpfreq                = 30;
config{1}.LFP.lpfilttype            = 'fir';

config{1}.EMG.hpfilter              = 'yes';
config{1}.EMG.hpfreq                = 10;
config{1}.EMG.bsfilter              = 'yes';
config{1}.EMG.bsfreq                = [49 51];
config{1}.EMG.bsfiltord             = 3;
config{1}.EMG.reref                 = 'no';
config{1}.EMG.rerefmethod           = [];
config{1}.EMG.refchannel            = [];

config{1}.cluster.name              = {'ied_all'};
config{1}.cluster.channel           = {'Fpz', 'Fz', 'Cz', 'Pz', 'Oz'}; % CHECK REFERENCE
config{1}.cluster.demean            = 'yes';
config{1}.cluster.reref             = 'yes';
config{1}.cluster.refmethod         = 'avg';
config{1}.cluster.refchannel        = config{1}.LFP.channel;
config{1}.cluster.resamplefs        = 256;
config{1}.cluster.baselinewindow    = [-0.7 -0.3];
config{1}.cluster.latency           = [-0.1 0.2];
config{1}.cluster.dbscan            = 'no';
config{1}.cluster.kmeans            = 'yes';
config{1}.cluster.kmedoids          = 'yes';
config{1}.cluster.N                 = 4; 
config{1}.cluster.align.latency     = [-0.3, 0.5];

config{1}.align.name                = {'Spike_wave','Spike','Polyspike_wave', 'Polyspike', 'ied_all'};
config{1}.align.channel             = {'Cz'};
config{1}.align.removeartefacts     = true;
config{1}.align.reref               = 'yes';
config{1}.align.refmethod           = 'avg';
config{1}.align.refchannel          = config{1}.LFP.channel;
config{1}.align.lpfilter            = 'yes';
config{1}.align.lpfreq              = 30;
config{1}.align.latency.Spike_wave          = [-0.3 0.5];
config{1}.align.latency.Spike               = [-0.3 0.5];
config{1}.align.latency.Polyspike_wave      = [-0.3 0.5];
config{1}.align.latency.Polyspike           = [-0.3 0.5];
config{1}.align.latency.ied_all             = [-0.3 0.5];
config{1}.align.reject              = 'BAD_cnt';

config{1}.EMG.envmethod             = 'rms';
config{1}.EMG.envparam              = 50;

%% Patient 1 

%general infos
config{1}                           = config{1};
config{1}.prefix                    = 'pat_12605-';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_12605');
config{1}.directorylist{1}          = {'EEG_120806'}; %dir = eeg file with all the electrodess

for ifield = string(config{1}.name)
    %electrode used for peak alignment and morphology measurements
%     config{1}.align.channel.(ifield)  = 'Cz';
%     config{1}.align.channel.(ifield)  = 'Cz';
    config{1}.EMG.(ifield)            = 'EMG1+';
end

config{1}.topoplot.channel = 'Cz';
config{1}.topoplot.toi = [-0.05 0.15];

end


