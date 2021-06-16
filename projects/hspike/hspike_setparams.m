function [config] = hspike_setparams(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [config] = hspike_setparams(varargin)
%
% This function outputs all the settings of a study, to be defined below
%
% Note the consideration of the operating system, because the pointers to
% the file server is dealt with differently. This could be different for
% you. You can force an OS with the input argument: 'pc' or 'unix'.
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('setting parameters');

pc = ispc;

% overwrite if os is forced
if nargin == 1
    if strcmp(varargin{1}, 'pc')
        pc = true;
    elseif strcmp(varargin{1}, 'unix')
        pc = false;
    else
        error('os not recognized');
    end
end

if pc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
    rootpath_epiloc   	= '\\lexport\iss01.epimicro\patients\epiloc\analyse';       
else
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
    rootpath_epiloc     = '/network/lustre/iss01/epimicro/patients/epiloc/analyse';
end


%% Patient 1

config{1}.prefix                        = '2711-';
config{1}.rawdir                        = fullfile(rootpath_data,     'pat_02711_1193', 'eeg');
config{1}.datasavedir                   = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.imagesavedir                  = fullfile(rootpath_analysis, 'images', 'hspike');
config{1}.locfile                       = fullfile(rootpath_epiloc,   'pat_02711_1193', 'pat_02711_1193_anatomical_localizations.csv');
config{1}.visible                       = 'on';

config{1}.name                          = {'Hspike'};
config{1}.muse.startmarker.Hspike       = "Hspike";
config{1}.muse.endmarker.Hspike         = "Hspike";
config{1}.muse.startmarker.template1    = "template1";
config{1}.muse.endmarker.template1      = "template1";
config{1}.muse.startmarker.template2    = "template2";
config{1}.muse.endmarker.template2      = "template2";
config{1}.muse.startmarker.template3    = "template3";
config{1}.muse.endmarker.template3      = "template3";
config{1}.muse.startmarker.template4    = "template4";
config{1}.muse.endmarker.template4      = "template4";
config{1}.muse.startmarker.template5    = "template5";
config{1}.muse.endmarker.template5      = "template5";
config{1}.muse.startmarker.template6    = "template6";
config{1}.muse.endmarker.template6      = "template6";
config{1}.muse.startmarker.combined1    = "combined1";
config{1}.muse.endmarker.combined1      = "combined1";
config{1}.muse.startmarker.combined2    = "combined2";
config{1}.muse.endmarker.combined2      = "combined2";
config{1}.muse.startmarker.combined3    = "combined3";
config{1}.muse.startmarker.combined3    = "combined3";
config{1}.muse.backupdir                = fullfile(rootpath_analysis, 'markerbackup');

config{1}.hyp.imagesavedir          = fullfile(rootpath_analysis, 'images', 'hspike');
config{1}.hyp.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{1}.hyp.markerdir             = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.hyp.micromedchannel       = 'F3p6';
config{1}.hyp.markers               = {'Hspike','template1','template2','template3','template4','template5','template6','combined1','combined2','combined3'};
config{1}.hyp.overwrite             = 'append';
config{1}.hyp.spikewindow           = 60;
config{1}.hyp.spikewindowoverlap    = 0.5;

config{1}.epoch.toi.Hspike          = [-0.5  1];
config{1}.epoch.pad.Hspike          = 1;
config{1}.epoch.toi.combined1       = [-0.5  1];
config{1}.epoch.pad.combined1       = 1;
config{1}.epoch.toi.combined2       = [-0.5  1];
config{1}.epoch.pad.combined2       = 1;
config{1}.epoch.toi.combined3       = [-0.5  1];
config{1}.epoch.pad.combined3       = 1;
config{1}.epoch.toi.template1       = [-0.5  1];
config{1}.epoch.pad.template1       = 1;
config{1}.epoch.toi.template2       = [-0.5  1];
config{1}.epoch.pad.template2       = 1;
config{1}.epoch.toi.template3       = [-0.5  1];
config{1}.epoch.pad.template3       = 1;
config{1}.epoch.toi.template4       = [-0.5  1];
config{1}.epoch.pad.template4       = 1;
config{1}.epoch.toi.template5       = [-0.5  1];
config{1}.epoch.pad.template5       = 1;
config{1}.epoch.toi.template6       = [-0.5  1];
config{1}.epoch.pad.template6       = 1;

config{1}.LFP.name                      = {'Hspike'};
config{1}.LFP.channel                   = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.LFP.hpfilter                  = 'no';
config{1}.LFP.hpfreq                    = 1;
config{1}.LFP.resamplefs                = 250;
config{1}.LFP.baseline                  = 'yes';
config{1}.LFP.baselinewindow.combined1  = [-1, -0.5];
config{1}.LFP.baselinewindow.combined2  = [-1, -0.5];
config{1}.LFP.baselinewindow.combined3  = [-1, -0.5];
config{1}.LFP.baselinewindow.Hspike     = [-1, -0.5];
config{1}.LFP.baselinewindow.template1  = [-1, -0.5];
config{1}.LFP.baselinewindow.template2  = [-1, -0.5];
config{1}.LFP.baselinewindow.template3  = [-1, -0.5];
config{1}.LFP.baselinewindow.template4  = [-1, -0.5];
config{1}.LFP.baselinewindow.template5  = [-1, -0.5];
config{1}.LFP.baselinewindow.template6  = [-1, -0.5];

config{1}.TFR.name                      = {'Hspike'};
config{1}.TFR.channel                   = 'all';
config{1}.TFR.keeptrials                = 'no';
config{1}.TFR.foi.combined1             = 1:1:100;
config{1}.TFR.foi.combined2             = 1:1:100;
config{1}.TFR.foi.combined3             = 1:1:100;
config{1}.TFR.foi.Hspike                = 1:1:100;
config{1}.TFR.foi.template1             = 1:1:100;
config{1}.TFR.foi.template2             = 1:1:100;
config{1}.TFR.foi.template3             = 1:1:100;
config{1}.TFR.foi.template4             = 1:1:100;
config{1}.TFR.foi.template5             = 1:1:100;
config{1}.TFR.foi.template6             = 1:1:100;

config{1}.TFR.t_ftimwin.combined1    = 5./config{1}.TFR.foi.combined1;
config{1}.TFR.t_ftimwin.combined2    = 5./config{1}.TFR.foi.combined2;
config{1}.TFR.t_ftimwin.combined3    = 5./config{1}.TFR.foi.combined3;
config{1}.TFR.t_ftimwin.Hspike       = 5./config{1}.TFR.foi.Hspike;
config{1}.TFR.t_ftimwin.template1    = 5./config{1}.TFR.foi.combined1;
config{1}.TFR.t_ftimwin.template2    = 5./config{1}.TFR.foi.combined2;
config{1}.TFR.t_ftimwin.template3    = 5./config{1}.TFR.foi.combined3;
config{1}.TFR.t_ftimwin.template4    = 5./config{1}.TFR.foi.combined1;
config{1}.TFR.t_ftimwin.template5    = 5./config{1}.TFR.foi.combined2;
config{1}.TFR.t_ftimwin.template6    = 5./config{1}.TFR.foi.combined3;

config{1}.TFR.toi.combined1          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.combined2          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.combined3          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.Hspike             = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template1          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template2          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template3          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template4          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template5          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);
config{1}.TFR.toi.template6          = config{1}.epoch.toi.Hspike(1) : 0.005 : config{1}.epoch.toi.Hspike(2);

config{1}.TFR.bl.combined1           = [-0.5, -0.2];
config{1}.TFR.bl.combined2           = [-0.5, -0.2];
config{1}.TFR.bl.combined3           = [-0.5, -0.2];
config{1}.TFR.bl.Hspike              = [-0.5, -0.2];
config{1}.TFR.bl.template1           = [-0.5, -0.2];
config{1}.TFR.bl.template2           = [-0.5, -0.2];
config{1}.TFR.bl.template3           = [-0.5, -0.2];
config{1}.TFR.bl.template4           = [-0.5, -0.2];
config{1}.TFR.bl.template5           = [-0.5, -0.2];
config{1}.TFR.bl.template6           = [-0.5, -0.2];

config{1}.align.name                = {'Hspike'};
config{1}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.align.reref               = 'yes';
config{1}.align.refmethod           = 'bipolar';
config{1}.align.latency.Hspike      = [-0.1 0.4];
config{1}.align.zerochannel         = 'HaT2_1-HaT2_2';

config{1}.cluster.name              = {'Hspike'};
config{1}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'}; % CHECK REFERENCE
config{1}.cluster.demean            = 'yes';
config{1}.cluster.reref             = 'yes';
config{1}.cluster.refmethod         = 'bipolar';
config{1}.cluster.resamplefs        = 250;
config{1}.cluster.baselinewindow    = [-0.5 -0.2];
config{1}.cluster.latency           = [-0.5 1];
config{1}.cluster.dbscan            = 'no';
config{1}.cluster.kmeans            = 'no';
config{1}.cluster.kmedoids          = 'yes';
config{1}.cluster.N                 = 6; % can be list, e.g. p2 3 4 6]
config{1}.cluster.align.latency     = [-0.2, 0.3];

config{1}.template.reref            = 'yes';
config{1}.template.refmethod        = 'bipolar';
config{1}.template.latency          = [-0.2, 0.5];
config{1}.template.resamplefs       = 250;
config{1}.template.threshold        = 2.7;

config{1}.circus.correct_chunk{1}               = true;
config{1}.circus.correct_chunk{2}               = true;
config{1}.circus.correct_chunk{3}               = true;
config{1}.circus.channel                        = {'mHaT2_1', 'mHaT2_3', 'mHaT2_4','mHaT2_6', 'mHaT2_7', 'mHaT2_8'};
config{1}.circus.reref                          = 'no';
config{1}.circus.refchan                        = '';
config{1}.circus.outputdir                      = 'SpykingCircus_new';
config{1}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'hspike', 'SpykingCircus.params');
config{1}.circus.params.detection.spike_thresh  = '6';
config{1}.circus.params.filtering.cut_off       = '300, auto';
config{1}.circus.params.filtering.remove_median = 'False';
config{1}.circus.params.clustering.max_elts     = '20000';
config{1}.circus.params.detection.peaks         = 'negative';
config{1}.circus.params.data.stream_mode        = 'mapping-file';
config{1}.circus.params.data.mapping_file       = 'filelist.txt';

config{1}.spike.name                = {'combined1', 'combined2'};
config{1}.spike.overlap             = [];                                    
config{1}.spike.slidestep           = [0.01, 0.01, 0.001];
config{1}.spike.toi.combined1       = [-0.5, 1];          
config{1}.spike.toi.combined2       = [-0.5, 1];           
config{1}.spike.toi.combined3       = [-0.5, 1];           
config{1}.spike.toi.Hspike          = [-0.5, 1];    
config{1}.spike.toi.template1       = [-0.5, 1];          
config{1}.spike.toi.template2       = [-0.5, 1];           
config{1}.spike.toi.template3       = [-0.5, 1];    
config{1}.spike.toi.template4       = [-0.5, 1];          
config{1}.spike.toi.template5       = [-0.5, 1];           
config{1}.spike.toi.template6       = [-0.5, 1];    

config{1}.spike.bl.combined1        = [-0.5, -0.2];
config{1}.spike.bl.combined2        = [-0.5, -0.2];
config{1}.spike.bl.combined3        = [-0.5, -0.2];
config{1}.spike.bl.Hspike           = [-0.5, -0.2];
config{1}.spike.bl.template1        = [-0.5, -0.2];
config{1}.spike.bl.template2        = [-0.5, -0.2];
config{1}.spike.bl.template3        = [-0.5, -0.2];
config{1}.spike.bl.template4        = [-0.5, -0.2];
config{1}.spike.bl.template5        = [-0.5, -0.2];
config{1}.spike.bl.template6        = [-0.5, -0.2];

config{1}.spike.pad.combined1        = 0.1;
config{1}.spike.pad.combined2        = 0.1;
config{1}.spike.pad.combined3        = 0.1;
config{1}.spike.pad.Hspike           = 0.1;
config{1}.spike.pad.template1        = 0.1;
config{1}.spike.pad.template2        = 0.1;
config{1}.spike.pad.template3        = 0.1;
config{1}.spike.pad.template4        = 0.1;
config{1}.spike.pad.template5        = 0.1;
config{1}.spike.pad.template6        = 0.1;

config{1}.spike.resamplefs.combined1 = 1000;
config{1}.spike.resamplefs.combined2 = 1000;
config{1}.spike.resamplefs.combined3 = 1000;
config{1}.spike.resamplefs.Hspike    = 1000;
config{1}.spike.resamplefs.template1 = 1000;
config{1}.spike.resamplefs.template2 = 1000;
config{1}.spike.resamplefs.template3 = 1000;
config{1}.spike.resamplefs.template4 = 1000;
config{1}.spike.resamplefs.template5 = 1000;
config{1}.spike.resamplefs.template6 = 1000;

config{1}.spike.pre                 = 0.001;
config{1}.spike.post                = 0.002;
config{1}.spike.baseline            = [-0.001 -0.0005];
config{1}.spike.ISIbins             = 0 : 0.0005 : 0.150; %in s
config{1}.spike.nrsdfbins           = 100;
config{1}.spike.psthbin.combined1   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.combined2   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.combined3   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.Hspike      = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template1   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template2   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template3   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template4   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template5   = 0.01; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.template6   = 0.01; % depends a lot on pattern, default is too large

config{1}.spike.sdftimwin.combined1 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.combined2 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.combined3 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.Hspike    = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template1 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template2 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template3 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template4 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template5 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.template6 = [-0.01 0.01]; % depends a lot on pattern, default is 0.05

config{1}.spikewin.windowsize       = 10; % to prevent too much data lost on overlaps with IEDs
config{1}.spikewin.windowoverlap    = 0.5;

config{1}.stats.toi.combined1       = [-0.5, 1];
config{1}.stats.toi.combined2       = [-0.5, 1];
config{1}.stats.toi.combined3       = [-0.5, 1];
config{1}.stats.toi.Hspike          = [-0.5, 1];
config{1}.stats.toi.template1       = [-0.5, 1];
config{1}.stats.toi.template2       = [-0.5, 1];
config{1}.stats.toi.template3       = [-0.5, 1];
config{1}.stats.toi.template4       = [-0.5, 1];
config{1}.stats.toi.template5       = [-0.5, 1];
config{1}.stats.toi.template6       = [-0.5, 1];

config{1}.stats.bl.combined1        = [-0.5 -0.2];
config{1}.stats.bl.combined2        = [-0.5 -0.2];
config{1}.stats.bl.combined3        = [-0.5 -0.2];
config{1}.stats.bl.Hspike           = [-0.5 -0.2];
config{1}.stats.bl.template1        = [-0.5 -0.2];
config{1}.stats.bl.template2        = [-0.5 -0.2];
config{1}.stats.bl.template3        = [-0.5 -0.2];
config{1}.stats.bl.template4        = [-0.5 -0.2];
config{1}.stats.bl.template5        = [-0.5 -0.2];
config{1}.stats.bl.template6        = [-0.5 -0.2];
config{1}.stats.alpha               = 0.025;

config{1}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined2';
                                        'template3', 'combined1';
                                        'template5', 'combined2';
                                        'template6', 'combined2'};
             
                                    

config{1}.plot.reref                = 'yes';
config{1}.plot.refmethod            = 'bipolar';                                  
config{1}.plot.name                 = {'combined1', 'combined2'};   

config{1}.plot.unit{1}              = [-1,  -1,  -1,  -1,  -1];                     
config{1}.plot.trial.template1{1}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template2{1}   = [1,    2,    3,    4,    5];
config{1}.plot.trial.template3{1}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template4{1}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template5{1}   = [10,    20,   30,   40,   50];
config{1}.plot.trial.template6{1}   = [100,  120,  130,  140,  150];

config{1}.plot.unit{2}              = [-1,  -1,  -1,  -1,  -1];                     
config{1}.plot.trial.template1{2}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template2{2}   = [1,    2,    3,    4,    5];
config{1}.plot.trial.template3{2}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template4{2}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template5{2}   = [10,    20,   30,   40,   50];
config{1}.plot.trial.template6{2}   = [100,  120,  130,  140,  150];

config{1}.plot.unit{3}              = [-1,  -1,  -1,  -1,  -1];                  
config{1}.plot.trial.template1{3}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template2{3}   = [1,    2,    3,    4,    5];
config{1}.plot.trial.template3{3}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template4{3}   = [100,  120,  130,  140,  150];
config{1}.plot.trial.template5{3}   = [10,    20,   30,   40,   50];
config{1}.plot.trial.template6{3}   = [100,  120,  130,  140,  150];


%% Patient 2
% Patterns more clear in bipolar reference. Consider this for template
% matching.
config{2}                           = config{1};
config{2}.prefix                    = '2718-';
config{2}.rawdir                    = fullfile(rootpath_data, 'pat_02718_1201', 'eeg');
config{2}.locfile                   = fullfile(rootpath_epiloc, 'pat_02718_1201', 'pat_02718_1201_anatomical_localizations.csv');
config{2}.hyp.micromedchannel       = 'HaT11';
config{2}.LFP.channel               = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{2}.align.channel             = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{2}.align.latency             = [-0.1, 0.4];
config{2}.align.zerochannel         = 'HmT2_1-HmT2_2';
config{2}.cluster.channel           = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{2}.cluster.align.latency     = [-0.1, 0.4];
config{2}.cluster.reref             = 'yes';
config{2}.cluster.refmethod         = 'bipolar';
config{2}.template.threshold        = 3;
config{2}.template.reref            = 'yes';
config{2}.template.refmethod        = 'bipolar';

config{2}.circus.correct_chunk{1}      = false;
config{2}.circus.correct_chunk{2}      = true;
config{2}.circus.correct_chunk{3}      = true;
config{2}.circus.channel            = {'mHaT1_7'};
config{2}.circus.reref              = 'no';
config{2}.circus.params.filtering.remove_median = 'False';
config{2}.circus.params.detection.spike_thresh  = '6';
config{2}.circus.params.detection.peaks         = 'negative';
config{2}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined1';
                                        'template4', 'combined1';
                                        'template6', 'combined1'};
config{2}.spike.name                = {'combined1'};
% config{2}.spike.name                = {'Hspike'};

config{2}.plot                      = [];
config{2}.plot.unit{1}              = [-1,  -1,  -1,  -1,  -1];                     
config{2}.plot.trial.template1{1}   = [100,  120,  130,  140,  150];
config{2}.plot.trial.template2{1}   = [1,    2,    3,    4,    5];
config{2}.plot.trial.template3{1}   = [11,   12,    13,    14,  15];
config{2}.plot.trial.template4{1}   = [100,  120,  130,  140,  150];
config{2}.plot.trial.template5{1}   = [10,    20,   30,   40,   50];
config{2}.plot.trial.template6{1}   = [100,  120,  130,  140,  150];

config{2}.plot.unit{2}              = [-1,  -1,  -1,  -1,  -1];                     
config{2}.plot.trial.template1{2}   = [100, 120,  130,  140,  150];
config{2}.plot.trial.template2{2}   = [1,    2,    3,    4,    5];
config{2}.plot.trial.template3{2}   = [1,    1,    1,    1,    1];
config{2}.plot.trial.template4{2}   = [100, 120,  130,  140,  150];
config{2}.plot.trial.template5{2}   = [10,    20,   30,   40,   50];
config{2}.plot.trial.template6{2}   = [100,  120,  130,  140,  150];

config{2}.plot.unit{3}              = [-1,  -1,  -1,  -1,  -1];                  
config{2}.plot.trial.template1{3}   = [100,  120,  130,  140,  150];
config{2}.plot.trial.template2{3}   = [1,    2,    3,    4,    5];
config{2}.plot.trial.template3{3}   = [11,   12,    13,    14,  15];
config{2}.plot.trial.template4{3}   = [1,    2,    3,    4,    5];
config{2}.plot.trial.template5{3}   = [10,    20,   30,   40,   50];
config{2}.plot.trial.template6{3}   = [100,  120,  130,  140,  150];

%% Patient 3
config{3}                           = config{1};
config{3}.prefix                    = '2660-';
config{3}.rawdir                    = fullfile(rootpath_data, 'pat_02660_1136', 'eeg');
config{3}.locfile                   = fullfile(rootpath_epiloc, 'pat_02660_1136', 'pat_02660_1136_anatomical_localizations.csv');
config{3}.hyp.micromedchannel       = 'Ha2g1';
config{3}.LFP.channel               = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Ha2g_5','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4','_Hm2g_5'};
config{3}.align.channel             = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Ha2g_5','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4','_Hm2g_5'};
config{3}.align.latency             = [-0.1 0.1];
config{3}.align.zerochannel         = 'Hm2g_2';
config{3}.cluster.channel           = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Ha2g_5','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4','_Hm2g_5'};
config{3}.cluster.align.latency     = [-0.2 0.5];
config{3}.cluster.reref             = 'no';
config{3}.cluster.refmethod         = 'bipolar';
config{3}.template.latency          = [-0.2, 0.5];
config{3}.template.threshold        = 3;
config{3}.template.reref            = 'no';
config{3}.template.refmethod        = 'bipolar';
config{3}.circus.params.filtering.remove_median = 'False';
config{3}.circus.channel            = {'mHa2g_3','mHa2g_4','mHa2g_7','mHa2g_8','mTBmd_1','mTBmd_2','mTBmd_4','mTBmd_5','mTBmd_6','mTBmd_7','mTBmd_8'}; % Changes over night! Night 1 would have some on ,'mHa2g_2 (current ref) when rereferencing to another; In night 2 ref changes to mHa2g_8; night 3 ref chanes to mHa3g_7
config{3}.circus.channelname        = {'mHa2g',  'mHa2g',  'mHa2g',  'mHa2g',  'mTBmd',  'mTBmd',  'mTBmd',  'mTBmd',  'mTBmd',  'mTBmd',  'mTBmd', };

config{3}.circus.params.detection.spike_thresh  = '7'; % TOO MUCH 50 HZ contamination -> TRY INCREASING THERSHOLD
config{3}.circus.params.filtering.remove_median = 'True';
% config{3}.circus.version                        = 'fieldtrip';

config{3}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined1';
                                        'template3', 'combined2';
                                        'template4', 'combined1';
                                        'template5', 'combined1';
                                        'template6', 'combined1'};
config{3}.spike.name                = { 'combined1', 'combined2'};
                                    
%% Patient 4
config{4}                           = config{1};
config{4}.prefix                    = '2680-';
config{4}.rawdir                    = fullfile(rootpath_data, 'pat_02680_1158', 'eeg');
config{4}.locfile                   = fullfile(rootpath_epiloc, 'pat_02680_1158', 'pat_02680_1158_anatomical_localizations.csv');
config{4}.hyp.micromedchannel       = 'HaT21';
config{4}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{4}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{4}.align.latency             = [-0.1 0.1];
config{4}.align.zerochannel         = 'HmT2_1';
config{4}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5'};
config{4}.cluster.align.latency     = [-0.2 0.5];
config{4}.cluster.reref             = 'no';
config{4}.cluster.refmethod         = 'bipolar';
config{4}.template.threshold        = 2.7;
config{4}.template.reref            = 'no';
config{4}.template.refmethod        = 'bipolar';
config{4}.circus.channel            = {'mHaT2_3','mHaT2_5','mHaT2_7'}; % might just as well remove amygdata electrodes as they are too noisy
config{4}.circus.params.detection.spike_thresh  = '9';
config{4}.circus.params.filtering.remove_median = 'False';
config{4}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined2';
                                        'template3', 'combined2';
                                        'template4', 'combined2';
                                        'template5', 'combined2';
                                        'template6', 'combined3'};
config{4}.spike.name                =  {'combined1', 'combined2', 'combined3'};

%% Patient 5
% some clear cases for merging in mLMI1 & mLMS2
config{5}                           = config{1};
config{5}.prefix                    = '2689-';
config{5}.rawdir                    = fullfile(rootpath_data, 'pat_02689_1168', 'eeg');
config{5}.locfile                   = fullfile(rootpath_epiloc, 'pat_02689_1168', 'pat_02689_1168_anatomical_localizations.csv');
config{5}.hyp.micromedchannel       = 'HmT31';
config{5}.LFP.channel               = {'_HmT3_1','_HmT3_2','_HmT3_3','_HmT3_4','_HmT3_5','_HpNI_1','_HpNI_2','_HpNI_3','_HpNI_4','_HpNI_5'};
config{5}.align.channel             = {'_HmT3_1','_HmT3_2','_HmT3_3','_HmT3_4','_HmT3_5','_HpNI_1','_HpNI_2','_HpNI_3','_HpNI_4','_HpNI_5'};
config{5}.align.zerochannel         = 'HmT3_2-HmT3_3';
config{5}.cluster.channel           = {'_HmT3_1','_HmT3_2','_HmT3_3','_HmT3_4','_HmT3_5','_HpNI_1','_HpNI_2','_HpNI_3','_HpNI_4','_HpNI_5'};
config{5}.cluster.align.latency     = [-0.2 0.5];
config{5}.cluster.reref             = 'yes';
config{5}.cluster.refmethod         = 'bipolar';
config{5}.template.threshold        = 3;
config{5}.template.reref            = 'yes';
config{5}.template.refmethod        = 'bipolar';
config{5}.circus.channel            = {'mHmT3_4','mHmT3_5','mHmT3_7', 'mLMI1_2','mLMI1_3','mLMI1_4','mLMI1_7', 'mLMS2_2'};
config{5}.circus.channelname        = {'mHmT3',  'mHmT3',  'mHmT3',   'mLMI1',  'mLMI1',  'mLMI1',  'mLMI1',   'mLMS2'};
config{5}.circus.params.detection.spike_thresh  = '6';
config{5}.circus.params.filtering.remove_median = 'False';
config{5}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined2';
                                        'template3', 'combined3';
                                        'template5', 'combined2'; };
config{5}.spike.name               =  {'combined1', 'combined2', 'combined3'};
                              
%% Patient 6
config{6}                           = config{1};
config{6}.prefix                    = '2651-'; %% REDO SPYKING CIRCUS AFTER REDOING ARTEFACTS
config{6}.rawdir                    = fullfile(rootpath_data, 'pat_02651_1127', 'eeg');
config{6}.locfile                   = fullfile(rootpath_epiloc, 'pat_02651_1127', 'pat_02651_1127_anatomical_localizations.csv');
config{6}.hyp.micromedchannel       = 'HmT21';
config{6}.LFP.channel               = {'_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5','_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{6}.align.channel             = {'_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5','_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{6}.align.zerochannel         = 'HmT2_1';
config{6}.cluster.channel           = {'_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HmT2_5','_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{6}.cluster.reref             = 'no';
config{6}.cluster.refmethod         = 'bipolar';
config{6}.template.threshold        = 3.2;
config{6}.template.threshold        = 3;
config{6}.template.reref            = 'no';
config{6}.template.refmethod        = 'bipolar';
config{6}.circus.channel            = {'mHaT2_6','mHaT2_7','mHaT2_8','mHmT2_3','mHmT2_4','mHmT2_5','mHmT2_6'};
config{6}.circus.channelname        = {'mHaT2','mHaT2','mHaT2','mHmT2','mHmT2','mHmT2','mHmT2'};
config{6}.circus.params.detection.spike_thresh  = '6';
config{6}.circus.params.filtering.remove_median = 'False';
config{6}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined1';
                                        'template3', 'combined1';
                                        'template4', 'combined1';
                                        'template5', 'combined1';
                                        'template6', 'combined1'};
config{6}.spike.name               =  {'combined1'};
                                  
%% Patient 7
config{7}                           = config{1};
config{7}.prefix                    = '2619-';
config{7}.rawdir                    = fullfile(rootpath_data, 'pat_02619_1078', 'eeg');
config{7}.locfile                   = fullfile(rootpath_epiloc, 'pat_02619_1078', 'pat_02619_1078_anatomical_localizations.csv');
config{7}.hyp.micromedchannel       = 'HaTB1';
config{7}.LFP.channel               = {'_HaTB_1','_HaTB_2','_HaTB_3','_HaTB_4','_HaTB_5'};
config{7}.align.channel             = {'_HaTB_1','_HaTB_2','_HaTB_3','_HaTB_4','_HaTB_5'};
config{7}.align.zerochannel         = 'HaTB_1';
config{7}.cluster.channel           = {'_HaTB_1','_HaTB_2','_HaTB_3','_HaTB_4','_HaTB_5'};
config{7}.cluster.reref             = 'no';
config{7}.cluster.refmethod         = 'bipolar';
config{7}.template.threshold        = 5;
config{7}.template.reref            = 'no';
config{7}.template.refmethod        = 'bipolar';
config{7}.circus.channel            = {'mAmT2_1','mAmT2_2','mAmT2_4','mAmT2_5','mAmT2_6'};
config{7}.circus.reref              = 'no';
config{7}.circus.params.detection.spike_thresh  = '7';
config{7}.circus.params.filtering.remove_median = 'False';
config{7}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined1';
                                        'template3', 'combined2';
                                        'template4', 'combined2';
                                        'template5', 'combined1';
                                        'template6', 'combined2'};
config{7}.spike.name                =  {'combined1', 'combined2'};

%% Patient NO MUA ON FIRST DAY - ONLY HYPNOGRAM ON FIRST DAY
config{8}                           = config{1};
config{8}.prefix                    = '2578-';
config{8}.rawdir                    = fullfile(rootpath_data, 'pat_02578_1036', 'eeg');
config{8}.epilocdir                 = fullfile(rootpath_epiloc, 'pat_02578_1036', 'pat_02578_1036_anatomical_localizations.csv');
config{8}.hyp.micromedchannel       = 'HaTB1';
config{8}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT1_1','_HmT1_2','_HmT1_3','_HmT1_4'};
config{8}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4'};
config{8}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT1_1','_HmT1_2','_HmT1_3','_HmT1_4'};
config{8}.cluster.reref             = 'no';
config{8}.cluster.refmethod         = 'bipolar';
config{8}.template.threshold        = 3;
config{8}.template.reref            = 'no';
config{8}.template.refmethod        = 'bipolar';
config{8}.circus.channel            = {'mAmT2_3'};
config{8}.circus.reref              = 'no';
config{8}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined2';
                                        'template4', 'combined2';
                                        'template5', 'combined1'};

config{9}                           = config{1};
config{9}.prefix                    = '2614-';
config{9}.rawdir                    = fullfile(rootpath_data, 'pat_02578_1036', 'eeg');
config{9}.epilocdir                 = fullfile(rootpath_epiloc, 'pat_02578_1036', 'pat_02578_1036_anatomical_localizations.csv');
config{9}.hyp.micromedchannel       = 'HaTB1';
config{9}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT1_1','_HmT1_2','_HmT1_3','_HmT1_4'};
config{9}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4'};
config{9}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT1_1','_HmT1_2','_HmT1_3','_HmT1_4'};
config{9}.cluster.reref             = 'no';
config{9}.cluster.refmethod         = 'bipolar';
config{9}.template.threshold        = 3;
config{9}.template.reref            = 'no';
config{9}.template.refmethod        = 'bipolar';
config{9}.circus.channel            = {'mAmT2_3'};
config{9}.circus.reref              = 'no';
config{9}.editmarkerfile.torename   = { 'template1', 'combined1';
                                        'template2', 'combined2';
                                        'template4', 'combined2';
                                        'template5', 'combined1'};

                                                                      
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
% %% Patient NO HYPNOGRAM
% config{}                           = config{1};
% config{}.prefix                    = '2599-';
% config{}.rawdir                    = fullfile(rootpath_data,'pat_02599_1057', 'eeg');
% config{}.hyp.micromedchannel       = 'HaT2';
% config{}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HpT2_1','_HpT2_2','_HpT2_3','_HpT2_4'};
% config{}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HpT2_1','_HpT2_2','_HpT2_3','_HpT2_4'};
% config{}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HmT2_1','_HmT2_2','_HmT2_3','_HmT2_4','_HpT2_1','_HpT2_2','_HpT2_3','_HpT2_4'};
% config{}.template.threshold        = 6;
% config{}.circus.channel            = {'mHaT2_6'};
% config{}.circus.reref              = 'no';
% config{}.cluster.reref             = 'no';
% config{}.cluster.refmethod         = 'bipolar';
% config{}.template.threshold        = 2.8;
% config{}.template.reref            = 'no';
% config{}.template.refmethod        = 'bipolar';

%% DATA

config{1}.directorylist             = [];
config{1}.directorylist{1}          =  {'02711_2019-04-17_12-29',...
                                        '02711_2019-04-17_14-29',...
                                        '02711_2019-04-17_16-29',...
                                        '02711_2019-04-17_18-29',...
                                        '02711_2019-04-17_20-29',...
                                        '02711_2019-04-17_22-29',...
                                        '02711_2019-04-18_00-29',...
                                        '02711_2019-04-18_02-29',...
                                        '02711_2019-04-18_04-29',...
                                        '02711_2019-04-18_06-29',...
                                        '02711_2019-04-18_08-29',...
                                        '02711_2019-04-18_10-29',...
                                        '02711_2019-04-18_11-04'};      % TIMECODES FIXED
config{1}.directorylist{2}          =  {'02711_2019-04-18_13-04',...
                                        '02711_2019-04-18_15-04',...
                                        '02711_2019-04-18_17-04',...
                                        '02711_2019-04-18_19-04',...
                                        '02711_2019-04-18_21-04',...
                                        '02711_2019-04-18_23-04',...
                                        '02711_2019-04-19_01-04',...    % TIMECODES FIXED
                                        '02711_2019-04-19_03-04',...
                                        '02711_2019-04-19_05-04',...
                                        '02711_2019-04-19_07-04',...
                                        '02711_2019-04-19_09-04',...
                                        '02711_2019-04-19_10-00',...
                                        '02711_2019-04-19_12-00'};
config{1}.directorylist{3}          =  {'02711_2019-04-19_14-00',...
                                        '02711_2019-04-19_16-00',...
                                        '02711_2019-04-19_18-00',...
                                        '02711_2019-04-19_20-00',...
                                        '02711_2019-04-19_22-00',...
                                        '02711_2019-04-20_00-00',...
                                        '02711_2019-04-20_02-00',...
                                        '02711_2019-04-20_04-00',...
                                        '02711_2019-04-20_06-00',...
                                        '02711_2019-04-20_08-00',...
                                        '02711_2019-04-20_10-00',...
                                        '02711_2019-04-20_12-00'};

config{2}.directorylist             = [];
config{2}.directorylist{1}          =  {'02718_2019-05-14_12-31',...
                                        '02718_2019-05-14_14-31',...
                                        '02718_2019-05-14_16-31',...
                                        '02718_2019-05-14_18-31',...
                                        '02718_2019-05-14_20-31',... % very artefacted
                                        '02718_2019-05-14_22-31',... % MISSING DATA
                                        '02718_2019-05-15_00-31',...
                                        '02718_2019-05-15_02-31',...
                                        '02718_2019-05-15_04-31',...
                                        '02718_2019-05-15_06-31',...
                                        '02718_2019-05-15_08-31',...
                                        '02718_2019-05-15_10-31',...
                                        '02718_2019-05-15_10-44',... % cortical burst of spikewaves at 3000s - mini seizure? % pretty artefacted, but not seen in relevant channels
                                        '02718_2019-05-15_11-52'};
config{2}.directorylist{2}          =  {'02718_2019-05-15_13-52',...
                                        '02718_2019-05-15_15-52',...
                                        '02718_2019-05-15_17-52',... % VERY ARTEFACTED
                                        '02718_2019-05-15_19-52',...
                                        '02718_2019-05-15_21-52',...
                                        '02718_2019-05-15_23-52',...
                                        '02718_2019-05-16_01-52',...
                                        '02718_2019-05-16_03-52',...
                                        '02718_2019-05-16_05-52',...
                                        '02718_2019-05-16_07-52',...
                                        '02718_2019-05-16_09-52',...
                                        '02718_2019-05-16_10-16',...
                                        '02718_2019-05-16_11-15'};   % very artefacted for most of the recording
config{2}.directorylist{3}          =  {'02718_2019-05-16_13-15',...
                                        '02718_2019-05-16_15-15',...
                                        '02718_2019-05-16_17-15',... % very artefacted
                                        '02718_2019-05-16_19-15',... % very artefacted
                                        '02718_2019-05-16_21-15',... % very artefacted for most of the recording
                                        '02718_2019-05-16_23-15',...
                                        '02718_2019-05-17_01-15',...
                                        '02718_2019-05-17_03-15',...
                                        '02718_2019-05-17_05-15',...
                                        '02718_2019-05-17_07-15',...
                                        '02718_2019-05-17_09-15',...
                                        '02718_2019-05-17_10-12'};

config{3}.directorylist             = [];
config{3}.directorylist{1}          = { '02660_2018-11-13_16-08',...
                                        '02660_2018-11-13_18-08',...
                                        '02660_2018-11-13_20-08',...
                                        '02660_2018-11-13_22-08',...
                                        '02660_2018-11-14_00-08',...
                                        '02660_2018-11-14_02-08',...
                                        '02660_2018-11-14_04-08',...
                                        '02660_2018-11-14_06-08',...
                                        '02660_2018-11-14_08-08',...
                                        '02660_2018-11-14_09-29',...
                                        '02660_2018-11-14_10-56',...
                                        '02660_2018-11-14_11-51',...
                                        '02660_2018-11-14_13-51',...
                                        '02660_2018-11-14_15-51'};
config{3}.directorylist{2}          = { '02660_2018-11-14_17-51',...
                                        '02660_2018-11-14_19-51',...
                                        '02660_2018-11-14_21-51',...
                                        '02660_2018-11-14_23-51',...
                                        '02660_2018-11-15_01-51',...
                                        '02660_2018-11-15_03-51',...
                                        '02660_2018-11-15_05-51',...
                                        '02660_2018-11-15_07-51',...
                                        '02660_2018-11-15_09-20',...
                                        '02660_2018-11-15_10-26',...
                                        '02660_2018-11-15_10-52',...
                                        '02660_2018-11-15_11-35',...
                                        '02660_2018-11-15_11-39',...
                                        '02660_2018-11-15_13-41',...
                                        '02660_2018-11-15_15-41'};
config{3}.directorylist{3}          = { '02660_2018-11-15_17-41',... %% LARGELY ARTEFACTED
                                        '02660_2018-11-15_19-41',...
                                        '02660_2018-11-15_21-41',...
                                        '02660_2018-11-15_23-41',... % missing data ShortITS
                                        '02660_2018-11-16_01-41',...
                                        '02660_2018-11-16_03-41',...
                                        '02660_2018-11-16_05-41',...
                                        '02660_2018-11-16_07-41',...
                                        '02660_2018-11-16_09-25',...
                                        '02660_2018-11-16_11-25',...
                                        '02660_2018-11-16_11-56',...
                                        '02660_2018-11-16_13-56',...
                                        '02660_2018-11-16_13-59',...
                                        '02660_2018-11-16_15-59'};

config{4}.directorylist             = [];
config{4}.directorylist{1}          = { '02680_2019-01-15_12-45'...
                                        '02680_2019-01-15_14-45'...
                                        '02680_2019-01-15_15-31'...
                                        '02680_2019-01-15_17-31'...
                                        '02680_2019-01-15_19-31'...
                                        '02680_2019-01-15_21-31'...
                                        '02680_2019-01-15_23-31'...
                                        '02680_2019-01-16_01-31'...
                                        '02680_2019-01-16_03-31'...
                                        '02680_2019-01-16_05-31'...
                                        '02680_2019-01-16_07-31'...
                                        '02680_2019-01-16_09-31'...
                                        '02680_2019-01-16_09-52'...
                                        '02680_2019-01-16_10-58'...
                                        '02680_2019-01-16_11-32'}; % ChWrongSize
config{4}.directorylist{2}          = { '02680_2019-01-16_13-32'...
                                        '02680_2019-01-16_15-32'...
                                        '02680_2019-01-16_17-32'...
                                        '02680_2019-01-16_19-32'... % ends in noise
                                        '02680_2019-01-16_21-32'... % full noise then loose signal (flatline)
                                        '02680_2019-01-16_23-32'... % no signal - flatline - looks hp filtered
                                        '02680_2019-01-17_01-32'... % no signal - flatline - looks hp filtered
                                        '02680_2019-01-17_03-32'... % no signal - flatline - looks hp filtered
                                        '02680_2019-01-17_05-32'... % no signal - flatline - looks hp filtered
                                        '02680_2019-01-17_07-32'... % no signal - flatline - looks hp filtered
                                        '02680_2019-01-17_09-32'... % signal back?
                                        '02680_2019-01-17_10-09'...
                                        '02680_2019-01-17_12-01'...
                                        '02680_2019-01-17_14-01'};
config{4}.directorylist{3}          = { '02680_2019-01-17_16-01'...
                                        '02680_2019-01-17_18-01'...
                                        '02680_2019-01-17_20-01'...
                                        '02680_2019-01-17_22-01'...
                                        '02680_2019-01-18_00-01'...
                                        '02680_2019-01-18_02-01'...
                                        '02680_2019-01-18_04-01'...
                                        '02680_2019-01-18_06-01'...
                                        '02680_2019-01-18_08-01'...
                                        '02680_2019-01-18_09-41'...
                                        '02680_2019-01-18_11-41'...
                                        '02680_2019-01-18_13-41'...
                                        '02680_2019-01-18_14-08'};

config{5}.directorylist             = [];
config{5}.directorylist{1}          = { '02689_2019-02-12_13-23'...
                                        '02689_2019-02-12_15-23'...
                                        '02689_2019-02-12_17-23'...
                                        '02689_2019-02-12_19-23'...
                                        '02689_2019-02-12_21-23'...
                                        '02689_2019-02-12_23-23'...
                                        '02689_2019-02-13_01-23'...
                                        '02689_2019-02-13_03-23'...
                                        '02689_2019-02-13_05-23'...
                                        '02689_2019-02-13_07-23'...
                                        '02689_2019-02-13_09-23'...
                                        '02689_2019-02-13_09-25'...
                                        '02689_2019-02-13_10-37'...
                                        '02689_2019-02-13_12-37'};
                                       % '02689_2019-02-13_10-32'... short
                                       % and full of artefacts
config{5}.directorylist{2}          = { '02689_2019-02-13_14-37'...
                                        '02689_2019-02-13_16-37'...
                                        '02689_2019-02-13_18-37'...
                                        '02689_2019-02-13_20-37'...
                                        '02689_2019-02-13_22-37'...
                                        '02689_2019-02-14_00-37'...
                                        '02689_2019-02-14_02-37'...
                                        '02689_2019-02-14_04-37'...
                                        '02689_2019-02-14_06-37'...
                                        '02689_2019-02-14_08-37'...
                                        '02689_2019-02-14_10-25'...
                                        '02689_2019-02-14_12-25'};
config{5}.directorylist{3}          = { '02689_2019-02-14_14-25'...
                                        '02689_2019-02-14_16-25'...
                                        '02689_2019-02-14_18-25'...
                                        '02689_2019-02-14_20-25'...
                                        '02689_2019-02-14_22-25'...
                                        '02689_2019-02-15_00-25'...
                                        '02689_2019-02-15_02-25'...
                                        '02689_2019-02-15_04-25'...
                                        '02689_2019-02-15_06-25'...
                                        '02689_2019-02-15_08-25'...
                                        '02689_2019-02-15_09-54'...
                                        '02689_2019-02-15_11-54'...
                                        '02689_2019-02-15_13-54'};

config{6}.directorylist             = [];
config{6}.directorylist{1}          = { '02651_2018-10-16_15-31'...
                                        '02651_2018-10-16_17-31'... % fully artefacted
                                        '02651_2018-10-16_19-31'... % fully artefacted
                                        '02651_2018-10-16_21-31'...
                                        '02651_2018-10-16_23-31'...
                                        '02651_2018-10-17_01-31'...
                                        '02651_2018-10-17_03-31'...
                                        '02651_2018-10-17_05-31'... % mostly artefacted
                                        '02651_2018-10-17_07-31'... % regularly loose wire
                                        '02651_2018-10-17_09-31'...
                                        '02651_2018-10-17_11-22'... % short and fully artefacted
                                        '02651_2018-10-17_11-46'...
                                        '02651_2018-10-17_13-46'};  % very regular clipping in mHaT2
config{6}.directorylist{2}          = { '02651_2018-10-17_15-46'...
                                        '02651_2018-10-17_17-46'... % fully artefacted
                                        '02651_2018-10-17_19-46'...
                                        '02651_2018-10-17_21-46'...
                                        '02651_2018-10-17_23-46'...
                                        '02651_2018-10-18_01-46'...
                                        '02651_2018-10-18_03-46'...
                                        '02651_2018-10-18_05-46'...
                                        '02651_2018-10-18_07-46'...
                                        '02651_2018-10-18_09-46'...
                                        '02651_2018-10-18_11-06'...
                                        '02651_2018-10-18_12-01'...
                                        '02651_2018-10-18_14-01'};
config{6}.directorylist{3}          = { '02651_2018-10-18_16-01'...
                                        '02651_2018-10-18_18-01'...
                                        '02651_2018-10-18_20-01'...
                                        '02651_2018-10-18_22-01'...
                                        '02651_2018-10-19_00-01'...
                                        '02651_2018-10-19_02-01'...
                                        '02651_2018-10-19_04-01'...
                                        '02651_2018-10-19_06-01'...
                                        '02651_2018-10-19_08-01'...
                                        '02651_2018-10-19_10-01'...
                                        '02651_2018-10-19_10-45'...
                                        '02651_2018-10-19_12-45'...
                                        '02651_2018-10-19_14-45'};

config{7}.directorylist             = [];
config{7}.directorylist{1}          = { '02619_2018-07-03_16-40'...
                                        '02619_2018-07-03_18-40'...
                                        '02619_2018-07-03_20-40'... % very artefacted - macro needs rereferencing (e.g. bipolar), micro constant movement artefacts
                                        '02619_2018-07-03_22-40'... % amygdala spike-waves starting
                                        '02619_2018-07-04_00-40'...
                                        '02619_2018-07-04_02-40'...
                                        '02619_2018-07-04_04-40'... % bursty neuron on mAmT2_6
                                        '02619_2018-07-04_06-40'...
                                        '02619_2018-07-04_08-40'...
                                        '02619_2018-07-04_12-45'... % pretty artefacted
                                        '02619_2018-07-04_14-45'};  % pretty artefacted
config{7}.directorylist{2}          = { '02619_2018-07-04_16-45'...
                                        '02619_2018-07-04_18-45'...
                                        '02619_2018-07-04_20-45'...
                                        '02619_2018-07-04_22-45'...
                                        '02619_2018-07-05_00-45'...
                                        '02619_2018-07-05_02-45'...
                                        '02619_2018-07-05_04-45'...
                                        '02619_2018-07-05_06-45'...
                                        '02619_2018-07-05_08-45'...
                                        '02619_2018-07-05_10-43'...
                                        '02619_2018-07-05_12-13'...
                                        '02619_2018-07-05_14-13'...
                                        '02619_2018-07-05_16-13'};
config{7}.directorylist{3}          = { '02619_2018-07-05_18-13'...
                                        '02619_2018-07-05_20-13'...
                                        '02619_2018-07-05_22-13'... % DiffSTartMacSync & short ITS
                                        '02619_2018-07-06_00-13'...
                                        '02619_2018-07-06_02-13'...
                                        '02619_2018-07-06_04-13'...
                                        '02619_2018-07-06_06-13'... % SEIZURES
                                        '02619_2018-07-06_08-13'...
                                        '02619_2018-07-06_10-13'...
                                        '02619_2018-07-06_10-47'...
                                        '02619_2018-07-06_11-11'...
                                        '02619_2018-07-06_11-57'...
                                        '02619_2018-07-06_12-30'...
                                        '02619_2018-07-06_14-30'...
                                        '02619_2018-07-06_16-30'};

config{8}.directorylist{1}          = { '02578_2018-02-13_12-11'...
                                        '02578_2018-02-13_14-11'...
                                        '02578_2018-02-13_16-11'...
                                        '02578_2018-02-13_18-11'...
                                        '02578_2018-02-13_20-11'...
                                        '02578_2018-02-13_22-11'...
                                        '02578_2018-02-14_00-11'...
                                        '02578_2018-02-14_02-11'...
                                        '02578_2018-02-14_04-11'...
                                        '02578_2018-02-14_06-11'...
                                        '02578_2018-02-14_08-11'...
                                        '02578_2018-02-14_09-24'...
                                        '02578_2018-02-14_09-41'...
                                        '02578_2018-02-14_09-56'...
                                        '02578_2018-02-14_11-56'};
                                        
config{8}.directorylist{2}          = { '02578_2018-02-14_13-56'...
                                        '02578_2018-02-14_15-56'...
                                        '02578_2018-02-14_17-56'...
                                        '02578_2018-02-14_19-56'...
                                        '02578_2018-02-14_21-56'...
                                        '02578_2018-02-14_23-56'...
                                        '02578_2018-02-15_01-56'...
                                        '02578_2018-02-15_03-56'...
                                        '02578_2018-02-15_05-56'...
                                        '02578_2018-02-15_07-56'...
                                        '02578_2018-02-15_09-56'...
                                        '02578_2018-02-15_10-29'...
                                        '02578_2018-02-15_10-52'...
                                        '02578_2018-02-15_11-51'};
                                        
config{8}.directorylist{3}          = { '02578_2018-02-15_13-51'...
                                        '02578_2018-02-15_15-51'...
                                        '02578_2018-02-15_17-51'...
                                        '02578_2018-02-15_19-51'...
                                        '02578_2018-02-15_21-51'...
                                        '02578_2018-02-15_23-51'...
                                        '02578_2018-02-16_01-51'...
                                        '02578_2018-02-16_03-51'...
                                        '02578_2018-02-16_05-51'...
                                        '02578_2018-02-16_07-51'...
                                        '02578_2018-02-16_09-51'...
                                        '02578_2018-02-16_10-45'...
                                        '02578_2018-02-16_10-58'...
                                        '02578_2018-02-16_12-58'};
                                    
                                    
% config{}.directorylist             = [];
% config{}.directorylist{1}          = { '02599_2018-04-24_12-36'...
%                                         '02599_2018-04-24_14-36'...
%                                         '02599_2018-04-24_16-36'...
%                                         '02599_2018-04-24_18-36'...
%                                         '02599_2018-04-24_20-36'...
%                                         '02599_2018-04-24_22-36'...
%                                         '02599_2018-04-25_00-36'...
%                                         '02599_2018-04-25_02-36'...
%                                         '02599_2018-04-25_04-36'...
%                                         '02599_2018-04-25_06-36'...
%                                         '02599_2018-04-25_08-36'...
%                                         '02599_2018-04-25_09-24'...
%                                         '02599_2018-04-25_11-24'...
%                                         '02599_2018-04-25_13-24'};
% config{}.directorylist{2}          = { '02599_2018-04-25_15-24'...
%                                         '02599_2018-04-25_17-24'...
%                                         '02599_2018-04-25_19-24'...
%                                         '02599_2018-04-25_21-24'...
%                                         '02599_2018-04-25_23-24'...
%                                         '02599_2018-04-26_01-24'...
%                                         '02599_2018-04-26_03-24'...
%                                         '02599_2018-04-26_05-24'...
%                                         '02599_2018-04-26_07-24'...
%                                         '02599_2018-04-26_09-24'...
%                                         '02599_2018-04-26_09-32'...
%                                         '02599_2018-04-26_11-32'...
%                                         '02599_2018-04-26_13-32'};
% config{}.directorylist{3}          = { '02599_2018-04-26_14-13'...
%                                         '02599_2018-04-26_14-53'...
%                                         '02599_2018-04-26_16-53'...
%                                         '02599_2018-04-26_18-53'...
%                                         '02599_2018-04-26_20-53'...
%                                         '02599_2018-04-26_22-53'...
%                                         '02599_2018-04-27_00-53'...
%                                         '02599_2018-04-27_02-53'...
%                                         '02599_2018-04-27_04-53'...
%                                         '02599_2018-04-27_06-53'...
%                                         '02599_2018-04-27_08-53'...
%                                         '02599_2018-04-27_09-23'...
%                                         '02599_2018-04-27_11-23'...
%                                         '02599_2018-04-27_13-23'};



%%
















% config{1}.hyp.notcontains         = {"ADStartLoss","ADEndLoss","TTL","StartRecord","StopRecord","NLXEvent","BAD"};

% config{1}.align.name                = {'Hspike','SpikeDetect'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.channel             = {'_HaT2_1','_HaT2_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.abs                 = {'no','no'};
% config{1}.align.method              = {'crawlback','min'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{1}.align.filter              = {'bp','bp'};
% config{1}.align.freq                = {[1, 10],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{1}.align.hilbert             = {'no','no'};
% config{1}.align.thresh              = [1, 0];
% config{1}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{1}        = [-0.05,  0.150];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{2}        = [-0.05,  0.05];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
%
% config{1}.circus.channel            = {'mHaT2_1','mHaT2_3','mHaT2_4','mHaT2_6','mHaT2_8'};
% config{1}.circus.reref              = 'yes';
% config{1}.circus.refchan            = 'mHaT2_2';
% config{1}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'hspike', 'SpykingCircus');
% config{1}.circus.hpfilter           = 'no';
% config{1}.circus.hpfreq             = 0;
% config{1}.circus.postfix            = '-1';
%
% config{1}.TFR.channel               = {'_HaT2_1','_HaT2_6'};
%
% config{1}.spike.slidestep           = 0.001;
% config{1}.spike.toispikerate{1}     = [-2, 1];
% config{1}.spike.resamplefs          = 1000;
% config{1}.spike.width               = 15;
% config{1}.spike.ISIbins             = 0 : 0.05 : 0.150;
%
% config{1}.stats.bltoi{1}            = [-0.5, -0.1];
% config{1}.stats.actoi{1}            = [0,    1];
% config{1}.stats.bltoi{2}            = [-0.5, -0.1];
% config{1}.stats.actoi{2}            = [0,    1];
% config{1}.stats.alpha               = 0.025;

% config{1}.spikedetect.LS            = 5;    % Left half-wave slope; default: 7
% config{1}.spikedetect.RS            = 5;    % Right half-wave slope; default: 7
% config{1}.spikedetect.TAMP          = 400;  % Total amplitude; default: 600
% config{1}.spikedetect.LD            = 1;    % Left half-wave duration; default = 10
% config{1}.spikedetect.RD            = 1;    % Right half-wave duration; default = 10
% config{1}.spikedetect.STDCoeff      = 4;    % Chebyshev inequality coefficient (distance from centre point or mean); default 4
% config{1}.spikedetect.SCALE         = 70;   % Scaling parameter
% config{1}.spikedetect.BlockSize     = 1;    % Data processing block size in minutes
% config{1}.spikedetect.TroughSearch  = 40;   % distance in ms to search for a trough on each side of a detected peak
% config{1}.spikedetect.FilterSpec    = [20; 50; 1; 35;];
%
% config{2}.align.name                = {'SpikeHaT1_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.align.channel             = {'_HaT1_3','_HaT1_3'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.align.abs                 = {'no','no'};
% config{2}.align.method              = {'nearest','max'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{2}.align.filter              = {'bp','bp'};
% config{2}.align.freq                = {[1, 40],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{2}.align.hilbert             = {'no','no'};
% config{2}.align.thresh              = [0, 0];
% config{2}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiactive{1}        = [-0.05,  0.150];                                         % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{2}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiactive{2}        = [-0.05,  0.05];                                          % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{2}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% %
% config{2}.circus.channel            = {'mHaT1_7'};
% config{2}.circus.reref              = 'no';
% config{2}.circus.refchan            = 'mHaT1_1';
% config{2}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'hspike', 'SpykingCircus');
% config{2}.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
% config{2}.circus.hpfreq             = 0; % even when not using
% config{2}.circus.postfix            = '-1'; % after using circus-gui-matlab's SAVE number
