function wrong_markers = verifymarkers(cfg, MuseStruct, ipart)

% cfg.verifymarkers.one_off
% cfg.verifymarkers.start
% cfg.verifymarkers.end

for idir = 1:length(MuseStruct{ipart})
    
    wrong_markers{idir}.prefix = cfg.prefix(1:end-1);
    wrong_markers{idir}.directory = cfg.directorylist{ipart}{idir};
 
    %% one off markers
    for markername = string(cfg.verifymarkers.one_off)
        markername_wrong = sprintf('%s__START__', markername);
        if isfield(MuseStruct{ipart}{idir}.markers,char(markername_wrong))
            wrong_markers{idir}.(markername_wrong) = MuseStruct{ipart}{idir}.markers.(markername_wrong).synctime';
        end
    end
        
    
    %% start/end markers
    for imarker = 1:size(cfg.verifymarkers.start, 2)

        MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker})          = ft_getopt(MuseStruct{ipart}{idir}.markers, cfg.verifymarkers.start{imarker}, []);
        MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker})            = ft_getopt(MuseStruct{ipart}{idir}.markers, cfg.verifymarkers.end{imarker}, []);
        
        MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime = ft_getopt(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}), 'synctime', []);
        MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime   = ft_getopt(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}), 'synctime', []);
        
        if length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime) ~= length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime)
            wrong_markers{idir}.(cfg.verifymarkers.start{imarker}) = sprintf('%d start markers', length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime));
            wrong_markers{idir}.(cfg.verifymarkers.end{imarker})   = sprintf('%d end markers', length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime));
        end
    end
    
    if numel(fieldnames(wrong_markers{idir})) > 2
        ft_warning('Wrong marker(s) found in %s, %s.\n', cfg.prefix(1:end-1), cfg.directorylist{ipart}{idir});
    else
        wrong_markers{idir}.All_markers_are_good = 'OK';
    end

end


