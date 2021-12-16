function MuseStruct = merge_MuseMarkers(cfg, MuseStruct)

% cfg.newmarker = 'new_name';
% cfg.markers_to_merge = {'marker1', 'marker2', 'marker3'};

for ipart = 1:size(MuseStruct, 2)
    for idir = 1:size(MuseStruct{ipart}, 2)
        synctime = [];
        clock = [];
        for markername = string(cfg.markers_to_merge)
            if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'synctime')
                continue
            end
            synctime = [synctime, MuseStruct{ipart}{idir}.markers.(markername).synctime];
            clock    = [clock, MuseStruct{ipart}{idir}.markers.(markername).clock];
        end
        synctime = sort(synctime);
        clock    = sort(clock);
        MuseStruct{ipart}{idir}.markers.(cfg.newmarker).synctime = synctime;
        MuseStruct{ipart}{idir}.markers.(cfg.newmarker).clock = clock;
    end
end