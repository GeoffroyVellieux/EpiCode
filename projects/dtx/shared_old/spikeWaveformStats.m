function 	stats = spikeWaveformStats(cfg, SpikeWaveforms, force)

% SPIKEWAVEFORMSTATS computes average and variance of the spike waveform.
% Then, the average is interpolated and the following parameters are
% computed : amplitude, halfwidth, peaktrough and trouhpeak.
% Input spike waveform data is get from readSpikeWaveforms.m.
%
% Use as :
%   stats = spikeWaveformStats(cfg, SpikeWaveforms, force)
%
% Note :
% - spikes must be peak-aligned to zero (automatically done by SpyKING-Circus)
% - spikes can be positive or negative
%
% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.


fname_out = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'spikewaveformstats.mat'));

if exist(fname_out, 'file') && force == false
    fprintf('Loading precomputed spike waveform stats\n');
    load(fname_out, 'stats');
    return
else
    fprintf('(re-) computing spike waveform stats\n');
end

% add external/intersections path if not already
mfile_name = mfilename('fullpath');
pathstr    = fileparts(fileparts(mfile_name));
pathCell   = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    onPath = any(strcmpi(fullfile(pathstr, ['external', filesep, 'intersections']), pathCell));
else
    onPath = any(strcmp(fullfile(pathstr, ['external', filesep, 'intersections']), pathCell));
end
if ~onPath
    addpath(fullfile(pathstr, ['external', filesep, 'intersections']));
end

for ipart = 1:size(SpikeWaveforms)
    for markername = string(fieldnames(SpikeWaveforms{ipart}))'
        ft_progress('init', 'text');
        for icluster = 1:size(SpikeWaveforms{ipart}.(markername), 2)

            ft_progress(0, 'Compute waveform for part %d, %s, cluster %d/%d', ipart, char(markername), icluster, size(SpikeWaveforms{ipart}.(markername), 2));
            ok = true;

            if ~isempty(SpikeWaveforms{ipart}.(markername){icluster})

                %average waveform
                evalc('waveformavg = ft_timelockanalysis([], SpikeWaveforms{ipart}.(markername){icluster});');%remove text output

                %interpolate the averaged waveform to have more precise values
                time_temp        = linspace(waveformavg.time(1),waveformavg.time(end),1000);
                avg_temp         = pchip(waveformavg.time,waveformavg.avg,time_temp);
                var_temp         = pchip(waveformavg.time,waveformavg.var,time_temp);
                waveformavg.time = time_temp;
                waveformavg.avg  = avg_temp;
                waveformavg.var  = var_temp;
                %scatter(waveformavg.time, waveformavg.avg, '.');

                %search if AP is positive or negative
                flip = sign(waveformavg.avg(nearest(waveformavg.time, 0)));

                %amplitude
                bl               = min(waveformavg.avg.*flip);
                [peak, peak_idx] = max(waveformavg.avg.*flip);
                amplitude.val    = peak-bl;
                amplitude.x      = waveformavg.time(peak_idx);
                amplitude.y      = waveformavg.avg(peak_idx);
                %scatter(amplitude.x, amplitude.y, 'xk');

                %halfwidth
                halfamp     = bl+(peak-bl)/2;
                x  = waveformavg.time;
                y1 = waveformavg.avg;
                y2 = ones(size(x)) .* halfamp .* flip;
                [x_intersect, y_intersect] = intersections(x,y1,x,y2,true);
                
                if length(x_intersect) >= 2 && any(x_intersect <0) && any(x_intersect >0)
                  
                    idx = find(x_intersect <0, 1, 'last');
                    halfwidth.val = diff(x_intersect([idx, idx+1]));
                    halfwidth.x   = x_intersect([idx, idx+1]);
                    halfwidth.y   = y_intersect([idx, idx+1]);
                    %scatter(halfwidth.x, halfwidth.y, 'xk');

                    % Find throughs
                    [Yneg,Xneg_temp] = findpeaks(-waveformavg.avg.*flip,waveformavg.time);
                    
                    if length(Xneg_temp) >= 2 && any(Xneg_temp-amplitude.x < 0) && any(Xneg_temp-amplitude.x > 0)
                        % Search first through before and after peak
                        [Xneg(1),x_idx] = max(Xneg_temp(Xneg_temp-amplitude.x < 0));
                        Xneg(2)         = min(Xneg_temp(Xneg_temp-amplitude.x > 0));
                        Yneg            = Yneg([x_idx, x_idx+1]);

                        peaktrough.val  = abs(amplitude.x-Xneg(1));
                        peaktrough.x    = [Xneg(1)  amplitude.x];
                        peaktrough.y    = [-Yneg(1)*flip amplitude.y];
                        troughpeak.val  = abs(amplitude.x-Xneg(2));
                        troughpeak.x    = [amplitude.x Xneg(2)];
                        troughpeak.y    = [amplitude.y -Yneg(2)*flip];
                        %scatter(peaktrough.x, peaktrough.y, 'xk');
                        %scatter(troughpeak.x, troughpeak.y, 'xk');
                        
                    else
                        ok = false;
                    end
                else
                    ok = false;
                end
            else
                ok = false;
            end %isempty

            %store for output
            
            if isempty(SpikeWaveforms{ipart}.(markername){icluster})
                stats{ipart}.(markername).label{icluster}          = [];
                stats{ipart}.(markername).waveformavg{icluster}    = [];
            else
                stats{ipart}.(markername).label{icluster}          = SpikeWaveforms{ipart}.(markername){icluster}.label{1};
                stats{ipart}.(markername).waveformavg{icluster}    = waveformavg;
            end
            if ok
                stats{ipart}.(markername).amplitude.val(icluster)  = amplitude.val;
                stats{ipart}.(markername).amplitude.x(icluster)    = amplitude.x;
                stats{ipart}.(markername).amplitude.y(icluster)    = amplitude.y;
                stats{ipart}.(markername).halfwidth.val(icluster)  = halfwidth.val;
                stats{ipart}.(markername).halfwidth.x(icluster,:)  = halfwidth.x;
                stats{ipart}.(markername).halfwidth.y(icluster,:)  = halfwidth.y;
                stats{ipart}.(markername).peaktrough.val(icluster) = peaktrough.val;
                stats{ipart}.(markername).peaktrough.x(icluster,:) = peaktrough.x;
                stats{ipart}.(markername).peaktrough.y(icluster,:) = peaktrough.y;
                stats{ipart}.(markername).troughpeak.val(icluster) = troughpeak.val;
                stats{ipart}.(markername).troughpeak.x(icluster,:) = troughpeak.x;
                stats{ipart}.(markername).troughpeak.y(icluster,:) = troughpeak.y;
            else
                stats{ipart}.(markername).amplitude.val(icluster)  = nan;
                stats{ipart}.(markername).amplitude.x(icluster)    = nan;
                stats{ipart}.(markername).amplitude.y(icluster)    = nan;
                stats{ipart}.(markername).halfwidth.val(icluster)  = nan;
                stats{ipart}.(markername).halfwidth.x(icluster,:)  = nan;
                stats{ipart}.(markername).halfwidth.y(icluster,:)  = nan;
                stats{ipart}.(markername).peaktrough.val(icluster) = nan;
                stats{ipart}.(markername).peaktrough.x(icluster,:) = nan;
                stats{ipart}.(markername).peaktrough.y(icluster,:) = nan;
                stats{ipart}.(markername).troughpeak.val(icluster) = nan;
                stats{ipart}.(markername).troughpeak.x(icluster,:) = nan;
                stats{ipart}.(markername).troughpeak.y(icluster,:) = nan;
            end
        end
        ft_progress('close');
    end
end

save(fname_out, 'stats', '-v7.3');
