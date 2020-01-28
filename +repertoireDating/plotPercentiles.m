function [avg_nDays, shift, span] = plotPercentiles(prctPerDay_avg, epochs, subEpochs,...
        epochs_forAvg, varargin)
    % PLOTPERCENTILES Plots repertoire dating percentiles.
    %
    % Inputs
    % ------
    % prctPerDay_avg: numel(prct) X numel(epochs) X numel(subEpochs)
    %
    % subEpochs: 1..nDivisions/epoch
    %
    % epochs: [e_0, e_0+1, ..., e_max]
    %
    % epochs_forAvg: subset of epochs
    %
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    %
    % nEpochs_in_avg: how many consecutive epochs to average (default: 2)
    %
    % percentiles: e.g. [5, 25, 50, 75, 95]
    %
    % Example
    % -------
    %
    % % Create random data with linear drift:
    % epochs = 40; 
    % subEpochs = 10; 
    % samples = 100; 
    % dim = 25;  
    % N = epochs * subEpochs * samples;
    % productionTime = 1:N;
    % epoch = ceil(productionTime/(subEpochs*samples));
    % subEpoch = floor(mod(productionTime-1, (subEpochs*samples))/samples) + 1;
    % data = randn(N, dim) + (1:N)'/N;
    % 
    % % Compute k-NN graph
    % NNids = knnsearch(data, data, 'K', 50);     % or use your preferred nearest neighbor searcher 
    % NNids = NNids(:, 2:end);                      % points should not be their own neighbors
    % 
    % % Compute and plot repertoire dating pecentiles 
    % [RPD, RPD_epoch, RPD_subEpoch] = repertoireDating.percentiles(NNids, epoch, subEpoch);
    % repertoireDating.plotPercentiles(RPD, RPD_epoch, RPD_subEpoch, 16:25);
    % 
    % % Compute and plot rendition percentiles for all datapoints from epoch 20
    % repertoireDating.renditionPercentiles(NNids, epoch, 'valid', epoch == 20, 'doPlot', true);
    % 
    % % Compute and plot rendition percentiles for all datapoints
    % repertoireDating.renditionPercentiles(NNids, epoch, 'doPlot', true);
    %
    % 
    % ---
    % Copyright (C) 2020 University Zurich, Sepp Kollmorgen
    % 
    % Reference (please cite):
    % Nearest neighbours reveal fast and slow components of motor learning.
    % Kollmorgen, S., Hahnloser, R.H.R.; Mante, V.
    % Nature (2020) doi:10.1038/s41586-019-1892-x
    % 
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU Affero General Public License for more details.
    %
    % You should have received a copy of the GNU Affero General Public License
    % along with this program (see LICENSE file).  If not, see <https://www.gnu.org/licenses/>.
    %
    % Repertoire-Dating on Github: <a href="matlab:web('https://github.com/skollmor/repertoireDating', '-browser')">https://github.com/skollmor/repertoireDating</a>
    % Dataspace on Github: <a href="matlab:web('https://github.com/skollmor/dspace', '-browser')">https://github.com/skollmor/dspace</a>

    
    S = repertoireDating.internal.parse_optional_arguments(varargin,...
        {'nEpochs_in_avg' 'percentiles', 'doPlot'},...
        {2, [5, 25, 50, 75, 95], true},...
        'repertoireDating.plotPercentiles');
    
    % Make sure subEpochs and epochs are row vectors
    if numel(subEpochs) ~= size(subEpochs, 2)
        subEpochs = reshape(subEpochs, 1, []);
    end
    if numel(epochs) ~= size(epochs, 2)
        epochs = reshape(epochs, 1, []);
    end
    
    if S.doPlot
        % Plot repertoire dating percentiles
        pfcn = repertoireDating.internal.bda_figure('', [5, 5], 1);
        pfcn(1, 1, '', 4, 5); hold all;
        for u = epochs
            for s = subEpochs
                for kk = 1:numel(S.percentiles)
                    if S.percentiles(kk) == 50
                        color = [1,0,0];
                    else
                        color = [0,0,0];
                    end
                    line(0.1 + 0.8*([s, s+1]-1)/max(subEpochs)+u,...
                        [1, 1] * (prctPerDay_avg(kk, epochs == u, subEpochs == s)), 'color', color);
                    if find(s == subEpochs) < numel(subEpochs)
                        line(0.1 + [0.8, 0.8] * (s)/max(subEpochs)+u,...
                            [(prctPerDay_avg(kk, epochs == u, find(subEpochs == s))),...
                            (prctPerDay_avg(kk, epochs == u, find(subEpochs == s)+1))], 'color', color); %#ok<FNDSB>
                    end
                end
            end
        end
        xlim([min(epochs), max(epochs)+1]);
        ylim([min(epochs), max(epochs)+1]);
        xlabel('production time (epochs)');
        ylabel('repertoire time (epochs)');
        set(gca, 'DataAspectRatio', [1 1 1]);
        
        %% Compute Averages
        pfcn(5, 1, '', 1, 4);
        
    end
    
    % cut out desired stage levels
    prctPerDay_avg_cutout = prctPerDay_avg(:, ismember(epochs, epochs_forAvg), :);
    avg_nDays = cell(1, S.nEpochs_in_avg);
    
    for jj = 1:S.nEpochs_in_avg
        avg_nDays{jj} = zeros(numel(S.percentiles), numel(subEpochs));
    end
    for kk = 1:(size(prctPerDay_avg_cutout, 2) - S.nEpochs_in_avg + 1)
        
        % Sum up repertoire dating percentiles
        for jj = 1:S.nEpochs_in_avg
            avg_nDays{jj} = avg_nDays{jj} + permute(prctPerDay_avg_cutout(:, kk+jj-1, :), [1, 3, 2]);
        end
        
    end
    med = NaN(1, S.nEpochs_in_avg);
    for jj = 1:S.nEpochs_in_avg
        avg_nDays{jj} = avg_nDays{jj}/(size(prctPerDay_avg_cutout, 2) - S.nEpochs_in_avg + 1);
        
        % compute medians for each averaged stage
        if ismember(50, S.percentiles)
            % take median over 50th percentile
            med(jj) = median(avg_nDays{jj}(S.percentiles == 50, :));
        else
            % take median over all percentiles
            med(jj) = median(reshape(avg_nDays{jj}(:, :), 1, []));
        end
    end
    
    % subtract overall mean of medians of 50th percentile from all averaged stages
    for jj = 1:S.nEpochs_in_avg
        avg_nDays{jj} = avg_nDays{jj} - mean(med);
    end
    
    if S.doPlot
        
        fprintf('Epochs in average: %i\n', size(prctPerDay_avg_cutout, 2));
        for u = 0:(S.nEpochs_in_avg-1)
            for s = subEpochs
                for kk = 1:numel(S.percentiles)
                    if S.percentiles(kk) == 50
                        color = [1,0,0];
                    else
                        color = [0,0,0];
                    end
                    
                    line(0.1 + 0.8*([s, s+1]-1)/max(subEpochs) + u,...
                        [1, 1] * avg_nDays{u+1}(kk, subEpochs == s), 'color', color);
                    
                    if find(s == subEpochs) < numel(subEpochs)
                        line(0.1 + [0.8, 0.8] * (s)/max(subEpochs) + u,...
                            [avg_nDays{u+1}(kk, subEpochs == s),
                            avg_nDays{u+1}(kk, find(subEpochs == s)+1)], 'color', color);
                    end
                end
            end
        end
        
        xlabel('production time (k, k+1)');
        ylabel('avg. repertoire time (mean subtracted)');
        %set(gcf, 'position', [14          69        1329        1035]);
        
        
        pfcn(5, 5, 'c');
        
    end
    
    if ismember(50, S.percentiles)
        ppi = find(S.percentiles == 50);
    else
        ppi = S.percentiles(floor(end/2));
    end
    
    
    
    % span is defined as (day1, late) - (day1, early)
    span = avg_nDays{1}(ppi, end) - avg_nDays{1}(ppi, 1);
    % shift is defined as (day2, early) - (day1, late)
    shift = avg_nDays{2}(ppi, 1) - avg_nDays{1}(ppi, end);
    
    if S.doPlot
        
        plot(shift, span, 'xk');
        xlabel(sprintf('shift (Pct: %.0f)', S.percentiles(ppi)));
        ylabel('span');
        mx = max(abs(shift), 2);
        my = max(abs(span), 2);
        xlim([-mx, mx]);
        ylim([-my, my]);
        repertoireDating.internal.bda_formatFigure(gcf, 1.5);
        
    end
end

