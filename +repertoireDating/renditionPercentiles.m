function RP = renditionPercentiles(NNids, epochs, varargin)
    % COMPUTEPERCENTILES This function computes repertoire dating percentiles for 
    % each rendition based on a nearest neighbor graph defined for a number of 
    % measurements/renditions (size(NNids, 1)) and measurement times (epochs). 
    %
    % If doPlot is true, rendition percentiles are plotted (sorted by the 
    % middle percentile; e.g. 3 for 5 given percentiles, 2 for 4 percentiles, etc.).
    % 
    % Inputs
    % ------
    % NNids: N x #Neighbours (integer)
    % - defines the k-NN graph
    % - entires of this matrix of indices must be in 1..N (N is the total number of datapoints)
    % - row j in this matrix contains the indices of the k nearest neighbours of datapoint j
    % - points should not be their own neighbors
    % - this input can have different dimension if an optional input 'valid' is specified
    % - (see below)
    % 
    % epochs:  N x 1 (double, e.g. day, ideally without gaps)
    % 
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    %
    % percentiles: 1 x P vector, default: [5, 25, 50, 75, 95].
    %
    % valid: N x 1 (logical)
    % - marks all datapoints used as querry points for the k-NN based analysis
    % - this does not restrict the neighbours of querry points but only the set of
    %   querry points itself
    % - can be empty (which is interpreted as true(N, 1))
    % - if sum(valid) < N and size(NNids, 1) == sum(valid), then the jth row of NNids
    %   contains the neighbours of datapoint valid_ids(j), where valid_ids = find(valid)
    %   Note: in that case, the entries of NNids are still interpreted as referring to 1..N
    %
    % doPlot: logical
    % - if true, rendition percentiles are plotted.
    %
    % Outputs
    % -------
    % RPD: numel(prct) x sum(valid)
    % - the repertoire dating percentiles for each rendition
    % - if no 'valid' parameter was provided: sum(valid) = N
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
    % Kollmorgen, S., Hahnloser, R.H.R. & Mante, V. Nearest neighbours reveal
    % fast and slow components of motor learning. Nature 577, 526-530 (2020).
    % https://doi.org/10.1038/s41586-019-1892-x
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
        {'percentiles', 'valid', 'doPlot', 'noiseLevel'},...
        {[5, 25, 50, 75, 95], [], false, 0.1},...
        'repertoireDating.renditionPercentiles');
    
    if isempty(S.valid)
        S.valid = true(size(NNids, 1), 1);
    end
       
    %% Compute percentiles
    % Obtain neighborhood labels for each valid rendition
    if sum(S.valid) == size(NNids, 1)
        nhLabels = epochs(NNids);
    else
        assert(size(NNids, 1) == numel(epochs));
        nhLabels = epochs(NNids(S.valid, :));
    end
    
    RP = prctile(nhLabels, S.percentiles, 2);
       
    %% Plot percentiles
    if S.doPlot
        pfcn = repertoireDating.internal.bda_figure('Rendition neighborhood label percentiles', [1, 1], 1);
        pfcn(1, 1, ''); % Create a figure panel
        [~, idx] = sort(RP(:, ceil(numel(S.percentiles)/2)));
        
        for kk = 1:numel(S.percentiles)
            if kk == ceil(numel(S.percentiles)/2)
                plot(RP(idx, kk)+randn(numel(idx), 1)*S.noiseLevel, 1:size(RP, 1), '.r', 'MarkerSize', 0.5);
            else
                plot(RP(idx, kk)+randn(numel(idx), 1)*S.noiseLevel, 1:size(RP, 1), '.k', 'MarkerSize', 0.5);
            end
        end
        xlabel('Neighborhood Label');
        ylabel(sprintf('Rendition Index (sorted by %ith percentile)', S.percentiles(ceil(numel(S.percentiles)/2))));
        title(['percentiles: ' num2str(S.percentiles)]);
        axis tight;
        repertoireDating.internal.bda_formatFigure(gcf, 1, 10, 'off');
    end
    
end

