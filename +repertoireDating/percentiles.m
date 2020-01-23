function [RPD, uEpochs, uSubEpochs] = percentiles(NNids, epochs, subEpochs, varargin)
    % PERCENTILES This function computes repertoire dating percentiles.
    % 
    % Inputs
    % ------
    % NNids: N x #Neighbours (integer)
    % - defines the k-NN graph
    % - entires of this matrix of indices must be in 1..N (N is the total number of datapoints)
    % - row j in this matrix contains the indices of the k nearest neighbours of datapoint j.
    % - points should not be their own neighbors
    %
    % epochs: N x 1 (double, e.g. day)
    %
    % subEpochs: N x 1 (double, e.g. hour within day, same binning in every epoch)
    % - if epochs are not further subdivided, subEpochs can be omitted (set to []) 
    %   in that case the output will have dimensions numel(prct) x numel(epochLevels) x 1
    %
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    % 
    % style: string 
    % - 'mode a': procedure used in the paper
    % - 'mean based': Uses the mean neighborhood production time and pools 
    %   over those means
    %
    % uEpochs: K1 x 1 (double)
    % - provide if a fixed set of epochs is to be used (even if
    %   some do not occur in the dataset)
    % 
    % uSubEpochs: K2 x 1 (double)
    % - provide if a fixed set of epochs is to be used (even if
    %   some do not occur in the dataset)
    %
    % #DEFINE-valid
    %
    % Outputs
    % -------
    % RPD: numel(prct) x numel(epochLevels), numel(subEpochLevels)
    % - the repertoire dating percentiles for each percentile/epoch/subepoch 
    % combination
    %
    % uEpochs: vector of values occuring in epochs
    %
    % uSubEpochs: vector of values occuring in subEpochs
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
    % along with this program (see LICENSE.md file).  If not, see <https://www.gnu.org/licenses/>.
    %
    % Repertoire-Dating on Github: <a href="matlab:web('https://github.com/skollmor/repertoireDating', '-browser')">https://github.com/skollmor/repertoireDating</a>
    % Dataspace on Github: <a href="matlab:web('https://github.com/skollmor/dspace', '-browser')">https://github.com/skollmor/dspace</a>

    
    % Parse optional arguments
    S = repertoireDating.internal.parse_optional_arguments(varargin,...
        {'percentiles', 'style', 'uEpochs', 'uSubEpochs', 'valid'},...
        {[5, 25, 50, 75, 95], 'mode a', [], [], []},...
        'repertoireDating.percentiles');
    
    % Check inputs
    if numel(epochs) ~= size(epochs, 1)
        epochs = reshape(epochs, [], 1);
    end
    if ~exist('subEpochs', 'var') || isempty(subEpochs)
        subEpochs = ones(numel(epochs), 1);
    else
        % make sure subEpochs is a column vector
        if numel(subEpochs) ~= size(subEpochs, 1)
            subEpochs = reshape(subEpochs, [], 1);
        end
    end
    if isempty(S.valid)
        S.valid = true(size(NNids, 1), 1);
    end
    if isempty(S.uEpochs)
        S.uEpochs = unique(epochs(S.valid));
    end
    if isempty(S.uSubEpochs)
        S.uSubEpochs = unique(subEpochs(S.valid));
    end
    nDatapoints = size(NNids, 1);
    assert(ismember(S.style, {'mode a', 'mean based'}));
    assert(size(epochs, 1) == nDatapoints);
    assert(size(subEpochs, 1) == nDatapoints);
        
    % Compute percentiles over pooled neighbourhood labels
    RPD = NaN(numel(S.percentiles), numel(S.uEpochs), numel(S.uSubEpochs));
    for epochId = 1:numel(S.uEpochs)
        valid_epoch = S.valid & epochs == S.uEpochs(epochId);
        if sum(valid_epoch) == 0
            continue;
        end
        for subEpochId = 1:numel(S.uSubEpochs)
            valid_subEpoch = valid_epoch & subEpochs == S.uSubEpochs(subEpochId); 
            if sum(valid_subEpoch) == 0
                continue;
            end
            % assemble the epochs (production times) for all neighbors of points 
            % within the current (querry) epoch and subEpoch
            nnEpochs = epochs(NNids(valid_subEpoch, :));       
            switch S.style
                case 'mode a'
                    RPD(:, epochId, subEpochId) = prctile(nnEpochs(:), S.percentiles);
                case 'mean based'
                    RPD(:, epochId, subEpochId) = prctile(nanmean(nnEpochs, 2), S.percentiles);
            end        
        end 
    end
    uEpochs = S.uEpochs;
    uSubEpochs = S.uSubEpochs; 
end

