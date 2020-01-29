function stratMM = stratifiedMixingMatrices(data, epochs, subEpochs, repertoireTimes, varargin)
    % STRATIFIEDMIXINGMATRIX This function computes an average stratified mixing 
    % matrix based on the given k-NN graph and labels. 
    %
    % Inputs
    % ------ 
    % data: N x D
    % - data matrix (rows are datapoints of dimension D)
    % - note: no dimensionality reduction is performed in this function (k-NN graphs
    %   are computed directly on the D dimensional data)
    %
    % epochs: N x 1 
    % 
    % subEpochs: N x 1
    %
    % repertoireTimes: N x 1
    %
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    % nStrata: (default: 5)
    % valid: N x 1, logical, (default: true(N, 1))
    % uEpochs: vector of unique epochs to consider
    % uSubEpochs: vector of unique sub-epochs to consider
    % nConsecutiveEpochs: 1 x 1, integer (default: 2)
    %
    % Output
    % -------
    % stratMM.allMMs                              Cell array that contains the outputs of 
    %                                             repertoireDating.mixingMatrix for each stratMM
    % stratMM.allMMs{j}.NH 
    % stratMM.allMMs{j}.counts                    Nn x Nn, Pooled neighbourhood labels
    % stratMM.allMMs{j}.countRatio                Nn x Nn, counts./NH
    % stratMM.allMMs{j}.log2CountRatio            Nn x Nn, log2(counts./NH)
    % stratMM.allMMs{j}.levels                    1 x Nn, Values corresponding to MM columns
    % stratMM.allMMs{j}.neighbourhoodLevels       1 x Nn, Values corresponding to MM rows
    % stratMM.allMMs{j}.levelSizes
    % stratMM.allMMs{j}.kNeighbours               1 x 1, Number of nearest neighbours
    % stratMM.allMMs{j}.epochs                    Epochs used for this strat-MM
    %
    % stratMM.levels: vector of the form 1:Ns
    % stratMM.strata                              assigns meaning to the rows and columns 
    %                                             of the MMs
    % stratMM.epochs                              -"- 
    % stratMM.subEpochs                           -"-
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
    % % Compute and plot mixing matrix 
    % RP = repertoireDating.renditionPercentiles(NNids, epoch, 'percentiles', 50);
    % % RP is the repertoire time (50th percentile)
    % stratMM = repertoireDating.stratifiedMixingMatrices(data, epoch, subEpoch, RP, 'uEpochs', 16:25);
    % % Avg the stratified mixing matrices
    % C = arrayfun(@(i) stratMM.allMMs{i}.log2CountRatio, 1:numel(stratMM.allMMs), 'uni', false);
    % avgStratMM = nanmean(cat(3, C{:}), 3);
    % % Visualize through MDS
    % repertoireDating.visualizeStratifiedMixingMatrix(avgStratMM, stratMM);
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
        {'nStrata', 'valid', 'uEpochs', 'uSubEpochs', 'nConsecutiveEpochs'},...
        {5, [], [], [], 2},...
        'repertoireDating.stratifiedMixingMatrices');
        
    % Check inputs
    epochs = reshape(epochs, [], 1);
    if ~exist('subEpochs', 'var') || isempty(subEpochs)
        subEpochs = ones(numel(epochs), 1);
    else
        % make sure subEpochs is a column vector
        subEpochs = reshape(subEpochs, [], 1);
    end
    if isempty(S.valid)
        S.valid = true(numel(epochs), 1);
    end
    if isempty(S.uEpochs)
        S.uEpochs = unique(epochs(S.valid));
    end
    if isempty(S.uSubEpochs)
        S.uSubEpochs = unique(subEpochs(S.valid));
    end
    
    % Stratify within each (epoch, subEpoch) pair
    stratifiedRepertoireTimes = NaN(numel(S.valid), 1);
    for ue = 1:numel(S.uEpochs)
        for use = 1:numel(S.uSubEpochs)
            valid = S.valid & epochs == S.uEpochs(ue) & subEpochs == S.uSubEpochs(use);
            [~, idx] = sort(repertoireTimes(valid), 1);
            
            ranks = NaN(numel(idx), 1);
            ranks(idx) = 1:numel(idx);
           
            stratifiedRepertoireTimes(valid) = ceil(S.nStrata * ranks / length(ranks));
        end
    end
    
    % Slice data for parallel processing
    jobs = cell(1, numel(S.uEpochs)-S.nConsecutiveEpochs+1);
    for ule = 1:numel(S.uEpochs)-S.nConsecutiveEpochs+1
        if unique(diff(S.uEpochs(ule:ule+S.nConsecutiveEpochs-1))) ~= 1
            continue;
        end
        valid = S.valid & ismember(epochs, S.uEpochs(ule:ule+S.nConsecutiveEpochs-1));
            
        jobs{ule}.data = data(valid, :);
        jobs{ule}.uEpochs = S.uEpochs(ule:ule+S.nConsecutiveEpochs-1);
        
        rt_ = stratifiedRepertoireTimes(valid, :); 
        [~, ~, epochs_] = unique(epochs(valid, :));
        subEpochs_ = subEpochs(valid, :);        
        
        % Create combined label
        jobs{ule}.labels = (rt_ - 1)*(S.nConsecutiveEpochs*numel(S.uSubEpochs)) +...
            (epochs_ - 1)*numel(S.uSubEpochs) + subEpochs_;
        jobs{ule}.strata = rt_;
        jobs{ule}.epochs = epochs_;
        jobs{ule}.subEpochs = subEpochs_;
    end
    levels = 1:S.nStrata*S.nConsecutiveEpochs*numel(S.uSubEpochs);
    
    allMMs = cell(1, numel(jobs));
    % Compute mixing matrices
    % parfor
    for ule = 1:numel(jobs)
        NNids_ = knnsearch(jobs{ule}.data, jobs{ule}.data, 'K', 50);     % or use your preferred nearest neighbor searcher 
        NNids_ = NNids_(:, 2:end);                                       % points should not be their own neighbors
        
        allMMs{ule} = repertoireDating.mixingMatrix(NNids_, jobs{ule}.labels, 'doPlot', false, 'levels', levels); 
        allMMs{ule}.uEpochs = jobs{ule}.uEpochs;
    end
    
    stratMM.levels = levels;
    stratMM.strata = ceil(levels/(numel(S.uSubEpochs)*S.nConsecutiveEpochs));
    stratMM.epochs = repmat(ceil( (1:numel(S.uSubEpochs)*S.nConsecutiveEpochs)/numel(S.uSubEpochs) ), 1, S.nStrata);
    stratMM.subEpochs = repmat(reshape(S.uSubEpochs, 1, []), 1, S.nStrata*S.nConsecutiveEpochs); 
    stratMM.allMMs = allMMs;
end

