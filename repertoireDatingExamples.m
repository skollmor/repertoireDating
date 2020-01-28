    
    % This script creates random example data and computes/plots 
    % repertoire dating statistics and mixing matrices.
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


%% Create random data
nEpochs = 40;                                          % number of epochs (e.g. days)
nSubEpochs = 5;                                        % number of sub epochs  (e.g. periods in day)
nSamples = 250;                                        % samples per subEpoch
dim = 25;                                              % data dimension
N = nEpochs * nSubEpochs * nSamples;                   % total number of datapoints
productionTime = (1:N)';                               % production time for each datapoint
epoch = ceil(productionTime/(nSubEpochs*nSamples));
subEpoch = floor(mod(productionTime-1, (nSubEpochs*nSamples))/nSamples) + 1;
data = randn(N, dim) + (1:N)'/N * 2;

%% Compute k-NN graph
NNids = knnsearch(data, data, 'K', 100);               % or use your preferred nearest neighbor searcher 
NNids = NNids(:, 2:end);                               % points should not be their own neighbors

%% Compute and plot repertoire dating pecentiles 
[RPD, RPD_epoch, RPD_subEpoch] = repertoireDating.percentiles(NNids, epoch, subEpoch);
repertoireDating.plotPercentiles(RPD, RPD_epoch, RPD_subEpoch, 16:25);

%% Compute and plot rendition percentiles for all datapoints from epoch 20
repertoireDating.renditionPercentiles(NNids, epoch, 'valid', epoch == 20, 'doPlot', true);
    
%% Compute mixing matrix and plot it
MM = repertoireDating.mixingMatrix(NNids, epoch, 'doPlot', true);

%% Compute stratified mixing matrix and plot it
RP = repertoireDating.renditionPercentiles(NNids, epoch, 'percentiles', 50);
% RP is the repertoire time (50th percentile)
stratMM = repertoireDating.stratifiedMixingMatrices(data, epoch, subEpoch, RP, 'uEpochs', 16:25);
% Avg the stratified mixing matrices
C = arrayfun(@(i) stratMM.allMMs{i}.log2CountRatio, 1:numel(stratMM.allMMs), 'uni', false);
avgStratMM = nanmean(cat(3, C{:}), 3);
% Visualize through MDS
repertoireDating.visualizeStratifiedMixingMatrix(avgStratMM, stratMM);

