function MM = mixingMatrix(NNids, labels, varargin)
    % MIXINGMATRIX This function computes a mixing matrix based on the given
    % k-NN graph and labels.
    %
    % Inputs
    % ------ 
    % NNids: N x #Neighbours (integer)
    % - defines the k-NN graph
    % - entires of this matrix of indices must be in 1..N (N is the total number of datapoints)
    % - row j in this matrix contains the indices of the k nearest neighbours of datapoint j.
    % - points should not be their own neighbors
    %
    % labels: N x 1 
    % - defines the labels mixing statistics are based on
    % - unique values are used as bins, no other binning takes place
    % 
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    %   --Fixed Unique Labels--
    %   Sometimes certain labels do not occur in the data but it it convenient to
    %   still include them in the mixing matrix. For these cases, the levels of the 
    %   label variable can be explicitly provided as a vector:
    %
    %   repertoireDating.mixingMatrix(_, 'levels', Lq)
    %
    %   --Different X and Y axes--
    %   The labels defining the x-axis of the mixing matrix and the labels defining the 
    %   y-axis can be indepent. Specify neighbourhoodLabels (y-axis) that differ from
    %   querryLabels (labels, x-axis) by using: 
    %
    %   repertoireDating.mixingMatrix(_, 'neighbourhoodLabels', labels);
    %
    %   Specify fixed bins by using: 
    %
    %   repertoireDating.mixingMatrix(_, 'neighbourhoodLevels', Ln)
    %
    %   --Plotting--
    %   To plot the computed mixing matrix use: 
    %   repertoireDating.mixingMatrix(_, 'doPlot', true);
    %
    % Output
    % -------
    % MM.NH                    Nn x Nq, Counts under a Null Hypothesis of permuted labels               
    % MM.counts                Nn x Nq, Pooled neighbourhood labels
    % MM.countRatio            Nn x Nq, counts./NH
    % MM.log2CountRatio        Nn x Nq, log2(counts./NH)
    %
    % MM.levels                1 x Nq, Values corresponding to MM columns
    % MM.neighbourhoodLevels   1 x Nn, Values corresponding to MM rows
    % MM.levelSizes
    % MM.kNeighbours           1 x 1, Number of nearest neighbours
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
    % MM = repertoireDating.mixingMatrix(NNids, epoch, 'doPlot', true, 'doMDS', true);
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

    
    % Create a struct from the list of optional arguments and their default values
    S = repertoireDating.internal.parse_optional_arguments(varargin,...
        {'levels', 'neighbourhoodLabels', 'neighbourhoodLevels', 'valid', 'doPlot', 'doMDS', 'doPlotMDS', 'minPointsPerLevel'},...
        {[], [], [], [], false, false, true, 10},...
        'repertoireDating.mixingMatrix');
        
    labels = reshape(labels, [], 1);
    if isempty(S.valid)
        S.valid = true(1, numel(labels));
    end
    if isempty(S.levels)
        S.levels = unique(labels);
    end
    if isempty(S.neighbourhoodLabels)
        S.neighbourhoodLabels = labels;
    end
    S.neighbourhoodLabels = reshape(S.neighbourhoodLabels, [], 1);
    if isempty(S.neighbourhoodLevels)
        S.neighbourhoodLevels = unique(S.neighbourhoodLabels);
    end
   
    LGnans = sum(isnan(NNids), 2);
    if sum(LGnans) > 0
        repertoireDating.internal.fprintf_orange('%i Neighborhoods contain NaN ids. (avg #NaNs: %.2f)\n', sum(LGnans > 0),...
            mean(LGnans(LGnans > 0)));
        return;
    end
    
    MM = [];
    [MM.counts, MM.NHcounts, MM.levelSizes, MM.neighbourhoodLevelSizes] =...
        mixingMatrix_outerproductNH(NNids, LGnans, S.levels, S.neighbourhoodLevels,...
        labels, S.neighbourhoodLabels, S.valid, S.minPointsPerLevel);
    
    MM.counts(MM.levelSizes < S.minPointsPerLevel, :) = NaN;
    MM.countRatio = MM.counts./MM.NHcounts;
    MM.log2CountRatio = log2(MM.countRatio);
    
    if S.doPlot
        titleStr = sprintf('Mix. Mat. - %i Pts; %i Nbs; %i Levels', sum(S.valid), size(NNids, 2), numel(S.levels));
        pfcn = repertoireDating.internal.bda_figure(titleStr, [1, 5], 1);
        pfcn(1, 1, '', 1, 4);
        imagesc([min(S.levels), max(S.levels)], [min(S.neighbourhoodLevels),...
            max(S.neighbourhoodLevels)], MM.log2CountRatio); axis xy; axis tight;
        colormap(parula(1000));
        cbar =  colorbar('Position', [0.91, 0.4, 0.025, 0.5]);
        ylabel(cbar, 'log2(#H/#NH)', 'interpreter', 'none');
        clims = get(gca, 'clim');
        set(gca, 'clim', max(abs(clims))*[-1, 1]);
        ax(1) = gca;
        pfcn(1, 5, '', 1, 1);
        title(sprintf('n=%i', sum(MM.levelSizes)));
        bar(S.levels, MM.levelSizes);
        xlim([min(S.levels), max(S.levels)]);
        ax(2) = gca;
        linkaxes(ax, 'x');
        xlabel('Levels', 'interpreter', 'none');
        repertoireDating.internal.bda_formatFigure(gcf, 1, 8);
    end
    
    if S.doMDS
        dissimilarity = -MM.log2CountRatio;
        dissimilarity = dissimilarity - min(dissimilarity(:));
        for k = 1:size(dissimilarity, 1)
            dissimilarity(k, k) = 0;
        end
        dissimilarity(isinf(dissimilarity(:))) = nanmax(dissimilarity(~isinf(dissimilarity(:))));
        dissimilarity =  dissimilarity +  dissimilarity';
        MM.MDS_dissimilarityMatrix = dissimilarity;
        opts = statset('MaxIter', 500);
        [MM.MDS_solution, MM.MDS_stress, MM.MDS_disparities] = mdscale(dissimilarity,...
            2, 'criterion', 'sstress', 'Options', opts);
   
        if S.doPlotMDS
            titleStr = sprintf('MDS of Mixing Matrix');
            pfcn2 = repertoireDating.internal.bda_figure(titleStr, [2, 2], 1);
            pfcn2(1,1,'');
            imagesc(dissimilarity);
            axis tight;
            colorbar('Position', [0.47, 0.8, 0.025, 0.1]);
            pfcn2(1,2,'');
            scatter(MM.MDS_solution(:, 1), MM.MDS_solution(:, 2), 40, S.levels, 'filled');
            plot(MM.MDS_solution(:, 1), MM.MDS_solution(:, 2), ':k');
            xlabel('MDS 1'); ylabel('MDS 2');
            title(sprintf('MDS stress=%f', MM.MDS_stress));
            xlims = get(gca, 'xlim');
            ylims = get(gca, 'ylim');
            ml = [-1, 1] * max(abs([xlims, ylims]));
            xlim(ml);
            ylim(ml);
            pfcn2(2, 1, '');
            ldDist = squareform(pdist(MM.MDS_solution));
            [~, ord] = sortrows([MM.MDS_disparities(:) dissimilarity(:)]);
            plot(dissimilarity(:),ldDist(:),'bo', ...
                dissimilarity(ord),MM.MDS_disparities(ord),'r.-');
            xlabel('HD Dissimilarity')
            ylabel('MDS Distance/Disparity')
            plot(dissimilarity(:), MM.MDS_disparities(:), '+k');
            legend({'MDS Distances' 'Disparity', 'Disparity'}, 'Location','NorthWest');
            pfcn2(2, 2, '');
            plot(MM.log2CountRatio(:), ldDist(:), '.k');
            xlabel('log2(#H/#NH)', 'interpreter', 'none');
            ylabel('MDS Distance');
            repertoireDating.internal.bda_formatFigure(gcf, 1, 8);
        end
    end 
end

function [counts, NHcounts, levelSizes, neighbourhoodLevelSizes] =...
        mixingMatrix_outerproductNH(LGids, LGnans, levels, targetBins, levelvar,...
        targetvar, pf, minLevelSize)
    
    counts = NaN(numel(targetBins), numel(levels));
    NHcounts = NaN(numel(targetBins), numel(levels));
    levelSizes = NaN(numel(levels), 1);
    neighbourhoodLevelSizes = NaN(numel(targetBins), 1);
    lgSize = size(LGids, 2);
    
    for k = 1:numel(targetBins)
        validIds = targetvar(pf) == targetBins(k) & LGnans == 0;
        neighbourhoodLevelSizes(k) = sum(validIds);
    end
    
    for k = 1:numel(levels)
        validIds = levelvar(pf) == levels(k) & LGnans == 0;
        levelSizes(k) = sum(validIds);
        if sum(validIds) < minLevelSize
            counts(:, k) = NaN(numel(targetBins), 1);
            NHcounts(:, k) = NaN(numel(targetBins), 1);
            NHcounts(:, k) = NaN(numel(targetBins), 1);
            continue;
        end
        A = LGids(validIds, :);
        A = targetvar(A);
        A = A(:);
        counts(:, k) = histc(A, targetBins);
        NHcounts(:, k) = lgSize * levelSizes(k) * neighbourhoodLevelSizes/sum(neighbourhoodLevelSizes);   
    end
end