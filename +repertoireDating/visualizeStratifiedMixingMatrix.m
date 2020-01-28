function [] = visualizeStratifiedMixingMatrix(avgStratMM, stratMM, varargin)
    % VISUALIZESTRATIFIEDMIXINGMATRIX This function visualizes a stratified mixing 
    % matrix using MDS 
    %
    % Inputs
    % ------ 
    % stratMM: structure returned by repertoireDating.stratifiedMixingMatrices
    % avgStratMM: (avg) stratification matrix to visualize
    %
    % Optional Inputs
    % ---------------
    % These inputs can be omitted or provided as name, value pairs.
    % 
    % mdsDimension: scalar
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
        {'mdsDimension'},...
        {10},...
        'repertoireDating.visualizeStratifiedMixingMatrix');
        
    sz = unique(size(avgStratMM));
    assert(numel(sz) == 1);
    
    epochs = stratMM.epochs;
    subEpochs = stratMM.subEpochs;
    strata = stratMM.strata;
    
    nStrata = numel(unique(strata));
    nSubEpochs = numel(unique(subEpochs));
    nEpochs = numel(unique(epochs));
    
    dissimilarity = -avgStratMM;
    dissimilarity = dissimilarity - min(dissimilarity(:));
    for k = 1:size(dissimilarity, 1)
        dissimilarity(k, k) = 0;
    end
    dissimilarity(isinf(dissimilarity(:))) = nanmax(dissimilarity(~isinf(dissimilarity(:))));
    dissimilarity =  dissimilarity +  dissimilarity';
    
    opts = statset('MaxIter', 1000);

    [Y, stress, disparities] = mdscale(dissimilarity, S.mdsDimension, 'criterion', 'sstress', 'Options', opts);
    
    Yproj = NaN(size(Y, 1), 4);
    for k = 1:size(Y, 1)
        Y_loo = Y;
        Y_loo(k, :) = NaN;
        basis = computeBasis(Y_loo);
        Yproj(k, :) = Y(k, :) * basis;
    end
    
    pfcn = repertoireDating.internal.bda_figure('Stratified Mix. Mat.', [2, 2], 1);
    pfcn(1, 1, 'a', 1, 1);
    imagesc(avgStratMM); axis xy;
    colormap(parula(1000));
    xlabel('Stratum x Epoch x Subepoch');
    ylabel('Stratum x Epoch x Subepoch');
    cbar =  colorbar('Position', [0.4, 0.75, 0.01, 0.1]);
    ylabel(cbar, 'avg. log2(#H/#NH)', 'interpreter', 'none');
    clims = get(gca, 'clim');
    set(gca, 'clim', max(abs(clims))*[-1, 1]);
    axis tight;
    pfcn(2, 1, 'b', 1, 1);
    scatter(Yproj(:, 1), Yproj(:, 2), 20, epochs + subEpochs/nSubEpochs/2); 
    connectPts(Yproj(:, 1), Yproj(:, 2));
    daspect([1, 1, 1]);
    xlabel('Slow 1');
    ylabel('Slow 2');
    pfcn(1, 2, 'c', 1, 1);
    scatter(Yproj(:, 1), Yproj(:, 3), 20, epochs + subEpochs/nSubEpochs/2);
    connectPts(Yproj(:, 1), Yproj(:, 3));
    daspect([1, 1, 1]);
    xlabel('Slow 1');
    ylabel('Within Epoch');
    pfcn(2, 2, 'd', 1, 1);
    scatter(Yproj(:, 1), Yproj(:, 4), 20, epochs + subEpochs/nSubEpochs/2);
    connectPts(Yproj(:, 1), Yproj(:, 4));
    daspect([1, 1, 1]);
    xlabel('Slow 1');
    ylabel('Across Epoch');
    repertoireDating.internal.bda_formatFigure(gcf, 1, 8);

    function basis = computeBasis(Y)
        % Compute Low-D directions
        Ystrata = cell2mat(arrayfun(@(i) nanmean(Y(strata==i, :), 1), (1:nStrata)', 'uni', false));
        [pcs, ~, ~] = pca(Ystrata);
        slow1 = pcs(:, 1); % pcs are ordered by variance
        
        if Ystrata(nStrata, :) * slow1 < Ystrata(1, :) * slow1
            slow1 = -slow1;
        end
        
        slow2 = pcs(:, 2);
        % 'Within day'
        YsubEpochs = cell2mat(arrayfun(@(i) nanmean(Y(subEpochs==i, :), 1), (1:nSubEpochs)', 'uni', false));
        withinEpoch = (YsubEpochs(end, :) - YsubEpochs(1, :))';
        % 'Across day'
        Yepochs = cell2mat(arrayfun(@(i) nanmean(Y(epochs==i, :), 1), (1:nEpochs)', 'uni', false));
        acrossEpochs = (Yepochs(end, :) - Yepochs(1, :))';
        
        %basis = orthogonalize([slow1 slow2 withinEpoch acrossEpochs]);
        
        % This is used in the paper:
        basis = orthogonalize([slow1 acrossEpochs withinEpoch slow2]);
        basis = basis(:, [1 4 3 2]);
    end
    
    function connectPts(x, y)
        for k_ = 1:nStrata
            for j_ = 1:nEpochs
                valid = strata == k_ & epochs == j_;
                if j_ == 1
                    plot(x(valid), y(valid), '-b');
                else
                    plot(x(valid), y(valid), '-r');
                end
            end
        end
    end
    
end

function v = orthogonalize(v)
    % orthogonalizes column vectors from left to right according to the GramÿSchmidt process
    k = size(v,2);
    assert(k>=2,'The input matrix must have more than one column.');
    for i = 1:1:k
        v(:, i) = v(:, i) / norm(v(:, i));
        for j = i+1:1:k
            v(:, j) = v(:, j) - proj(v(:, i),v(:, j));
        end
    end
end

function w = proj(u, v)
    % projects vector v on vector u
    w = (dot(v, u) / dot(u, u)) * u;
end


