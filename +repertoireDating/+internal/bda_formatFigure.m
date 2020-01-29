function bda_formatFigure( fh, lineWidth, fontSize, gridstate )
    % BDA_FORMATFIGURE Function to modify figure appearence. 
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



    defaultWidth = 1200;
    defaultFontSize = 13;
    
    if nargin < 2 || isempty(lineWidth)
        lineWidth = 2;
    end
    if nargin < 3 || isempty(fontSize)
        fontSize = 13;
    end
    if nargin < 4 || isempty(gridstate)
        gridstate = 'on';
    end
    defaultFontSize = fontSize;
    
    %     op = get(fh, 'Position');
    %     set(fh, 'Position', [20, 80, defaultWidth, op(4)*defaultWidth/op(3)]);
    set(fh, 'Color', 'w');
    ch = get(fh, 'Children');
    
    for h = ch'
        try
            set(h, 'fontsize', defaultFontSize);
            xl = get(h, 'xlabel');
            cellfun(@(q) set(q, 'fontsize', defaultFontSize), {xl});
            yl = get(h, 'ylabel');
            cellfun(@(q) set(q, 'fontsize', defaultFontSize), {yl});
            tl = get(h, 'title');
            cellfun(@(q) set(q, 'fontsize', defaultFontSize), {tl});
            
            %ls = get(h, 'children');
            for hl = findall(h, 'Type', 'line')
                try
                    set(hl, 'linewidth', lineWidth);
                    %set(hl, 'linesmoothing', 'on');
                catch
                   fprintf('N');
                end
            end
            
            set(h, 'box', 'on');
            set(h, 'xgrid', gridstate);
            set(h, 'ygrid', gridstate);
        catch
        end
        
    end
    
end

