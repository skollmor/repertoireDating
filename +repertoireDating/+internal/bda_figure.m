function [ panelfcn, figureId ] = bda_figure( name, design, aspectratioYoverX, widthOnPaper,...
        figureId, horizontalSpacing, verticalSpacing)
    % BDA_FIGURE Funtion to simplify the creation of multi-panel figures.
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

    
    defaultWidth = 800;
    
    % Make sure figure size lies within screen dimensions
    maxSize = get(groot, 'Screensize');
    % Adjust figure width < screen width
    defaultWidth = ~(defaultWidth > maxSize(3)) * defaultWidth +...
        (defaultWidth > maxSize(3)) * maxSize(3);
    % Adjust figure height < screen height-taskbar/window border (approx. 110px)
    defaultWidth = ~(defaultWidth*aspectratioYoverX > maxSize(4)-110) * defaultWidth +...
        (defaultWidth*aspectratioYoverX > maxSize(4)-110) * ((maxSize(4)-110)/aspectratioYoverX);
    %disp(['bda_figure window: height: ' num2str(defaultWidth*aspectratioYoverX) ' width: ' num2str(defaultWidth)]);
    %clear maxSize
    
    function h = panel_fcn(x, y, panelLabel, xe, ye)
        if nargin == 2
            panelNumber = x-1;
            panelLabel = y;
            y = floor(panelNumber/design(1))+1;
            x = rem(panelNumber, design(1))+1;
        end
        if nargin < 5
            xe = 1;
            ye = 1;
        end
        figure(figureId);
        h = repertoireDating.internal.subaxis.subaxis(design(2), design(1), x, y, xe, ye, 'SpacingHoriz', horizontalSpacing,...
            'SpacingVert', verticalSpacing);
        
        text(-0.09 * design(1)/xe, 1 + 0.07/ye, panelLabel, 'Units', 'normalized',...
            'Fontsize', 15, 'fontweight', 'bold'); hold all;
    end
    panelfcn = @panel_fcn;

    screenPos = [20, 80, defaultWidth, aspectratioYoverX*defaultWidth];
    if nargin < 7 || isempty(verticalSpacing)
        verticalSpacing = 0.08;
    end
    
    if nargin < 6 || isempty(horizontalSpacing)
        horizontalSpacing = 0.08;
    end
    
    if nargin < 5 || isempty(figureId)
        figureId = figure();
    else
        figureId = figure(figureId);
        clf;
    end
    if nargin < 4 || isempty(widthOnPaper)
        widthOnPaper = 17.35;
    end
    
    
    set(figureId, 'Name', name, 'Position', screenPos); 
    set(figureId, 'PaperUnits', 'centimeters');
    set(figureId, 'PaperPosition', [1, 1, widthOnPaper, aspectratioYoverX*widthOnPaper]);
    
end

