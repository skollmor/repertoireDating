function S = parse_optional_arguments(args, argNames, argDefaults, identifier)
    %PARSE_OPTIONAL_ARGUMENTS args should be a cell array with the structure 
    % {'valName1', val1, 'valName2', val2}.
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

    
    S = [];
    for k = 1:numel(argNames)
        S.(argNames{k}) = argDefaults{k};
    end
    
    for k = 1:2:numel(args)-1
        if ~isstring(args{k}) && ~ischar(args{k})
            fprintf_orange(['%s: Optional arguments need to be specified in the form '...
            '..., ''argName'', argValue, ''argName2'', argValue2.\nThe name for '...
            'optional argument #%i is not a string:'], identifier, k);
            %disp(args{k});
            fprintf('\nThis argument will be ignored.\n');
            continue;
        end
        
        if ~ismember(args{k}, argNames)
            fprintf_orange(['%s: Optional argument name %i (%s) is not among the allowed argument '...
            'names.\nAllowed argument names are: [%s\b].\nThis argument will be ignored.'],...
            identifier, k, args{k}, cell2mat(cellfun(@(s) [s ' '], argNames, 'uni', false)));
            continue;
        end
       
        S.(args{k}) = args{k+1};
    end
   
end

function fprintf_orange(varargin)
    fprintf(['[\b' varargin{1} ']\b'], varargin{2:end});
end