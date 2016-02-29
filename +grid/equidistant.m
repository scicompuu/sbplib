% Creates a cartesian grid of dimension length(m).
% over the doman xlim, ylim, ...
% Examples:
%   g = grid.equidistant([mx, my], xlim, ylim)
%   g = grid.equidistant([10, 15], {0,1}, {0,2})
function g = equidistant(m, varargin)
    if length(m) ~= length(varargin)
        error('grid:equidistant:NonMatchingParameters','The number of provided dimensions do not match.')
    end

    for i = 1:length(m)
        if ~iscell(varargin{i}) || numel(varargin{i}) ~= 2
           error('grid:equidistant:InvalidLimits','The limits should be cell arrays with 2 elements.');
        end

        if varargin{i}{1} > varargin{i}{2}
            error('grid:equidistant:InvalidLimits','The elements of the limit must be increasing.');
        end
    end

    X = {};
    h = [];
    for i = 1:length(m)
        [X{i}, h(i)] = util.get_grid(varargin{i}{:},m(i));
    end

    g = grid.Cartesian(X{:});
    g.h = h;
end