% Creates a curvilinear grid of dimension length(m).
% over the logical domain xi_lim, eta_lim, ...
% If all limits are ommited they are set to {0,1}.
% Examples:
%   g = grid.equidistantCurvilinear(mapping, [m_xi, m_eta])
%   g = grid.equidistantCurvilinear(mapping, [m_xi, m_eta], xi_lim, eta_lim)
%   g = grid.equidistantCurvilinear(mapping, [10, 15], {0,1}, {0,1})
function g = equidistantCurvilinear(mapping, m, varargin)
    if isempty(varargin)
        varargin = repmat({{0,1}}, [1 length(m)]);
    end

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

    g = grid.Curvilinear(mapping, X{:});
    g.logic.h = h;
end