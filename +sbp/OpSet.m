classdef (Abstract) OpSet
    properties (Abstract)
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end

    methods (Abstract)

    end

    methods (Abstract, Static)
        lambda = smallestGrid()
    end
end