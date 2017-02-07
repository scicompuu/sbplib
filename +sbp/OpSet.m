classdef (Abstract) OpSet
    properties (Abstract)
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
        x % Grid
    end

    methods (Abstract)
        % Returns a string representation of the type of operator.
        str = string(obj)
    end

end