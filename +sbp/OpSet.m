classdef (Abstract) OpSet
    properties (Abstract)
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
        x % Grid
    end

end