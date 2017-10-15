classdef (Abstract) InterpOps
    properties (Abstract)
        % C and F may refer to coarse and fine, but it's not a must.
        IC2F % Interpolation operator from "C" to "F"
        IF2C % Interpolation operator from "F" to "C"

    end

    methods (Abstract)
        % Returns a string representation of the type of operator.
        str = string(obj)
    end

end