classdef (Abstract) InterpOps
    properties (Abstract)
        Iu2v % Interpolation operator(s) from "u" to "v"
        Iv2u % Interpolation operator(s) from "v" to "u"
    end

    methods (Abstract)
        % Returns a string representation of the type of operator.
        str = string(obj)
    end

end