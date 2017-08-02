% Calls the constructor of an object.
% Might be usefull to call the constructor of a subclass object in the superclass
function obj = callConstructor(subclassObj, varargin)
    fun = str2func(class(subclassObj));
    obj = fun(varargin{:});
end
