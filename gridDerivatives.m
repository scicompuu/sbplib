% Calculates derivatives in all directions of function.
% Uses multi-dim vector form.
%   gridDerivatives(F,Dx,Dy)
%   gridDerivatives(F,Dx,Dy,Dz)
function varargout = gridDerivatives(F, varargin)
    assert(length(varargin) == ndims(F));

    switch ndims(F)
        case 2
            varargout{1} = (varargin{1}*F')';
            varargout{2} = varargin{2}*F;
        otherwise
            error('Not implemented for ndims(F) = %d',ndims(F));
    end
end