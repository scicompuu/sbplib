% Accepts data points (t_i, f_i) and returns a Curve object,
% using spline interpolation.
% The spline curve is parametrized with the arc length parametrization
% to facilitate better grids.
%
% t_data 	- vector
% f_data 	- vector
% C 	- curve object
function C = dataSpline(t_data, f_data)

	assert(length(t_data)==length(f_data),'Vectors must be same length');
	m_data = length(t_data);

	pp_g = spapi(4, t_data, f_data); % equivalent to g = spapi(aptknt(t_data, 4), t_data, f_data)
	% or  (not sure what the difference is?!)
	% g = spapi(optknt(t_data, 4), t_data, f_data)
	pp_gp = fnder(g);

	g = @(t)fnval(pp_g, t);
	pp_gp = @(t)fnval(pp_gp, t);

	C = parametrization.Curve(g, gp);

	% Reparametrize with arclength parametrization
	C = C.arcLengthParametrization(m_data);
end
