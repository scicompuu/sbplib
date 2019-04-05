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

	pp_g = spapi(4, t_data, f_data);
	pp_gp = fnder(pp_g);

	% Reparametrize with parameter s from 0 to 1 to use Curve class
    tmin = min(t_data);
    tmax = max(t_data);
    t = @(s) tmin + s*(tmax-tmin);
    dt_ds = @(s) 0*s + (tmax-tmin);
	g = @(s) [t(s); fnval(pp_g, t(s))];
	gp = @(s) [dt_ds(s); fnval(pp_gp, t(s)).*dt_ds(s)];

	% Create Curve object and reparametrize with arclength parametrization
	C = parametrization.Curve(g, gp);
	C = C.arcLengthParametrization(m_data);
end
