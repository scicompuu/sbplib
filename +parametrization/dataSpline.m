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

	% Create spline interpolant
	f = parametrization.Curve.spline(t_data, f_data);

	% Reparametrize with a parameter s in [0, 1]
	tmin = min(t_data);
	tmax = max(t_data);
	t = @(s) tmin + s*(tmax-tmin);

	% Create parameterized curve
	g = @(s) [t(s); f(t(s))];

	% Compute numerical derivative of curve using twice as many points as in data set
	m = 2*m_data;
	ops = sbp.D2Standard(m, {0, 1}, 6);
	gp = parametrization.Curve.numericalDerivative(g, ops.D1);

	% Create curve object
	C = parametrization.Curve(g, gp);

	% Reparametrize with arclength parametrization
	C = C.arcLengthParametrization(m_data);

	% To avoid nested function calls, evaluate curve and compute final spline.
	tv = linspace(0, 1, m_data);
	gv = C.g(tv);
	g = parametrization.Curve.spline(tv, gv);
	C = parametrization.Curve(g);

end
