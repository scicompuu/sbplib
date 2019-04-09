% dataSpline calculates a Curve through the points f_i using cubic spline interpolation.
% The spline curve is parametrized with the arc length parametrization
% to facilitate better grids.
%
% f 	- m x D matrix of m points in D dimensions
function C = dataSpline(f)
	m = size(f, 1);

	t = linspace(0,1,m);

	pp_g = spapi(4, t, f');
	pp_gp = fnder(pp_g);

	g  = @(t) fnval(pp_g, t);
	gp = @(t) fnval(pp_gp, t);

	C = parametrization.Curve(g, gp);
	C = C.arcLengthParametrization();
end
