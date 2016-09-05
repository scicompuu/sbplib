classdef triang_interp
    properties
        g1, g2 ,g3  % Curves encirling the tirangle in the positive direction.
        A,B,C  % The corners of the triangle
        Sa, Sb, Sc % Mappings from square with different sides collapsed
    end

    methods
        function o = triang_interp(g1,g2,g3)
            o.g1 = g1;
            o.g2 = g2;
            o.g3 = g3;
            o.A = g1(0);
            o.B = g2(0);
            o.C = g3(0);
            o.Sa = parametrization.triang_interp.square_to_triangle_interp(g2,g3,g1);
            o.Sb = parametrization.triang_interp.square_to_triangle_interp(g3,g1,g2);
            o.Sc = parametrization.triang_interp.square_to_triangle_interp(g1,g2,g3);
        end


        function show(o,N)
            % Show the mapped meridians of the triangle.
            % Might be used for the barycentric coordinates.
            ma = @(t)o.Sa(1/2,1-t);
            mb = @(t)o.Sb(1/2,1-t);
            mc = @(t)o.Sc(1/2,1-t);

            na = @(t)o.Sa(t,1/2);

            ka = @(t)(o.g1(1-t)+o.g2(t))/2;

            h = parametrization.plot_curve(ma);
            h.Color = Color.blue;
            h = parametrization.plot_curve(mb);
            h.Color = Color.blue;
            h = parametrization.plot_curve(mc);
            h.Color = Color.blue;

            h = parametrization.plot_curve(na);
            h.Color = Color.red;

            h = parametrization.plot_curve(ka);
            h.Color = Color.red;

            [a(1),a(2)] = ma(1/3);
            [b(1),b(2)] = mb(1/3);
            [c(1),c(2)] = mc(1/3);

            d = ka(1-1/3);


            parametrization.label_pt(a,b,c,d);


            % t = linspace(0,1,N);
            % for i = 1:N
            %     sa = @(s)o.Sa(s,t(i));
            %     sb = @(s)o.Sb(s,t(i));
            %     sc = @(s)o.Sc(s,t(i));

            %     h = parametrization.plot_curve(sa);
            %     h.Color = Color.blue;
            %     h = parametrization.plot_curve(sb);
            %     h.Color = Color.blue;
            %     h = parametrization.plot_curve(sc);
            %     h.Color = Color.blue;
            % end

            h = parametrization.plot_curve(o.g1);
            h.LineWidth = 2;
            h.Color = Color.red;

            h = parametrization.plot_curve(o.g2);
            h.LineWidth = 2;
            h.Color = Color.red;

            h = parametrization.plot_curve(o.g3);
            h.LineWidth = 2;
            h.Color = Color.red;

        end


    end

    methods(Static)
        % Makes a mapping from the unit square to a triangle by collapsing
        % one of the sides of the squares to a corner on the triangle
        % The collapsed side is mapped to the corner oposite to g1.
        % This is done such that for S(s,t), S(s,1) = g1(s)
        function S = square_to_triangle_interp(g1,g2,g3)
            corner = parametrization.line_segment(g3(0),g3(0));
            S = parametrization.transfinite_interp(corner,g3,f(g1),f(g2))

            % Function to flip a curve
            function h = f(g)
                h = @(t)g(1-t);
            end
        end
    end

end

% % Return a mapping from u.v to x,y of the domain encircled by g1 g2 g3 in the the positive direction. created be using transfinite interpolation.
% function S = triang_interp(g1,g2,g3)
%     A = g1(0)
%     B = g2(0)
%     C = g3(0)

%     function [x,y] = S_fun(u,v)
%         w = sqrt((u-1)^2+v^2)/sqrt(2); % Parameter for g3
%         v = v*(1-u-v)*g1(u) + u*(1-u-v)*g2(v) + u*v*g3(w) ...
%             +(1-u)*(1-v)*A+u*(1-v)*B + (1-u)*v*C;
%         x = v(1);
%         y = v(2);
%     end
%     S = @S_fun;
% end



% function subsref(obj,S)
%       if ~all(isnumeric(S.subs{:}))
%         error('Only supports calling object with number')
%       end
%       if numel(S.subs{:}) > 1
%         disp('You''ve called the object with more than one argument');
%       else
%         disp(['You called the object with argument = ',num2str(S.subs{:})]);
%       end
%     end