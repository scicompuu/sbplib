classdef Curve
    properties
        g
        gp
        transformation
    end

    methods
        %TODO:
        % Concatenation of curves
        % Subsections of curves
        % Stretching of curve parameter - done for arc length.
        % Curve to cell array of linesegments

        % Returns a curve object.
        %   g -- curve parametrization for parameter between 0 and 1
        %  gp -- parametrization of curve derivative
        function obj = Curve(g,gp,transformation)
            default_arg('gp',[]);
            default_arg('transformation',[]);
            p_test = g(0);
            assert(all(size(p_test) == [2,1]), 'A curve parametrization must return a 2x1 vector.');

            if ~isempty(transformation)
                transformation.base_g = g;
                transformation.base_gp = gp;
                [g,gp] = parametrization.Curve.transform_g(g,gp,transformation);
            end

            obj.g = g;
            obj.gp = gp;
            obj.transformation = transformation;

        end

        function n = normal(obj,t)
            assert(~isempty(obj.gp),'Curve has no derivative!');
            deriv = obj.gp(t);
            normalization = sqrt(sum(deriv.^2,1));
            n = [-deriv(2,:)./normalization; deriv(1,:)./normalization];
        end


        % Plots a curve g(t) for 0<t<1, using n points. Returns a handle h to the plotted curve.
        %   h = plot_curve(g,n)
        function h = plot(obj,n)
            default_arg('n',100);

            t = linspace(0,1,n);
            p = obj.g(t);
            h = line(p(1,:),p(2,:));
        end

        function h= plot_normals(obj,l,n,m)
            default_arg('l',0.1);
            default_arg('n',10);
            default_arg('m',100);
            t_n = linspace(0,1,n);

            normals = obj.normal(t_n)*l;

            n0 = obj.g(t_n);
            n1 = n0 + normals;

            h = line([n0(1,:); n1(1,:)],[n0(2,:); n1(2,:)]);
            set(h,'Color',Color.red);
            obj.plot(m);
        end

        function h= show(obj,name)
            default_arg('name', '')
            p = obj.g(1/2);
            n = obj.normal(1/2);
            p = p + n*0.1;

            % Add arrow


            if ~isempty(name)
                h = text(p(1),p(2),name);
                h.HorizontalAlignment = 'center';
                h.VerticalAlignment = 'middle';
            end

            h = obj.plot();
        end
            % Shows curve with name and arrow for direction.

        % Length of arc from parameter t0 to t1 (which may be vectors).
        % Computed using derivative.
        function L = arcLength(obj,t0,t1)
            assert(~isempty(obj.gp),'Curve has no derivative!');
            speed = @(t) sqrt(sum(obj.gp(t).^2,1));
            L =  integral_vec(speed,t0,t1);
        end

        % Creates the arc length parameterization of a curve.
        %    N -- number of points used to approximate the arclength function
        function curve = arcLengthParametrization(obj,N)
            default_arg('N',100);
            assert(~isempty(obj.gp),'Curve has no derivative!');

            % Construct arcLength function using splines
            tvec = linspace(0,1,N);
            arcVec = obj.arcLength(0,tvec);

            % t as a function of arcLength. Monotonicity-preserving cubic splines.
            tFunc = @(arcLen) pchip(arcVec,tvec,arcLen);

            L = obj.arcLength(0,1);
            arcPar = @(s) tFunc(s*L);

            % New function and derivative
            g_new = @(t)obj.g(arcPar(t));
            gp_new = @(t) normalize(obj.gp(arcPar(t)));
            curve = parametrization.Curve(g_new,gp_new);

        end

        % how to make it work for methods without returns
        function p = subsref(obj,S)
            %Should i add error checking here?
            %Maybe if you want performance you fetch obj.g and then use that
            switch S(1).type
                case '()'
                    p = obj.g(S.subs{1});
                % case '.'

                    % p = obj.(S.subs);
                otherwise
                    p = builtin('subsref',obj,S);
                    % error()
            end
        end


        %% TRANSFORMATION OF A CURVE
        function D = reverse(C)
            % g = C.g;
            % gp = C.gp;
            % D = parametrization.Curve(@(t)g(1-t),@(t)-gp(1-t));
            D = C.transform([],[],-1);
        end

        function D = transform(C,A,b,flip)
            default_arg('A',[1 0; 0 1]);
            default_arg('b',[0; 0]);
            default_arg('flip',1);
            if isempty(C.transformation)
                g  = C.g;
                gp = C.gp;
                transformation.A = A;
                transformation.b = b;
                transformation.flip = flip;
            else
                g  = C.transformation.base_g;
                gp = C.transformation.base_gp;
                A_old = C.transformation.A;
                b_old = C.transformation.b;
                flip_old = C.transformation.flip;

                transformation.A = A*A_old;
                transformation.b = A*b_old + b;
                transformation.flip = flip*flip_old;
            end

            D = parametrization.Curve(g,gp,transformation);

        end

        function D = translate(C,a)
            g = C.g;
            gp = C.gp;

            % function v = g_fun(t)
            %     x = g(t);
            %     v(1,:) = x(1,:)+a(1);
            %     v(2,:) = x(2,:)+a(2);
            % end

            % D = parametrization.Curve(@g_fun,gp);

            D = C.transform([],a);
        end

        function D = mirror(C, a, b)
            assertSize(a,[2,1]);
            assertSize(b,[2,1]);

            g = C.g;
            gp = C.gp;

            l = b-a;
            lx = l(1);
            ly = l(2);


            % fprintf('Singular?\n')

            A = [lx^2-ly^2 2*lx*ly; 2*lx*ly ly^2-lx^2]/(l'*l);

            % function v = g_fun(t)
            %     % v = a + A*(g(t)-a)
            %     x = g(t);

            %     ax1 = x(1,:)-a(1);
            %     ax2 = x(2,:)-a(2);
            %     v(1,:) = a(1)+A(1,:)*[ax1;ax2];
            %     v(2,:) = a(2)+A(2,:)*[ax1;ax2];
            % end

            % function v = gp_fun(t)
            %     v = A*gp(t);
            % end

            % D = parametrization.Curve(@g_fun,@gp_fun);

            % g = A(g-a)+a = Ag - Aa + a;
            b = - A*a + a;
            D = C.transform(A,b);

        end

        function D = rotate(C,a,rad)
            assertSize(a, [2,1]);
            assertSize(rad, [1,1]);
            g = C.g;
            gp = C.gp;


            A = [cos(rad) -sin(rad); sin(rad) cos(rad)];

            % function v = g_fun(t)
            %     % v = a + A*(g(t)-a)
            %     x = g(t);

            %     ax1 = x(1,:)-a(1);
            %     ax2 = x(2,:)-a(2);
            %     v(1,:) = a(1)+A(1,:)*[ax1;ax2];
            %     v(2,:) = a(2)+A(2,:)*[ax1;ax2];
            % end

            % function v = gp_fun(t)
            %     v = A*gp(t);
            % end

            % D = parametrization.Curve(@g_fun,@gp_fun);


             % g = A(g-a)+a = Ag - Aa + a;
            b = - A*a + a;
            D = C.transform(A,b);
        end
    end

    methods (Static)

        % Computes the derivative of g: R -> R^2 using an operator D1
        function gp_out = numericalDerivative(g,D1)
            m = length(D1);
            t = linspace(0,1,m);
            gVec = g(t)';
            gpVec = (D1*gVec)';

            gp1_fun = spline(t,gpVec(1,:));
            gp2_fun = spline(t,gpVec(2,:));
            gp_out = @(t) [gp1_fun(t);gp2_fun(t)];
        end

        function obj = line(p1, p2)

            function v = g_fun(t)
                v(1,:) = p1(1) + t.*(p2(1)-p1(1));
                v(2,:) = p1(2) + t.*(p2(2)-p1(2));
            end

            function v = g_fun_deriv(t)
                v(1,:) = t.*0 + (p2(1)-p1(1));
                v(2,:) = t.*0 + (p2(2)-p1(2));
            end

            obj = parametrization.Curve(@g_fun, @g_fun_deriv);
        end

        function obj = circle(c,r,phi)
            default_arg('phi',[0; 2*pi])
            default_arg('c',[0; 0])
            default_arg('r',1)

            function v = g_fun(t)
                w = phi(1)+t*(phi(2)-phi(1));
                v(1,:) = c(1) + r*cos(w);
                v(2,:) = c(2) + r*sin(w);
            end

            function v = g_fun_deriv(t)
                w = phi(1)+t*(phi(2)-phi(1));
                v(1,:) = -(phi(2)-phi(1))*r*sin(w);
                v(2,:) =  (phi(2)-phi(1))*r*cos(w);
            end

            obj = parametrization.Curve(@g_fun,@g_fun_deriv);
        end

        function obj = bezier(p0, p1, p2, p3)
            function v = g_fun(t)
                v(1,:) = (1-t).^3*p0(1) + 3*(1-t).^2.*t*p1(1) + 3*(1-t).*t.^2*p2(1) + t.^3*p3(1);
                v(2,:) = (1-t).^3*p0(2) + 3*(1-t).^2.*t*p1(2) + 3*(1-t).*t.^2*p2(2) + t.^3*p3(2);
            end

            function v = g_fun_deriv(t)
                v(1,:) = 3*(1-t).^2*(p1(1)-p0(1)) + 6*(1-t).*t*(p2(1)-p1(1)) + 3*t.^2*(p3(1)-p2(1));
                v(2,:) = 3*(1-t).^2*(p1(2)-p0(2)) + 6*(1-t).*t*(p2(2)-p1(2)) + 3*t.^2*(p3(2)-p2(2));
            end

            obj = parametrization.Curve(@g_fun,@g_fun_deriv);
        end


        function [g_out,gp_out] = transform_g(g,gp,tr)
            A = tr.A;
            b = tr.b;
            flip = tr.flip;

            function v = g_fun_noflip(t)
                % v = A*g + b
                x = g(t);

                v(1,:) = A(1,:)*x+b(1);
                v(2,:) = A(2,:)*x+b(2);
            end

            function v = g_fun_flip(t)
                % v = A*g + b
                x = g(1-t);

                v(1,:) = A(1,:)*x+b(1);
                v(2,:) = A(2,:)*x+b(2);
            end


            switch flip
                case 1
                    g_out  = @g_fun_noflip;
                    gp_out = @(t)A*gp(t);
                case -1
                    g_out  = @g_fun_flip;
                    gp_out = @(t)-A*gp(1-t);
            end
        end

    end
end

function g_norm = normalize(g0)
    g1 = g0(1,:);
    g2 = g0(2,:);
    normalization = sqrt(sum(g0.^2,1));
    g_norm = [g1./normalization; g2./normalization];
end

function I = integral_vec(f,a,b)
% Wrapper around the built-in function integral that
% handles multiple limits.

    Na = length(a);
    Nb = length(b);
    assert(Na == 1 || Nb == 1 || Na==Nb,...
        'a and b must have same length, unless one is a scalar.');

    if(Na>Nb);
        I = zeros(size(a));
        for i = 1:Na
            I(i) = integral(f,a(i),b);
        end
    elseif(Nb>Na)
        I = zeros(size(b));
        for i = 1:Nb
            I(i) = integral(f,a,b(i));
        end
    else
        I = zeros(size(b));
        for i = 1:Nb
            I(i) = integral(f,a(i),b(i));
        end
    end
end

% Returns a function handle to the spline.
function f = spline(tval,fval,spline_order)
    default_arg('spline_order',4);
    [m,~] = size(tval);
    assert(m==1,'Need row vectors.');

    f_spline = spapi( optknt(tval,spline_order), tval, fval );
    f = @(t) fnval(f_spline,t);
end
