classdef Curve
    properties
        g
        gp
        transformation
        arc_length_fun 
        
        % arc_length_fun: 
        % function handle. arc_length(s)=distance 
        % between t=0 and t=s.
    end
    methods
        %TODO:
        % Errors or FD if there is no derivative function added.
        % -semi-done

        % Concatenation of curves
        % Subsections of curves
        % Stretching of curve paramter - semi-done
        % Curve to cell array of linesegments

        % Should supply either derivative or a difference operator D1
        function obj = Curve(g,gp,transformation,D1,arc_length_fun)
            default_arg('gp',[]);
            default_arg('transformation',[]);
            default_arg('D1',[]);
            default_arg('arc_length_fun',[]);
            p_test = g(0);
            assert(all(size(p_test) == [2,1]), 'A curve parametrization must return a 2x1 vector.');

            if(isempty(gp) && isempty(D1));
                % Should be error instead of disp.
                disp(['You really should supply either the exact derivative ',...
                'or a suitable difference operator to compute an approximation.']);
            elseif(isempty(gp))
                gp = grid.Curve.numerical_derivative(g,D1);
            end
            
            if ~isempty(transformation)
                transformation.base_g = g;
                transformation.base_gp = gp;
                [g,gp] = grid.Curve.transform_g(g,gp,transformation);
            end

            if(isempty(arc_length_fun))
                if(~isempty(D1)); N = length(D1); % Same accuracy as for deriv.
                else N = 101; % Best way to let user choose?
                end
                tvec = linspace(0,1,N);
                arc_vec = grid.Curve.arc_length(gp,0,tvec);
                arc_length_fun = grid.Curve.spline(tvec,arc_vec);
            end
            
            obj.g = g;
            obj.gp = gp;
            obj.transformation = transformation;
            obj.arc_length_fun = arc_length_fun;
            

        end

        % Made up length calculation!! Science this before actual use!!
        % Calculates the length of the curve. Makes sure the longet segment used is shorter than h_max.
        function [L,t] = curve_length(C,h_max)
            default_arg('h_max',0.001);
            g = C.g;
            h = 0.1;
            m = 1/h+1;
            t = linspace(0,1,m);

            [p,d] = get_d(t,g);

            while any(d>h_max)
                I = find(d>h_max);

                % plot(p(1,:),p(2,:),'.')
                % waitforbuttonpress

                new_t = [];
                for i = I
                    new_t(end +1) = (t(i)+t(i+1))/2;
                end
                t = [t new_t];
                t = sort(t);

                [p,d] = get_d(t,g);
            end

            L = sum(d);

            function [p,d] = get_d(t,g)
                n = length(t);

                p = zeros(2,n);
                for i = 1:n
                    p(:,i) = g(t(i));
                end

                d = zeros(1,n-1);
                for i = 1:n-1
                    d(i) = norm(p(:,i) - p(:,i+1));
                end
            end
        end

        function n = normal(obj,t)
            deriv = obj.gp(t);
            normalization = sqrt(sum(deriv.^2,1));
            n = [-deriv(2,:)./normalization; deriv(1,:)./normalization];
        end


        % Plots a curve g(t) for 0<t<1, using n points. Returns a handle h to the plotted curve.
        %   h = plot_curve(g,n)
        function h = plot(obj,n,marker)
            default_arg('n',100);
            default_arg('marker','line')

            t = linspace(0,1,n);

            p = obj.g(t);

            switch marker
                case 'line'
                    h = line(p(1,:),p(2,:));
                otherwise
                    h = plot(p(1,:),p(2,:),marker);
            end
        end
        
        % Plots the derivative gp(t) for 0<t<1, using n points. Returns a handle h to the plotted curve.
        %   h = plot_curve(gp,n)
        function h = plot_derivative(obj,n,marker)
            default_arg('n',100);
            default_arg('marker','line')

            t = linspace(0,1,n);

            p = obj.gp(t);

            switch marker
                case 'line'
                    h = line(p(1,:),p(2,:));
                otherwise
                    h = plot(p(1,:),p(2,:),marker);
            end
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
            p = obj.g(1/2);
            n = obj.normal(1/2);
            p = p + n*0.1;

            % Add arrow

            h = text(p(1),p(2),name);
            h.HorizontalAlignment = 'center';
            h.VerticalAlignment = 'middle';

            obj.plot();
        end
            % Shows curve with name and arrow for direction.

        
        function curve = stretch_parameter(obj,type)
            default_arg('type','arc_length');
            switch type
                % Arc length parameterization.
                case 'arc_length'
                    arcLength = obj.arc_length_fun;
                    arcPar = @(t) util.fzero_vec(@(s)arcLength(s) - t*arcLength(1),[0-10*eps,1+10*eps]);
                    g_new = @(t)obj.g(arcPar(t));
                    gp_old = obj.gp;
                    gp_new = @(t) normalize(gp_old(arcPar(t)));
                    
                    arc_len_new = @(t) t;
                    curve = grid.Curve(g_new,gp_new,[],[],arc_len_new);
                    
                otherwise
                    error('That stretching is not implemented.');
            end
            
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
            % D = grid.Curve(@(t)g(1-t),@(t)-gp(1-t));
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

            D = grid.Curve(g,gp,transformation);

        end

        function D = translate(C,a)
            g = C.g;
            gp = C.gp;

            % function v = g_fun(t)
            %     x = g(t);
            %     v(1,:) = x(1,:)+a(1);
            %     v(2,:) = x(2,:)+a(2);
            % end

            % D = grid.Curve(@g_fun,gp);

            D = C.transform([],a);
        end

        function D = mirror(C, a, b)
            assert_size(a,[2,1]);
            assert_size(b,[2,1]);

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

            % D = grid.Curve(@g_fun,@gp_fun);

            % g = A(g-a)+a = Ag - Aa + a;
            b = - A*a + a;
            D = C.transform(A,b);

        end

        function D = rotate(C,a,rad)
            assert_size(a, [2,1]);
            assert_size(rad, [1,1]);
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

            % D = grid.Curve(@g_fun,@gp_fun);


             % g = A(g-a)+a = Ag - Aa + a;
            b = - A*a + a;
            D = C.transform(A,b);
        end
    end

    methods (Static)
        
        % Length of arc from parameter t0 to t1 (which may be vectors).
        % Computed using derivative.
        function L = arc_length(deriv,t0,t1)
            speed = @(t) sp(deriv(t));
            
            function s = sp(deriv)
                s = sqrt(sum(deriv.^2,1));
            end
            L =  util.integral_vec(speed,t0,t1);
        end
        
        function gp_out = numerical_derivative(g,D1)
            m = length(D1); L = 1; % Assume curve parameter from 0 to 1.
            t = linspace(0,L,m); 
            g = g(t)';
            gp = (D1*g)';

            gp1_fun = grid.Curve.spline(t,gp(1,:));
            gp2_fun = grid.Curve.spline(t,gp(2,:));
            gp_out = @(t) [gp1_fun(t);gp2_fun(t)];
            
        end
        
        % Returns a function handle to the spline.
        function f = spline(tval,fval,spline_order)
            default_arg('spline_order',4);
            [m,~] = size(tval);
            assert(m==1,'Need row vectors.');
            
            % make vectors longer to be safe slightly beyond edges.
            dt0 = tval(2)-tval(1); dt1 = tval(end)-tval(end-1);
            df0 = fval(2)-fval(1); df1 = fval(end)-fval(end-1);
            tval = [tval(1)-dt0, tval, tval(end)+dt1];
            fval = [fval(1)-df0, fval, fval(end)+df1];
            
            f_spline = spapi( optknt(tval,spline_order), tval, fval );
            f = @(t) fnval(f_spline,t);
        end
        
        function obj = line(p1, p2)

            function v = g_fun(t)
                v(1,:) = p1(1) + t.*(p2(1)-p1(1));
                v(2,:) = p1(2) + t.*(p2(2)-p1(2));
            end
            g = @g_fun;

            obj = grid.Curve(g);
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

            obj = grid.Curve(@g_fun,@g_fun_deriv);
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

            obj = grid.Curve(@g_fun,@g_fun_deriv);
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
    g1 = g0(1,:); g2 = g0(2,:);
    normalization = sqrt(sum(g0.^2,1));
    g_norm = [g1./normalization; g2./normalization];
end

    
    
