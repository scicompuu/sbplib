classdef Ti
    properties
        gs % {4}Curve
        S  % FunctionHandle(u,v)
    end

    methods
        % TODO function to label boundary names.
        %  function to find largest and smallest delta h in the grid. Maybe shouldnt live here
        function obj = Ti(C1,C2,C3,C4)
            obj.gs = {C1,C2,C3,C4};

            g1 = C1.g;
            g2 = C2.g;
            g3 = C3.g;
            g4 = C4.g;

            A = g1(0);
            B = g2(0);
            C = g3(0);
            D = g4(0);

            function o = S_fun(u,v)
                x1 = g1(u);
                x2 = g2(v);
                x3 = g3(1-u);
                x4 = g4(1-v);
                o1 = (1-v).*x1(1,:) + u.*x2(1,:) + v.*x3(1,:) + (1-u).*x4(1,:) ...
                    -((1-u)*(1-v).*A(1,:) + u*(1-v).*B(1,:) + u*v.*C(1,:) + (1-u)*v.*D(1,:));
                o2 = (1-v).*x1(2,:) + u.*x2(2,:) + v.*x3(2,:) + (1-u).*x4(2,:) ...
                    -((1-u)*(1-v).*A(2,:) + u*(1-v).*B(2,:) + u*v.*C(2,:) + (1-u)*v.*D(2,:));

                o = [o1;o2];
            end

            obj.S = @S_fun;
        end

        function [X,Y] = map(obj,u,v)
            default_arg('v',u);

            if isscalar(u)
                u = linspace(0,1,u);
            end

            if isscalar(v)
                v = linspace(0,1,v);
            end

            S = obj.S;

            nu = length(u);
            nv = length(v);

            X = zeros(nv,nu);
            Y = zeros(nv,nu);

            u = rowVector(u);
            v = rowVector(v);

            for i = 1:nv
                p = S(u,v(i));
                X(i,:) = p(1,:);
                Y(i,:) = p(2,:);
            end
        end

        function h = plot(obj,nu,nv)
            S = obj.S;

            default_arg('nv',nu)

            u = linspace(0,1,nu);
            v = linspace(0,1,nv);

            m = 100;

            X = zeros(nu+nv,m);
            Y = zeros(nu+nv,m);


            t = linspace(0,1,m);
            for i = 1:nu
                p = S(u(i),t);
                X(i,:) = p(1,:);
                Y(i,:) = p(2,:);
            end

            for i = 1:nv
                p = S(t,v(i));
                X(i+nu,:) = p(1,:);
                Y(i+nu,:) = p(2,:);
            end

            h = line(X',Y');
        end


        function h = show(obj,nu,nv)
            default_arg('nv',nu)
            S = obj.S;

            if(nu>2 || nv>2)
                h_grid = obj.plot(nu,nv);
                set(h_grid,'Color',[0 0.4470 0.7410]);
            end

            h_bord = obj.plot(2,2);
            set(h_bord,'Color',[0.8500 0.3250 0.0980]);
            set(h_bord,'LineWidth',2);
        end


        % TRANSFORMATIONS
        function ti = translate(obj,a)
            gs = obj.gs;

            for i = 1:length(gs)
                new_gs{i} = gs{i}.translate(a);
            end

            ti = grid.Ti(new_gs{:});
        end

        % Mirrors the Ti so that the resulting Ti is still left handed.
        %  (Corrected by reversing curves and switching e and w)
        function ti = mirror(obj, a, b)
            gs = obj.gs;

            new_gs = cell(1,4);

            new_gs{1} = gs{1}.mirror(a,b).reverse();
            new_gs{3} = gs{3}.mirror(a,b).reverse();
            new_gs{2} = gs{4}.mirror(a,b).reverse();
            new_gs{4} = gs{2}.mirror(a,b).reverse();

            ti = grid.Ti(new_gs{:});
        end

        function ti = rotate(obj,a,rad)
            gs = obj.gs;

            for i = 1:length(gs)
                new_gs{i} = gs{i}.rotate(a,rad);
            end

            ti = grid.Ti(new_gs{:});
        end

        function ti = rotate_edges(obj,n);
            new_gs = cell(1,4);
            for i = 0:3
                new_i = mod(i - n,4);
                new_gs{new_i+1} = obj.gs{i+1};
            end
            ti = grid.Ti(new_gs{:});
        end
    end

    methods(Static)
        function obj = points(p1, p2, p3, p4)
            g1 = grid.Curve.line(p1,p2);
            g2 = grid.Curve.line(p2,p3);
            g3 = grid.Curve.line(p3,p4);
            g4 = grid.Curve.line(p4,p1);

            obj = grid.Ti(g1,g2,g3,g4);
        end

        function label(varargin)
            if nargin == 2 && ischar(varargin{2})
                label_impl(varargin{:});
            else
                for i = 1:length(varargin)
                    label_impl(varargin{i},inputname(i));
                end
            end


            function label_impl(ti,str)
                S = ti.S;

                pc = S(0.5,0.5);

                margin = 0.1;
                pw = S(  margin,      0.5);
                pe = S(1-margin,      0.5);
                ps = S(     0.5,   margin);
                pn = S(     0.5, 1-margin);


                ti.show(2,2);
                grid.place_label(pc,str);
                grid.place_label(pw,'w');
                grid.place_label(pe,'e');
                grid.place_label(ps,'s');
                grid.place_label(pn,'n');
            end
        end
    end
end