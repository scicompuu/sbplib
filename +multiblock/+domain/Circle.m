classdef Circle < multiblock.DefCurvilinear
    properties
        r, c

        hs
        r_arc
        omega
    end

    methods
        function obj = Circle(r, c, hs)
            default_arg('r', 1);
            default_arg('c', [0; 0]);
            default_arg('hs', 0.435);


            % alpha = 0.75;
            % hs = alpha*r/sqrt(2);

            % Square should not be a square, it should be an arc. The arc radius
            % is chosen so that the three angles of the meshes are all equal.
            % This gives that the (half)arc opening angle of should be omega = pi/12
            omega = pi/12;
            r_arc = hs*(2*sqrt(2))/(sqrt(3)-1); %  = hs* 1/sin(omega)
            c_arc = c - [(1/(2-sqrt(3))-1)*hs; 0];

            cir = parametrization.Curve.circle(c,r,[-pi/4 pi/4]);

            c2 = cir(0);
            c3 = cir(1);

            s1 = [-hs; -hs];
            s2 = [ hs; -hs];
            s3 = [ hs;  hs];
            s4 = [-hs;  hs];

            sp2 = parametrization.Curve.line(s2,c2);
            sp3 = parametrization.Curve.line(s3,c3);

            Se1 = parametrization.Curve.circle(c_arc,r_arc,[-omega, omega]);
            Se2 = Se1.rotate(c,pi/2);
            Se3 = Se2.rotate(c,pi/2);
            Se4 = Se3.rotate(c,pi/2);


            S = parametrization.Ti(Se1,Se2,Se3,Se4).rotate_edges(-1);

            A = parametrization.Ti(sp2, cir, sp3.reverse, Se1.reverse);
            B = A.rotate(c,1*pi/2).rotate_edges(-1);
            C = A.rotate(c,2*pi/2).rotate_edges(-1);
            D = A.rotate(c,3*pi/2).rotate_edges(0);

            blocks = {S,A,B,C,D};
            blocksNames = {'S','A','B','C','D'};

            conn = cell(5,5);
            conn{1,2} = {'e','w'};
            conn{1,3} = {'n','s'};
            conn{1,4} = {'w','s'};
            conn{1,5} = {'s','w'};

            conn{2,3} = {'n','e'};
            conn{3,4} = {'w','e'};
            conn{4,5} = {'w','s'};
            conn{5,2} = {'n','s'};

            boundaryGroups = struct();
            boundaryGroups.E = multiblock.BoundaryGroup({2,'e'});
            boundaryGroups.N = multiblock.BoundaryGroup({3,'n'});
            boundaryGroups.W = multiblock.BoundaryGroup({4,'n'});
            boundaryGroups.S = multiblock.BoundaryGroup({5,'e'});
            boundaryGroups.all = multiblock.BoundaryGroup({{2,'e'},{3,'n'},{4,'n'},{5,'e'}});

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blocksNames);

            obj.r     = r;
            obj.c     = c;
            obj.hs    = hs;
            obj.r_arc = r_arc;
            obj.omega = omega;
        end

        function ms = getGridSizes(obj, m)
            m_S = m;

            % m_Radial
            s = 2*obj.hs;
            innerArc = obj.r_arc*obj.omega;
            outerArc = obj.r*pi/2;
            shortSpoke = obj.r-s/sqrt(2);
            x = (1/(2-sqrt(3))-1)*obj.hs;
            longSpoke =  (obj.r+x)-obj.r_arc;
            m_R = parametrization.equal_step_size((innerArc+outerArc)/2, m_S, (shortSpoke+longSpoke)/2);

            ms = {[m_S m_S], [m_R m_S], [m_S m_R], [m_S m_R], [m_R m_S]};
        end
    end
end
