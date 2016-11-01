classdef Ti3D
    properties
        gs % {6}Surfaces
        V  % FunctionHandle(XI,ETA,ZETA)
    end
    
    methods
        % TODO function to label boundary names.
        %  function to find largest and smallest delta h in the grid. Maybe shouldnt live here
        function obj = Ti3D(CW,CE,CS,CN,CB,CT)
            obj.gs = {CE,CW,CS,CN,CB,CT};
            
            gw = CW.g;
            ge = CE.g;
            gs = CS.g;
            gn = CN.g;
            gb = CB.g;
            gt = CT.g;
            
            function o = V_fun(XI,ETA,ZETA)
                XI=XI';
                ETA=ETA';
                ZETA=ZETA';
                
                one=0*ETA+1;
                zero=0*ETA;
                
                Sw = gw(ETA,(1-ZETA));
                Se = ge((1-ETA),(1-ZETA));
                Ss = gs(XI,ZETA);
                Sn = gn((1-XI),(1-ZETA));
                Sb = gb((1-XI),ETA);
                St = gt(XI,ETA);
                
                Ewt = gw(ETA,zero);
                Ewb = gw(ETA,one);               
                Ews = gw(zero,1-ZETA);
                Ewn = gw(one,1-ZETA);
                Eet = ge(1-ETA,zero);
                Eeb = ge(1-ETA,one);
                Ees = ge(one,1-ZETA);
                Een = ge(zero,1-ZETA);
                Enb = gn(1-XI,one);
                Ent = gn(1-XI,zero);
                Est = gs(XI,one);
                Esb = gs(XI,zero);
                
                Cwbs = gw(zero,one);
                Cwbn = gw(one,one);
                Cwts = gw(zero,zero);
                Cwtn = gw(one,zero);
                Cebs = ge(one,one);
                Cebn = ge(zero,one);
                Cets = ge(one,zero);
                Cetn = ge(zero,zero);
                
                
                X1 = (1-XI).*Sw(1,:,:) + XI.*Se(1,:,:);
                X2 = (1-ETA).*Ss(1,:,:) + ETA.*Sn(1,:,:);
                X3 = (1-ZETA).*Sb(1,:,:) + ZETA.*St(1,:,:);
                
                X12 = (1-XI).*(1-ETA).*Ews(1,:,:) + (1-XI).*ETA.*Ewn(1,:,:) + XI.*(1-ETA).*Ees(1,:,:) + XI.*ETA.*Een(1,:,:);
                X13 = (1-XI).*(1-ZETA).*Ewb(1,:,:) + (1-XI).*ZETA.*Ewt(1,:,:) + XI.*(1-ZETA).*Eeb(1,:,:) + XI.*ZETA.*Eet(1,:,:);
                X23 = (1-ETA).*(1-ZETA).*Esb(1,:,:) + (1-ETA).*ZETA.*Est(1,:,:) + ETA.*(1-ZETA).*Enb(1,:,:) + ETA.*ZETA.*Ent(1,:,:);
                
                X123 = (1-XI).*(1-ETA).*(1-ZETA).*Cwbs(1,:,:) + (1-XI).*(1-ETA).*ZETA.*Cwts(1,:,:) + (1-XI).*ETA.*(1-ZETA).*Cwbn(1,:,:) + ...
                    (1-XI).*ETA.*ZETA.*Cwtn(1,:,:) + XI.*(1-ETA).*(1-ZETA).*Cebs(1,:,:) + XI.*(1-ETA).*ZETA.*Cets(1,:,:) + ...
                    XI.*ETA.*(1-ZETA).*Cebn(1,:,:) + XI.*ETA.*ZETA.*Cetn(1,:,:);
                
                X = X1 + X2 + X3 - X12 - X13 - X23 + X123;
                
                
                Y1 = (1-XI).*Sw(2,:,:) + XI.*Se(2,:,:);
                Y2 = (1-ETA).*Ss(2,:,:) + ETA.*Sn(2,:,:);
                Y3 = (1-ZETA).*Sb(2,:,:) + ZETA.*St(2,:,:);
                
                Y12 = (1-XI).*(1-ETA).*Ews(2,:,:) + (1-XI).*ETA.*Ewn(2,:,:) + XI.*(1-ETA).*Ees(2,:,:) + XI.*ETA.*Een(2,:,:);
                Y13 = (1-XI).*(1-ZETA).*Ewb(2,:,:) + (1-XI).*ZETA.*Ewt(2,:,:) + XI.*(1-ZETA).*Eeb(2,:,:) + XI.*ZETA.*Eet(2,:,:);
                Y23 = (1-ETA).*(1-ZETA).*Esb(2,:,:) + (1-ETA).*ZETA.*Est(2,:,:) + ETA.*(1-ZETA).*Enb(2,:,:) + ETA.*ZETA.*Ent(2,:,:);
                
                Y123 = (1-XI).*(1-ETA).*(1-ZETA).*Cwbs(2,:,:) + (1-XI).*(1-ETA).*ZETA.*Cwts(2,:,:) + (1-XI).*ETA.*(1-ZETA).*Cwbn(2,:,:) + ...
                    (1-XI).*ETA.*ZETA.*Cwtn(2,:,:) + XI.*(1-ETA).*(1-ZETA).*Cebs(2,:,:) + XI.*(1-ETA).*ZETA.*Cets(2,:,:) + ...
                    XI.*ETA.*(1-ZETA).*Cebn(2,:,:) + XI.*ETA.*ZETA.*Cetn(2,:,:);
                
                Y = Y1 + Y2 + Y3 - Y12 - Y13 - Y23 + Y123;
                
                
                Z1 = (1-XI).*Sw(3,:,:) + XI.*Se(3,:,:);
                Z2 = (1-ETA).*Ss(3,:,:) + ETA.*Sn(3,:,:);
                Z3 = (1-ZETA).*Sb(3,:,:) + ZETA.*St(3,:,:);
                
                Z12 = (1-XI).*(1-ETA).*Ews(3,:,:) + (1-XI).*ETA.*Ewn(3,:,:) + XI.*(1-ETA).*Ees(3,:,:) + XI.*ETA.*Een(3,:,:);
                Z13 = (1-XI).*(1-ZETA).*Ewb(3,:,:) + (1-XI).*ZETA.*Ewt(3,:,:) + XI.*(1-ZETA).*Eeb(3,:,:) + XI.*ZETA.*Eet(3,:,:);
                Z23 = (1-ETA).*(1-ZETA).*Esb(3,:,:) + (1-ETA).*ZETA.*Est(3,:,:) + ETA.*(1-ZETA).*Enb(3,:,:) + ETA.*ZETA.*Ent(3,:,:);
                
                Z123 = (1-XI).*(1-ETA).*(1-ZETA).*Cwbs(3,:,:) + (1-XI).*(1-ETA).*ZETA.*Cwts(3,:,:) + (1-XI).*ETA.*(1-ZETA).*Cwbn(3,:,:) + ...
                    (1-XI).*ETA.*ZETA.*Cwtn(3,:,:) + XI.*(1-ETA).*(1-ZETA).*Cebs(3,:,:) + XI.*(1-ETA).*ZETA.*Cets(3,:,:) + ...
                    XI.*ETA.*(1-ZETA).*Cebn(3,:,:) + XI.*ETA.*ZETA.*Cetn(3,:,:);
                
                Z = Z1 + Z2 + Z3 - Z12 - Z13 - Z23 + Z123;
                o = [X;Y;Z];
            end
            
            obj.V = @V_fun;
        end
        
        %Should be rewritten so that the input is xi eta zeta 
        function [X,Y,Z] = map(obj,XI,ETA,ZETA)
            
            V = obj.V;
            
            p = V(XI,ETA,ZETA);
            X = p(1,:)';
            Y = p(2,:)';
            Z = p(3,:)';
            
        end
        
        %         function h = plot(obj,nu,nv)
        %             S = obj.S;
        %
        %             default_arg('nv',nu)
        %
        %             u = linspace(0,1,nu);
        %             v = linspace(0,1,nv);
        %
        %             m = 100;
        %
        %             X = zeros(nu+nv,m);
        %             Y = zeros(nu+nv,m);
        %
        %
        %             t = linspace(0,1,m);
        %             for i = 1:nu
        %                 p = S(u(i),t);
        %                 X(i,:) = p(1,:);
        %                 Y(i,:) = p(2,:);
        %             end
        %
        %             for i = 1:nv
        %                 p = S(t,v(i));
        %                 X(i+nu,:) = p(1,:);
        %                 Y(i+nu,:) = p(2,:);
        %             end
        %
        %             h = line(X',Y');
        %         end
        %
        %
        %         function h = show(obj,nu,nv)
        %             default_arg('nv',nu)
        %             S = obj.S;
        %
        %             if(nu>2 || nv>2)
        %                 h_grid = obj.plot(nu,nv);
        %                 set(h_grid,'Color',[0 0.4470 0.7410]);
        %             end
        %
        %             h_bord = obj.plot(2,2);
        %             set(h_bord,'Color',[0.8500 0.3250 0.0980]);
        %             set(h_bord,'LineWidth',2);
        %         end
        %
        %
        %         % TRANSFORMATIONS
        %         function ti = translate(obj,a)
        %             gs = obj.gs;
        %
        %             for i = 1:length(gs)
        %                 new_gs{i} = gs{i}.translate(a);
        %             end
        %
        %             ti = grid.Ti(new_gs{:});
        %         end
        %
        %         % Mirrors the Ti so that the resulting Ti is still left handed.
        %         %  (Corrected by reversing curves and switching e and w)
        %         function ti = mirror(obj, a, b)
        %             gs = obj.gs;
        %
        %             new_gs = cell(1,4);
        %
        %             new_gs{1} = gs{1}.mirror(a,b).reverse();
        %             new_gs{3} = gs{3}.mirror(a,b).reverse();
        %             new_gs{2} = gs{4}.mirror(a,b).reverse();
        %             new_gs{4} = gs{2}.mirror(a,b).reverse();
        %
        %             ti = grid.Ti(new_gs{:});
        %         end
        %
        %         function ti = rotate(obj,a,rad)
        %             gs = obj.gs;
        %
        %             for i = 1:length(gs)
        %                 new_gs{i} = gs{i}.rotate(a,rad);
        %             end
        %
        %             ti = grid.Ti(new_gs{:});
        %         end
        %
        %         function ti = rotate_edges(obj,n);
        %             new_gs = cell(1,4);
        %             for i = 0:3
        %                 new_i = mod(i - n,4);
        %                 new_gs{new_i+1} = obj.gs{i+1};
        %             end
        %             ti = grid.Ti(new_gs{:});
        %         end
        %     end
        %
        %     methods(Static)
        %         function obj = points(p1, p2, p3, p4)
        %             g1 = grid.Curve.line(p1,p2);
        %             g2 = grid.Curve.line(p2,p3);
        %             g3 = grid.Curve.line(p3,p4);
        %             g4 = grid.Curve.line(p4,p1);
        %
        %             obj = grid.Ti(g1,g2,g3,g4);
        %         end
        %
        %         function label(varargin)
        %             if nargin == 2 && ischar(varargin{2})
        %                 label_impl(varargin{:});
        %             else
        %                 for i = 1:length(varargin)
        %                     label_impl(varargin{i},inputname(i));
        %                 end
        %             end
        %
        %
        %             function label_impl(ti,str)
        %                 S = ti.S;
        %
        %                 pc = S(0.5,0.5);
        %
        %                 margin = 0.1;
        %                 pw = S(  margin,      0.5);
        %                 pe = S(1-margin,      0.5);
        %                 ps = S(     0.5,   margin);
        %                 pn = S(     0.5, 1-margin);
        %
        %
        %                 ti.show(2,2);
        %                 grid.place_label(pc,str);
        %                 grid.place_label(pw,'w');
        %                 grid.place_label(pe,'e');
        %                 grid.place_label(ps,'s');
        %                 grid.place_label(pn,'n');
        %             end
 %                end
    end
end