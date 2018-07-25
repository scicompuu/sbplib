function [stencil_F2C,BC_F2C,HcU,HfU] = intOpAWW_orders_6to6_ratio_2to1_accC2F3_accF2C4
%INT_ORDERS_6TO6_RATIO_2TO1_ACCC2F3_ACCF2C4_STENCIL_13_BC_3_17
%    [STENCIL_F2C,BC_F2C,HCU,HFU] = INT_ORDERS_6TO6_RATIO_2TO1_ACCC2F3_ACCF2C4_STENCIL_13_BC_3_17

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    21-May-2018 15:35:47

stencil_F2C = [-5.95703125e-3,3.0./5.12e2,3.57421875e-2,-2.5e1./5.12e2,-8.935546875e-2,7.5e1./2.56e2,3.17e2./5.12e2,7.5e1./2.56e2,-8.935546875e-2,-2.5e1./5.12e2,3.57421875e-2,3.0./5.12e2,-5.95703125e-3];
if nargout > 1
    BC_F2C = reshape([5.233890618131365e-1,-1.594459192645467e-2,3.532688727637403e-2,8.021234957689208e-1,3.90683205173562e-1,-1.73222951632239e-1,-8.87662686483442e-2,2.90091235796637e-1,-1.600356115709148e-1,-1.044025375027475e-1,2.346179009198368e-1,6.091329306528956e-1,1.561275522703128e-1,-1.168382445709856e-1,1.040987887887311,-2.387061036980731e-1,1.363504965974361e-1,1.082611928255256e-1,-5.745310654054326e-1,3.977694249198785e-1,-9.376217911402619e-1,3.554518646054656e-2,5.609157787396987e-5,-1.564625311232018e-1,6.107907974027401e-1,-3.786608696698368e-1,7.02265869951125e-1,2.054294270642538e-1,-1.302300112378257e-1,2.478941407690889e-1,-2.657085326191479e-1,1.568761445933572e-1,-2.632906518005349e-1,-1.228047556139644e-1,7.182193248980271e-2,-1.1291238242346e-1,5.811258780405158e-2,-3.466364400805378e-2,5.683680203338252e-2,1.781337575097077e-2,-1.030572999042704e-2,1.577197067502767e-2,-1.443802115768554e-2,8.391316308374261e-3,-1.29534117369052e-2,-1.548018627738296e-3,8.794183904936319e-4,-1.298961407229805e-3,1.573818938200601e-3,-8.940753636685258e-4,1.320610764016968e-3],[3,17]);
end
if nargout > 2
    t2 = [3.159490740740741e-1,1.390393518518519,6.275462962962963e-1,1.240509259259259,9.116898148148148e-1,1.013912037037037];
    HcU = t2;
end
if nargout > 3
    HfU = t2;
end
