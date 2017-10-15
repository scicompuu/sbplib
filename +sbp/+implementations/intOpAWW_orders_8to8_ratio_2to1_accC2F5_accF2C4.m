function [stencil_F2C,BC_F2C,HcU,HfU] = intOpAWW_orders_8to8_ratio_2to1_accC2F5_accF2C4
%INT_ORDERS_8TO8_RATIO_2TO1_ACCC2F5_ACCF2C4_STENCIL_17_BC_6_24
%    [STENCIL_F2C,BC_F2C,HCU,HFU] = INT_ORDERS_8TO8_RATIO_2TO1_ACCC2F5_ACCF2C4_STENCIL_17_BC_6_24

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    11-Jul-2017 14:44:32

stencil_F2C = [-2.564929780864058e-4,-1.220703125e-3,2.051943824691247e-3,1.1962890625e-2,-7.181803386419364e-3,-5.9814453125e-2,1.436360677283873e-2,2.99072265625e-1,4.820454915339516e-1,2.99072265625e-1,1.436360677283873e-2,-5.9814453125e-2,-7.181803386419364e-3,1.1962890625e-2,2.051943824691247e-3,-1.220703125e-3,-2.564929780864058e-4];
if nargout > 1
    BC_F2C = reshape([3.84889400178856e-1,1.112426550601988e-1,-1.318495369488137,1.887814023562389e-1,-4.11247635928463e-1,2.655099795914844e-2,1.028610922995518,2.364211968142455e-1,2.059181912003073,-4.340414328006023e-1,1.075497633269875,-7.409800036086324e-2,7.37838637890288e-2,1.306633125422226e-2,8.451322719188257e-1,-1.210055485680529e-1,2.636024797025984e-1,-1.701872129868639e-2,-4.635230469768725e-1,6.090754537930958e-1,-1.489759296753512,4.8674037213944e-1,-1.179464341635484,7.944500441079773e-2,-7.260695877364954e-2,7.032243845450474e-2,-3.47295282501933e-2,1.20392464519063e-1,-2.651360013733321e-1,1.767343065829761e-2,-2.096982321858537e-2,4.637607722376962e-3,9.018016569036531e-1,2.906245591263948e-1,-3.442417302951498e-1,2.096340897620274e-2,-2.874137204380543e-3,2.384228040892216e-3,-2.026040077382159e-2,2.573962718116355e-1,1.054557409115141e-2,-3.708993639747955e-3,1.992191138619206e-2,-1.470164138115057e-2,2.155634694638224e-2,1.340780622608816e-1,7.925832983348031e-1,-4.875561424378325e-2,-2.350572186645021e-4,3.076716550868981e-4,-8.141271258942936e-3,3.404241369471586e-3,1.185477795471203,8.70934314005432e-3,-1.262976611986137e-2,8.032247391372543e-3,1.316324784899156e-3,-4.086253158542394e-2,8.015278252366783e-1,2.175869750865401e-1,4.007238901530734e-4,3.864273806524533e-4,-1.508600083692204e-2,3.153609486577761e-3,3.042286208824603e-3,3.82047294861268e-1,4.397753338251498e-2,-2.563454837101723e-2,1.41743994040605e-3,7.608105232124472e-2,-5.977718736752081e-1,3.20050013977378e-1,6.802094242804826e-2,-4.194020001568932e-2,7.25315420916367e-2,9.186386275429062e-2,-6.573062309340446e-1,1.385772992079256e-1,1.775827558199515e-2,-1.127406881785768e-2,1.811068993199809e-2,2.791591769899442e-2,-1.928604231232766e-1,6.083089919559692e-3,-3.221756216941285e-2,1.711212226954443e-2,2.901092273409817e-2,-5.378021102520328e-2,2.992550025960758e-1,-4.666662872498935e-2,-3.465064705533973e-2,2.133105869707025e-2,-4.756361192707918e-2,-3.854942718050275e-2,2.54971564582897e-1,-3.334541871752157e-2,-1.305394187273881e-3,4.611114655140197e-4,8.244321726072152e-3,-4.278000168944293e-3,2.105663050394826e-2,-1.821301944614694e-3,5.395937659538634e-3,-2.783571877605622e-3,-8.060157292374778e-3,1.022596015716632e-2,-5.706697079870192e-2,7.93956741906313e-3,1.155407428597162e-4,2.111811854597317e-4,-7.982283350430192e-3,2.363443780394039e-3,-9.897869059620054e-3,1.174758160329132e-3,-1.714734530349539e-3,1.229820838718235e-3,-7.4133437108253e-3,-5.029932466958482e-4,6.964815291000956e-3,-1.280740277254815e-3,-2.947086158334019e-4,2.091452223113523e-4,-1.20812830748558e-3,-1.054095401829176e-4,1.280105034963581e-3,-2.33676862216909e-4,-4.642133749083303e-4,4.859542003098669e-4,-6.379322438624295e-3,1.046455654592667e-3,-2.762766278515372e-3,2.407540045947593e-4,-1.891185255825058e-4,1.906094692056257e-4,-2.389371525102409e-3,3.700689438675072e-4,-9.076857917948225e-4,7.171015964436335e-5,8.001159359856251e-4,-7.732304515235331e-4,9.164656932421384e-3,-1.312190264641192e-3,2.858518569557242e-3,-1.845518711255185e-4],[6,24]);
end
if nargout > 2
    t2 = [2.948906761778786e-1,1.525720623897707,2.57452876984127e-1,1.798113701499118,4.127080577601411e-1,1.278484623015873,9.232955798059965e-1,1.009333860859158];
    HcU = t2;
end
if nargout > 3
    HfU = t2;
end
