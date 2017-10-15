function [stencil_F2C,BC_F2C,HcU,HfU] = intOpAWW_orders_6to6_ratio_2to1_accC2F4_accF2C3
%INT_ORDERS_6TO6_RATIO_2TO1_ACCC2F4_ACCF2C3_STENCIL_13_BC_4_17
%    [STENCIL_F2C,BC_F2C,HCU,HFU] = INT_ORDERS_6TO6_RATIO_2TO1_ACCC2F4_ACCF2C3_STENCIL_13_BC_4_17

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    11-Jul-2017 14:49:13

stencil_F2C = [1.074218750000866e-3,3.0./5.12e2,-6.445312500005198e-3,-2.5e1./5.12e2,1.6113281250013e-2,7.5e1./2.56e2,4.785156249999827e-1,7.5e1./2.56e2,1.6113281250013e-2,-2.5e1./5.12e2,-6.445312500005198e-3,3.0./5.12e2,1.074218750000866e-3];
if nargout > 1
    BC_F2C = reshape([1.0./2.0,1.257360194105895e-16,-6.179984855389442e-17,3.122529826246464e-18,6.876076086160102e-1,1.5e1./3.2e1,-3.461879841386957e-1,3.502577439820879e-2,3.099721991997423e-3,2.228547003766734e-1,9.36365299936193e-3,-3.157910466040993e-3,-1.057891432064606e-1,2.355631659452257e-1,6.070385985337517e-1,-4.847496617839159e-2,-4.809232246873087e-3,5.154694068408974e-3,7.049223649022339e-1,1.01674934980955e-2,3.460117842882258e-2,-4.99662149858484e-2,5.210650763094795e-1,2.233592875303227e-1,-2.239024779745777e-3,4.254668119952138e-4,9.213872590841121e-3,3.910524351558031e-1,-9.048711215107312e-2,8.597375629526316e-2,-3.468219752858717e-1,3.250682992395968e-1,-1.309107466664176e-1,1.199249722305999e-1,-4.071552384959323e-1,1.474417123938745e-1,-3.028638773902806e-2,3.046016554982085e-2,-1.081315128181481e-1,1.120705238850529e-2,7.977722870035704e-2,-6.772545887787698e-2,2.002704188375896e-1,-5.427183680490626e-2,6.524845570188395e-2,-5.637610037043228e-2,1.711235764016971e-1,-4.20364690240717e-2,4.639548341457899e-3,-4.055884445499834e-3,1.258630924936465e-2,-2.776451559295031e-3,-1.023784366070731e-2,8.817108288936754e-3,-2.659664906860907e-2,7.14445034054856e-3,-1.435466025443095e-3,1.240274208150634e-3,-3.764718680840621e-3,1.02871629501858e-3,1.032012418491475e-3,-8.794183904932091e-4,2.597922814458923e-3,-6.571159498039887e-4,1.892022767233693e-4,-1.612267049239958e-4,4.762858493182703e-4,-1.204712574642361e-4],[4,17]);
end
if nargout > 2
    t2 = [3.159490740740741e-1,1.390393518518519,6.275462962962963e-1,1.240509259259259,9.116898148148148e-1,1.013912037037037];
    HcU = t2;
end
if nargout > 3
    HfU = t2;
end
