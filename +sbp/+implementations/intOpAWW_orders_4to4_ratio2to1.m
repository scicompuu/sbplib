function [IC2F,IF2C,Hc,Hf] = IntOp_orders_4to4_ratio2to1(mc,hc,ACC)

% ACC is a string.
% ACC = 'C2F' creates IC2F with one order of accuracy higher than IF2C.
% ACC = 'F2C' creates IF2C with one order of accuracy higher than IC2F.
ratio = 2;
mf = ratio*mc-1; 
hf = hc/ratio;

switch ACC
    case 'F2C'
        [stencil_F2C,BC_F2C,HcU,HfU] = ...
        sbp.implementations.intOpAWW_orders_4to4_ratio_2to1_accC2F2_accF2C3;
    case 'C2F'
        [stencil_F2C,BC_F2C,HcU,HfU] = ...
        sbp.implementations.intOpAWW_orders_4to4_ratio_2to1_accC2F3_accF2C2;
end

stencil_width = length(stencil_F2C);
stencil_hw = (stencil_width-1)/2;
[BC_rows,BC_cols] = size(BC_F2C);

%%% Norm matrices %%%
Hc = speye(mc,mc);
HcUm = length(HcU); 
Hc(1:HcUm,1:HcUm) = spdiags(HcU',0,HcUm,HcUm);
Hc(mc-HcUm+1:mc,mc-HcUm+1:mc) = spdiags(rot90(HcU',2),0,HcUm,HcUm);
Hc = Hc*hc;

Hf = speye(mf,mf);
HfUm = length(HfU);
Hf(1:HfUm,1:HfUm) = spdiags(HfU',0,HfUm,HfUm);
Hf(mf-length(HfU)+1:mf,mf-length(HfU)+1:mf) = spdiags(rot90(HfU',2),0,HfUm,HfUm);
Hf = Hf*hf;
%%%%%%%%%%%%%%%%%%%%%%

%%% Create IF2C from stencil and BC
IF2C = sparse(mc,mf);
for i = BC_rows+1 : mc-BC_rows
    IF2C(i,ratio*i-1+(-stencil_hw:stencil_hw)) = stencil_F2C; %#ok<SPRIX>
end
IF2C(1:BC_rows,1:BC_cols) = BC_F2C;
IF2C(end-BC_rows+1:end,end-BC_cols+1:end) = rot90(BC_F2C,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create IC2F using symmetry condition %%%%
IC2F = Hf\IF2C.'*Hc;