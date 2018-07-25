function [IC2F,IF2C] = intOpMC_orders_2to2_ratio2to1(mc)

mf = 2*(mc-1) + 1;

stencil_F2C = [1.0./4.0,1.0./2.0,1.0./4.0];
stencil_width = length(stencil_F2C);
stencil_halfwidth = (stencil_width-1)/2;

BC_F2C = [1.0./2.0,1.0./2.0];

Hc = speye(mc,mc);
Hc(1,1) = 1/2;
Hc(end,end) = 1/2;

Hf = 1/2*speye(mf,mf);
Hf(1,1) = 1/4;
Hf(end,end) = 1/4; 

IF2C = sparse(mc,mf);
[BCrows, BCcols] = size(BC_F2C);
IF2C(1:BCrows, 1:BCcols) = BC_F2C;
IF2C(mc-BCrows+1:mc, mf-BCcols+1:mf) = rot90(BC_F2C,2);

for i = BCrows+1 : mc-BCrows
	IF2C(i,(2*i-stencil_halfwidth-1) :(2*i+stencil_halfwidth-1))...
		 = stencil_F2C;
end

IC2F = Hf\(IF2C'*Hc);






