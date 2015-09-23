function [Ylm,theta,phi]=sphHarmBankmGrid(l,m,tv,pv)
%sphHarmBankmGrid spherical harmonic of degree l and order m

[YB,theta,phi]=sphHarmBankGrid(l,tv,pv);
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
