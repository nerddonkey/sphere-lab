function [Ylm,theta,phi]=sphereHarmBankm(l,m,tt,pp)
%sphereHarmBankm spherical harmonic of degree l and order m
% pulls out the component from sphereHarmBank
% same interface as function [Ylm,theta,phi]=sphereHarm(l,m,tt,pp)

[YB,theta,phi]=sphereHarmBank(l,tt,pp);
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
