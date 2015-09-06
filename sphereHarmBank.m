function [Yl,theta,phi]=sphereHarmBank(l,tt,pp)
%sphereHarmBank spherical harmonic of degree l and all non-negative orders m
% see notes.tex/pdf for maths; seems robust for l to 2000

if isvector(tt) && isvector(pp) && l>=0
	[theta,phi]=ndgrid(tt,pp);
	Sl=legendre(l,cos(tt),'sch');
	Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
	Ql=sqrt((2*l+1)/(8*pi));
	for m=0:l
		Slm=Sl(m+1,:)'; % pull out S_l^m row and transpose to column
		Yl{m+1}=(-1)^m*Ql*kron(ones(1,length(pp)),Slm).*exp(1j*m*phi); % cell array
	end
else
	error('@@ invalid inputs or out-of-range')
end

function [Ylm,theta,phi]=sphereHarmBankm(l,m,tt,pp)
%sphereHarmBankm spherical harmonic of degree l and order m via sphereHarmBank

[YB,theta,phi]=sphereHarmBank(l,tt,pp);
Ylm=YB{abs(m)+1};
if m<0
   Ylm=(-1)^m*conj(Ylm);
end
