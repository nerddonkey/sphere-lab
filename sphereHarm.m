function [Ylm,theta,phi]=sphereHarm(l,m,tt,pp)
%sphereHarm spherical harmonic of degree l and order m
% see notes.tex/pdf for maths; seems robust for l to 2000

if isvector(tt) && isvector(pp) && abs(m)<=l && l>=0
	[theta,phi]=ndgrid(tt,pp);
	Sl=legendre(l,cos(tt),'sch');
	Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
	Ql=sqrt((2*l+1)/(8*pi));
	Slm=Sl(abs(m)+1,:)'; % pull out S_l^|m| row and transpose to column
	Ylm=(-1)^m*Ql*kron(ones(1,length(pp)),Slm).*exp(1j*abs(m)*phi);
	if m<0
		Ylm=(-1)^m*conj(Ylm);
	end
else % one day can add code to deal with mesh data
	error('@@ invalid inputs or out-of-range')
end
