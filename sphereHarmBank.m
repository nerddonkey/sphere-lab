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
else % one day can add code to deal with mesh data
	error('@@ invalid inputs or out-of-range')
end
