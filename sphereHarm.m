function [Ylm,theta,phi]=sphereHarm(l,m,tt,pp)
%sphereHarm spherical harmonic of degree l and order m
% see notes.tex/pdf for maths; seems robust for l to 2000
	if isvector(tt) && isvector(pp) && abs(m)<=l && l>=0
		tt=tt(:)'; pp=pp(:)'; % ensure they are row vectors
		[theta,phi]=ndgrid(tt,pp);
		Sl=legendre(l,cos(tt),'sch');
		Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment
		Ql=sqrt((2*l+1)/(8*pi));
		Slm=Sl(abs(m)+1,:)'; % pull out S_l^|m|
		Ylm=(-1)^m*Ql*kron(ones(size(pp)),Slm).*exp(1j*abs(m)*phi);
		if m<0
			Ylm=(-1)^m*conj(Ylm);
		end
	else
		error('@@ invalid inputs or out-of-range')
	end
end

% to do: for efficiency a more general version
% function [Yl,theta,phi]=sphereHarmBank(l,tt,pp)
% where Yl is multidimensional
% Yl= [Yl(0,theta,phi) ... Yl(m,theta,phi) ... Yl(l,theta,phi)]
% somewhat mimicking the legendre call
