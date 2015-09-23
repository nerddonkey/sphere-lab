function YB=sphHarmBank(l,theta,phi)
%sphHarm

if ~isequal(size(theta),size(phi))
   error('@@ sphHarmBank: theta and phi need to be the same size')
end

Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf
YB=zeros(size(theta)); % defines shape for output
S_max=numel(theta);
for s=1:S_max % walk thru theta,phi matrices as vectors
	Sl=legendre(l,cos(theta(s)),'sch'); % vector in m, scalar in theta
   for m=0:l
   	Slm=Sl(m+1,1); % pull out S_l^m row scalar
      Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
   	YB{m+1}=(-1)^m*Ql*Slm*exp(1j*m*phi(s)); % cell array
   end
end


function Ylm=sphHarmBankm(l,m,theta,phi)
%sphHarmBankm spherical harmonic of degree l and order m

YB=sphHarmBank(l,theta,phi)
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
