function [F,theta,phi]=ishtGrid2(w,tv,pv,useProgressBar)
%ishtGrid2

addpath code/

if nargin<4
	useProgressBar=false;
end

% pad w with zeros if its length is not a square
L_max=ceil(sqrt(length(w)))-1;
N_tot=(L_max+1)^2;
if length(w)<N_tot % pad with zeros to make it a squared length
	w(N_tot)=0;
end

%% Create the output mesh (theta,phi) and allocate and zero the output (F)
[theta,phi]=ndgrid(tv,pv);
F=zeros(size(theta));

if useProgressBar~=0; progressbar('ishtGrid2'); end

%% Perform the inverse spherical harmonic transform
for l=0:L_max
	Sl=legendre(l,cos(tv),'sch');
	Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
	Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf
	for m=0:l
		n=l*(l+1)+m;  n1=l*(l+1)-m; % (7.39)
		wlm=w(n+1);  wlm1=w(n1+1); % the weight for degree l and order m
		Slm=Sl(m+1,:)'; % pull out S_l^m (evaluated on theta vector tt)
		Ylm=(-1)^m*Ql*kron(ones(size(pv)),Slm).*exp(1j*m*phi);
		F=F+wlm*Ylm; % m contribution
		if m~=0 % add the -m contribution
			F=F+wlm1*(-1)^m*conj(Ylm);
		end
	end
	if useProgressBar~=0; progressbar((l+1)^2/N_tot); end;
end
