function [F,theta,phi]=ishtRectGridFFT(w,tv,pv,useProgressBar,doReal)
%ishtRectGridFFT

if nargin<4
	useProgressBar=0;
end
if nargin<5
	doReal=0
end
pv=pv(:)'; % needs to be a row vector for later outer product
% pad w with zeros if its length is not a square
L_max=ceil(sqrt(length(w)))-1;
N_tot=(L_max+1)^2;
if length(w)<N_tot % pad with zeros to make it a squared length
	w(N_tot)=0;
end

%% Create the output mesh (theta,phi) and allocate and zero the output (F)
[theta,phi]=ndgrid(tv,pv);
F=zeros(size(theta));

if useProgressBar~=0; progressbar('ishtRectGridFFT'); end

if 0
	tic
	for m=-L_max:L_max
		Gm=zeros(size(theta));
		pmo=1;
		if m>0; pmo=(-1)^m; end;
		for l=abs(m):L_max
			Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf
			n=l*(l+1)+m; % (7.39)
			wlm=w(n+1);
			Sl=legendre(l,cos(tv),'sch');
			Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
			Slm=Sl(abs(m)+1,:)'; % pull out S_l^|m| (evaluated on theta vector tt)
			Yl0=pmo*Ql*kron(ones(size(pv)),Slm);
			Gm=Gm+wlm*Yl0;
		end
		F=F+Gm.*exp(1j*m*phi);
		if useProgressBar~=0; progressbar((m+L_max+1)/(2*L_max+1)); end;
	end
	toc
	return
end

%% Perform the inverse spherical harmonic transform
tic
%expphi=exp(1j*phi);
for l=0:L_max
	Sl=legendre(l,cos(tv),'sch');
	Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf

	% m=0
	n=l*(l+1);
	wlm=w(n+1);
	Slm=sqrt(2)*Sl(1,:)'; % m=0 adjustment, see note.tex/pdf
	Ylm=Ql*Slm*ones(size(pv));
	F=F+wlm*Ylm; % m=0 contribution

	if doReal~=0
		for m=1:l
			n=l*(l+1)+m; % (7.39)
			wlm=w(n+1);
			Slm=Sl(m+1,:)';
			Ylm=(-1)^m*Ql*Slm*exp(1j*m*pv);
			F=F+2*real(wlm*Ylm);
		end
	else
		for m=1:l
			n=l*(l+1)+m;  n1=l*(l+1)-m; % (7.39)
			wlm=w(n+1);  wlm1=w(n1+1); % the weight for degree l and order m
			%if wlm==0 && wlm1==0; continue; end % skip if zero
			Slm=Sl(m+1,:)'; % pull out S_l^m (evaluated on theta vector tv)
			%Ylm=(-1)^m*Ql*kron(ones(size(pv)),Slm).*expphi.^m;
			%Ylm=(-1)^m*Ql*kron(ones(size(pv)),Slm).*exp(1j*m*phi);
			Ylm=(-1)^m*Ql*kron(ones(size(pv)),Slm).*kron(ones(size(tv')),exp(1j*m*pv));
			%Ylm=(-1)^m*Ql*Slm*exp(1j*m*pv); % outer product; pv is a row vector
			F=F+wlm*Ylm+wlm1*(-1)^m*conj(Ylm);
		end
	end
	if useProgressBar~=0; progressbar((l+1)^2/N_tot); end;
end
toc
