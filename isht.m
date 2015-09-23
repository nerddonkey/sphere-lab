function F=isht(w,theta,phi,useProgressBar)
%isht
% Y_l^m = (-1)^m NSl P_l^m(sch) e^m\phi

addpath code/

if nargin<4
	useProgressBar=0;
end

%% pad w with zeros if its length is not a square
L_max=ceil(sqrt(length(w)))-1;
N_tot=(L_max+1)^2;
if length(w)<N_tot % pad with zeros to make it a squared length
	w(N_tot)=0;
end

S_max=numel(theta);
F=zeros(size(theta)); % mesh answer to be accumulated

if useProgressBar~=0;  progressbar('isht');  end

for s=1:S_max % walk thru theta,phi matrices as vectors
	for l=0:L_max
		Sl=legendre(l,cos(theta(s)),'sch');
		Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
		Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf
		for m=0:l
			n=l*(l+1)+m;  n1=l*(l+1)-m; % (7.39)
			wlm=w(n+1);  wlm1=w(n1+1); % the weight for degree l and order m
			Slm=Sl(m+1,1); % pull out scalar S_l^m (evaluated at single theta)
			Ylm=(-1)^m*Ql*Slm.*exp(1j*m*phi(s)); % scalar
			F(s)=F(s)+wlm*Ylm;
			if m~=0
				F(s)=F(s)+wlm1*(-1)^m*conj(Ylm);
			end
		end
	end
	if useProgressBar~=0;  progressbar(s/S_max);  end;
end
