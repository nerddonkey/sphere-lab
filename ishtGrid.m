function [F,theta,phi]=ishtGrid(w,tv,pv,useProgressBar)
%ishtGrid

if nargin<4
	useProgressBar=0;
end

%% pad w with zeros if its length is not a square
L_max=ceil(sqrt(length(w)))-1;
N_tot=(L_max+1)^2;
if length(w)<N_tot % pad with zeros to make it a squared length
	w(N_tot)=0;
end

if useProgressBar~=0; progressbar('ishtGrid'); end

for l=0:L_max
	[YB,theta,phi]=sphHarmBankGrid(l,tv,pv);
	if l==0
		F=zeros(size(theta)); % mesh answer to be accumulated
	end
	for m=0:l
		n=l*(l+1)+m;  n1=l*(l+1)-m; % (7.39)
		wlm=w(n+1);  wlm1=w(n1+1); % the weight for degree l and order m
		Ylm=YB{m+1};
		F=F+wlm*Ylm; % m contribution
		if m~=0 % add the -m contribution
			F=F+wlm1*(-1)^m*conj(Ylm);
		end
	end
	if useProgressBar~=0; progressbar((l+1)^2/N_tot); end;
end
