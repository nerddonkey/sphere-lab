function [D]=SlepianDH(L_max,tv,pv,R_mask,useProgressBar)
%SlepianDH generate the Hermitian Slepian D matrix for mesh region
%  \int_R I_R(\unit{x}) Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})
%  L_max is the highest non-zero degree
%  tv and pv define the integration region R
%  R_mask irregular region mask for I_R(\unit{x}) with 0-1 matrix on tv-pv grid

addpath code/ % progressbar

%% Defaults
if nargin<5
	useProgressBar=0;
end
if nargin < 4
	R_mask=0; % no mask; trapSphereMaskedR will bypass applying mask
end

%% Pre-allocate D
N_tot=(L_max+1)^2;
D=zeros(N_tot,N_tot);

%% Number of elements to compute for Hermitian
total=N_tot*(N_tot+1)/2;
if useProgressBar~=0; count=0; progressbar('SlepianDH overall','inner loop'); end;

%% Populate the Hermitian D using (8.27)
for r=1:N_tot % rows of D
	if useProgressBar~=0 && count~=0; progressbar([],0); end;
	l=floor(sqrt(r-1)); m=r-1-l*(l+1); % (7.40) here n=r-1
	f=sphHarmGrid(l,m,tv,pv);
	D(r,r)=trapSphereMaskedR(f.*conj(f),tv,pv,R_mask); % diagonal element
	if useProgressBar~=0; count=count+1; progressbar(count/total,[]); end
	for c=r+1:N_tot % columns of D
		p=floor(sqrt(c-1)); q=c-1-p*(p+1); % (7.40) here n=c-1
		g=sphHarmGrid(p,q,tv,pv);
		D(r,c)=trapSphereMaskedR(g.*conj(f),tv,pv,R_mask); % (8.27)
		D(c,r)=conj(D(r,c)); % Hermitian symmetric element
		if useProgressBar~=0; count=count+1; progressbar(count/total,(c-r)/(N_tot-r)); end
	end
	if useProgressBar~=0; progressbar(count/total,[]); end;
end
